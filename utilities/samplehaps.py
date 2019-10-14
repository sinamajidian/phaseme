


"""
Credit: https://github.com/vibansal/IntegratedPhasing/blob/master/encodereads.py 

Vikas Bansal, Integrating read-based and population-based phasing for dense and accurate haplotyping of individual genomes, Bioinformatics 2019.

Edited by Sina Majidian


"""

import os,sys,subprocess,random
import argparse





#SHAPEIT="/mnt/LTR_userdata/majid001/a/shapeit2/bin/shapeit";  ## set this path 

SHAPEIT="shapeit";  ## set this path 


## first line of output file has list of variant-ids (relative to VCF file) 
def sample_haplotypes(VCF,nsamples):
	HAPfile = VCF + '.hapsamples'; 
	snptable = {}; variantid=0;

	VCF_file = open(VCF+'.vcf');
	for line in VCF_file:
		if line[0]== '#': continue;
		var = line.split();
		chrom = var[0]; position = int(var[1]); allele1 = var[3]; alleles = var[4].split(','); allele2 = alleles[0];
		genotype = var[9].split(':');
		variantid +=1;
		snptable[(chrom,position,allele1)] = [genotype[0],1,variantid]; 
	VCF_file.close();


	if nsamples > 1: F1 = open(HAPfile,'w'); # new output file

	for r in xrange(nsamples):
		seed = random.randint(1,10000000);
		outfile = VCF + '.phased'; logfile = VCF + '.shapeit.log';
		samplehap = SHAPEIT + ' -convert' + ' --seed ' + str(r) + ' --input-graph ' + VCF + '.graph' + ' --output-sample ' + outfile + ' -L ' + logfile;
		#print >>sys.stderr,samplehap
		subprocess.call(samplehap,shell=True);
		hap1 = []; hap2 = [];
		File = open(outfile + '.haps','r');
		for line in File: 
			var = line.strip().split(); h1 = var[5]; h2 = var[6]; 
			if h1 == h2: continue; ## homozoygous variant # output of shapeit all seed are consistent, if var1 is homo. it is in all seeds
			if r ==0: 
				varid = snptable[(var[0],int(var[2]),var[3])][2]
				print >>F1, varid, ## variant-id 
			hap1.append(h1); hap2.append(h2); 
		File.close();
		if r ==0: print >>F1,'\n',
		if nsamples > 1: print >>F1, ''.join(hap1); ## only write one of the two haplotypes
		if (r+1)==10: print >>sys.stderr,'sample',r+1;
		if (r+1)==50: print >>sys.stderr,'sample',r+1;
		if (r+1)%100==0: print >>sys.stderr,'sample',r+1;
	if nsamples > 1: F1.close(); 
	#subprocess.call('rm -f ' + outfile + '.* ' + logfile,shell=True);


if len(sys.argv) < 3: 
	print >>sys.stderr, "program to sample haplotypes from SHAPEIT2 HMM ";
	print >>sys.stderr, "set path to shapeit before running this program, the shapeit graph file also needs to exist ";
	print >>sys.stderr, "usage: python samplehaps.py VCFfile number_haps(integer)";
	sys.exit();
sample_haplotypes(sys.argv[1],int(sys.argv[2])); 

