
 #!/usr/bin/env python3

from sys import argv
import subprocess
import random
	
# credit 
#https://github.com/vibansal/IntegratedPhasing/blob/master/samplehaps.py
#https://github.com/vibansal/IntegratedPhasing/blob/master/encodereads.py





def run_shapeit_graph(chr_vcf, data_1000G_address, chrom): # run per chromosom

	subprocess.call("mkdir files_"+str(chrom), shell=True)	

	shapeit_check= SHAPEIT+" -check --input-vcf "+chr_vcf+" -R "+data_1000G_address+"1000GP_Phase3_chr"+str(chrom)+".hap.gz "+data_1000G_address+"1000GP_Phase3_chr"+str(chrom)+".legend.gz  "+data_1000G_address+"1000GP_Phase3.sample --output-log files_"+str(chrom)+"/shapeit_check >> files_"+str(chrom)+"/samples.log"

	subprocess.call(shapeit_check, shell=True)

	shapeit_generate_graph= SHAPEIT+" --input-vcf  "+chr_vcf+" -R "+data_1000G_address+"1000GP_Phase3_chr"+str(chrom)+".hap.gz "+data_1000G_address+"1000GP_Phase3_chr"+str(chrom)+".legend.gz "+data_1000G_address+"1000GP_Phase3.sample  -M "+data_1000G_address+"genetic_map_chr"+str(chrom)+"_combined_b37.txt  --output-log files_"+str(chrom)+"/shapeit_graph --output-graph "+chr_vcf[:-4]+".graph --exclude-snp  files_"+str(chrom)+"/shapeit_check.snp.strand.exclude >> files_"+str(chrom)+"/samples.log"

	subprocess.call(shapeit_generate_graph, shell=True)

	return 1



		

def read_haplotype_sample(haplotype_sample_address):
	
	sample_file = open(haplotype_sample_address,'r')
	
	var_pos_list = []				 # genomic position of variants
	haplotype1 = []				 # hetrozygous variant only
	
	for line in sample_file: 
		
		line_splitted = line.strip().split(); 
		
		chrom = line_splitted[0]
		var_pos = line_splitted[2]
		
		allele1 = line_splitted[5]
		allele2 = line_splitted[6]
		
		if allele1 != allele2:
			haplotype1.append(allele1)
			var_pos_list.append(var_pos)
			
	return var_pos_list, haplotype1

	
def pairwise(haplotype1_samples, nsamples, i, j):
	
	identical_phasing = 0
	opposite_phasing = 0
	
	for sample_i in range(nsamples):
		
		if haplotype1_samples[sample_i][i] == '0' and haplotype1_samples[sample_i][j] == '0': identical_phasing +=1
		elif haplotype1_samples[sample_i][i] == '1' and haplotype1_samples[sample_i][j] == '1': identical_phasing +=1
			
		elif haplotype1_samples[sample_i][i] == '0' and haplotype1_samples[sample_i][j] == '1': opposite_phasing +=1 
		elif haplotype1_samples[sample_i][i] == '1' and haplotype1_samples[sample_i][j] == '0': opposite_phasing +=1 
			
	return identical_phasing, opposite_phasing



def extract_pairs(haplotype1_samples):

	num_samples =  len(haplotype1_samples)
	num_variants = len(haplotype1_samples[0])

	pairs =[];

	for i in range(num_variants):

		for j in range(i+1, min(num_variants, i+NEIGHBOURS)):

			identical_phasing, opposite_phasing= pairwise(haplotype1_samples, num_samples, i, j)

			f = float(identical_phasing+0.5)/(identical_phasing+opposite_phasing+1)

			if f > THRESH: 
				pairs.append([i,j,'identical']);  

			elif 1.0-f > THRESH: 
				pairs.append([i,j,'opposite']); 
	

	pairs_sorted = sorted(pairs, key=lambda item: item[0])

	return pairs_sorted


  
	

def report_pairs(file_pair_address, pairs, var_pos_list):
	
	file_pairs= open(file_pair_address, 'w');  # output

	for pair in pairs:

		line_out = str(var_pos_list[pair[0]])+'\t'+str(var_pos_list[pair[1]])+'\t'+pair[2]+'\n'
		file_pairs.write(line_out)
	
	file_pairs.close()
	return 1





def sample_haplotype_graph(input_graph, num_samples):

	subprocess.call("mkdir files_"+str(chrom)+"/samples", shell=True)	
	for sample_i in range(num_samples):


		seed = random.randint(1,10000000)
		shapeit_sample_hap_graph= SHAPEIT+" -convert --seed "+str(seed)+" --input-graph "+str(input_graph)+" --output-sample files_"+str(chrom)+"/samples/sample_"+str(sample_i)+" -L files_"+str(chrom)+"/samples/sample_"+str(sample_i)+">> files_"+str(chrom)+"/samples.log"

		subprocess.call(shapeit_sample_hap_graph, shell=True)
		
		
	return 1


def read_haplotype_samples(num_samples):
	haplotype1_samples = []

	for sample_i in range(num_samples):
		
	
		
		haplotype_sample_address= "files_"+str(chrom)+"/samples/sample_"+str(sample_i)+".haps"		
		
		var_pos_list, haplotype1  = read_haplotype_sample(haplotype_sample_address)

		haplotype1_samples.append(haplotype1)

		if sample_i>0:
			if var_pos_list_pre!=var_pos_list:
				print("inconsistency in shapeit output, sample index",sample_i)
				exit(1)
		
		var_pos_list_pre=var_pos_list
		
				

	return haplotype1_samples, var_pos_list




if __name__ == "__main__":

	
	SHAPEIT="/home/ssm/Documents/phaseme/shapeit"


	NEIGHBOURS = 10			   # number of neighbour variants to be checked 
	THRESH=0.90       		  # the extent of  between samples 
	num_samples = 500			# number that we sample the haplotype graph (output of shapeit)


	
	# input_vcf = "out.vcf"  # argv[1]
	# data_1000G_address = "data/1000GP_Phase3/" # argv[2]
	#downlaod_1000g= "wget https://mathgen.stats.ox.ac.uk/impute/1000GP_Phase3.tgz; wget https://mathgen.stats.ox.ac.uk/impute/1000GP_Phase3_chrX.tgz "
	#subprocess.call(downlaod_1000g, shell=True)
	# run_shapeit_graph(input_vcf, data_1000G_address)

	data_1000G_address = "/home/ssm/Documents/phaseme/data1/1000g/"



	print("In case of error, please check the log files as well.")
	

	for chrom in range(1,23):
		print("Working on chr ",chrom)
	
		
		input_vcf = "ccs/ccs.vcf" 
		chr_vcf = "ccs/"+str(chrom)+"/"+str(chrom)+".vcf"		






		## grep "#" input.vcf | sed 's/Type=Integer,Description="Phase set identifier">/Type=String,Description="Phase set in which this variant falls">/' > ont/22/22.vcf
		## grep "^22\b" input.vcf | grep -v "\.:\.:\.:\." >> ont/22/22.vcf

		subprocess.call("mkdir ccs/"+str(chrom), shell=True)
		#subprocess.call("grep \"#\" "+input_vcf+" | sed \'s/Type=Integer,Description=\"Phase set identifier\"/Type=String,Description=\"Phase set in which this variant falls\"/\' >"+chr_vcf, shell=True)	

		subprocess.call("grep \"#\" "+input_vcf+">"+chr_vcf, shell=True)	
		
		subprocess.call("grep \"^"+str(chrom)+"\\b\" "+input_vcf+" | grep -v \"\\.:\\.:\\.:\\.\">>"+chr_vcf, shell=True)	
		
		
	
		run_shapeit_graph(chr_vcf, data_1000G_address, chrom)
		haplotype_graph = chr_vcf[:-4]+".graph" 
		print("haplotype graph is generated: "+haplotype_graph)


		sample_haplotype_graph(haplotype_graph, num_samples)



		print(str(num_samples)+" haplotype samples are generated from haplotype graph .")
		
		haplotype1_samples, var_pos_list = read_haplotype_samples(num_samples)
		print("samples are read.")
		

		pairs = extract_pairs(haplotype1_samples)
	
		file_pair_address = chr_vcf[:-4]+"_pairs_"+str(num_samples)+"_"+str(THRESH)+".txt"
		report_pairs(file_pair_address, pairs, var_pos_list)
		print(str(len(pairs))+" pairs are reported in "+file_pair_address)

	






