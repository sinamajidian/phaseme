
 #!/usr/bin/env python3

from sys import argv
import subprocess
import random
    
# credit 
#https://github.com/vibansal/IntegratedPhasing/blob/master/samplehaps.py
#https://github.com/vibansal/IntegratedPhasing/blob/master/encodereads.py






def run_shapeit(input_vcf, data_1000G_address):

	vcf_chr = input_vcf[:-4]



	# check for ungenotype
	# check for vcf not to be empty 
	


	chrom = 22	
	
	# for chrom in range(1,24):


	## split input vcf per chromosome


	# vcf_chr = input_vcf[:-4]+str(chrom)+".vcf"
	# split_vcf = "grep # "+input_vcf+">"+vcf_chr+"; " grep "+str(chrom)+" >>"+vcf_chr
	# subprocess.call(split_vcf, shell=True)



	shapeit_check= SHAPEIT+" -check --input-vcf "+vcf_chr+".vcf -R "+1000G_address+"1000GP_Phase3_chr"+str(chrom)+".hap.gz "+1000G_address+"1000GP_Phase3_chr"+str(chrom)+".legend.gz  "+1000G_address+"1000GP_Phase3.sample --output-log files/log_shapeit_check"

	subprocess.call(shapeit_check, shell=True)

	shapeit_generate_graph= SHAPEIT+" --input-vcf  "+vcf_chr+" -R "+1000G_address+"1000GP_Phase3_chr"+str(chrom)+".hap.gz "+1000G_address+"1000GP_Phase3_chr"+str(chrom)+".legend.gz "+1000G_address+"1000GP_Phase3.sample  -M "+1000G_address+"genetic_map_chr"+str(chrom)+"_combined_b37.txt  --output-log files/log_shapeit_graph --output-graph "+vcf_chr+".graph --exclude-snp  "+vcf_chr+".snp.strand.exclude"

	subprocess.call(shapeit_generate_graph, shell=True)

	 # return  output of subprocess


	return 1



        

def read_haplotype_sample(haplotype_sample_address):
    
    sample_file = open(haplotype_sample_address,'r')
    
    var_pos_list = []                 # genomic position of variants
    haplotype1 = []                 # hetrozygous variant only
    
    for line in sample_file: 
        
        line_splitted = line.strip().split(); 
        
        chrom =  line_splitted[0]
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

        pairs_sorted=sorted(pairs, key=lambda item: item[0])

    return pairs_sorted


  
    

def report_pairs(file_pair_address, pairs, var_pos_list):
    
    file_pairs= open(file_pair_address, 'w');  # output

    for pair in pairs:

        line_out = str(var_pos_list[pair[0]])+'\t'+str(var_pos_list[pair[1]])+'\t'+pair[2]+'\n'
        file_pairs.write(line_out)
    
    file_pairs.close()
    return 1





def sample_haplotype_graph(input_graph, num_samples):
	haplotype1_samples = []
	
	for sample_i in range(num_samples):

		subprocess.call("mkdir files", shell=True)
		seed = random.randint(1,10000000)
		shapeit_sample_hap_graph= SHAPEIT+" -convert --seed "+str(seed)+" --input-graph "+str(input_graph)+" --output-sample files/hap_sample_"+str(sample_i)+" -L files/log_sample_"+str(sample_i)

		subprocess.call(shapeit_sample_hap_graph, shell=True)
	
		haplotype_sample_address= "files/hap_sample_"+str(sample_i)+".haps"
		var_pos_list, haplotype1  = read_haplotype_sample(haplotype_sample_address)

		haplotype1_samples.append(haplotype1)
	return haplotype1_samples, var_pos_list










if __name__ == "__main__":

	
	SHAPEIT="/home/ssm/Downloads/phaseme/shapeit"

	NEIGHBOURS = 10                 # number of neighbour variants to be checked 

	THRESH=0.98

	num_samples = 500  # number that we sample the haplotype graph (output of shapeit)





	# input_vcf = "out.vcf"  # argv[1]
	
	# data_1000G_address = "data/1000GP_Phase3/" # argv[2]

	
	# run_shapeit(input_vcf, data_1000G_address)

	
	haplotype_graph = "out.graph" # "chr22.graph"

	haplotype1_samples, var_pos_list = sample_haplotype_graph(haplotype_graph, num_samples)
    
	pairs = extract_pairs(haplotype1_samples)
	
	
	file_pair_address = 'pairs_500.txt' 
	report_pairs(file_pair_address, pairs, var_pos_list)




