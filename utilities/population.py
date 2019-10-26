
 #!/usr/bin/env python3

import subprocess
import random
    
# credit 
#https://github.com/vibansal/IntegratedPhasing/blob/master/samplehaps.py
#https://github.com/vibansal/IntegratedPhasing/blob/master/encodereads.py



SHAPEIT="/home/ssm/Downloads/phaseme/shapeit"



THRESH=0.98

        

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
    W = 10                  # number of neighbour variants to be checked 

    for i in range(num_variants):

        for j in range(i+1, min(num_variants, i+W)):

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

if __name__ == "__main__":


	num_samples = 500

	#var_pos_list, allele1_list_samples= extract_samples()

	haplotype1_samples = []
	
	for sample_i in range(num_samples):

		subprocess.call("mkdir files", shell=True)
		seed = random.randint(1,10000000)
		shapeit_sample_hap_graph= SHAPEIT+" -convert --seed "+str(seed)+" --input-graph chr22_haplotype.graph --output-sample files/hap_sample_"+str(sample_i)+" -L files/log_sample_"+str(sample_i)

		subprocess.call(shapeit_sample_hap_graph, shell=True)
	
		haplotype_sample_address= "files/hap_sample_"+str(sample_i)+".haps"
		var_pos_list, haplotype1  = read_haplotype_sample(haplotype_sample_address)

		haplotype1_samples.append(haplotype1)

    
	pairs = extract_pairs(haplotype1_samples)
	
	file_pair_address = 'pairs_500.txt' 
	report_pairs(file_pair_address, pairs, var_pos_list)




