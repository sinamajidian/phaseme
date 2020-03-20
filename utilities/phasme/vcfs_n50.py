



from sys import argv
import numpy as np



def read_vcf_file_with_truth(vcf_file_address):    
    
    """
    Reading the vcf file
    
    input: vcf file
    outputs: 
            lines_list: list of string. each string is a line of phased vcf file.
            var_pos_het_list: genomic position of phased hetrozygous variants in the vcf file. 
            line_number_het_list: list of line numbers in the vcf file that are phased hetrozygous variant (needed in phasing)
            id_blocks: list of ids of phas blocks
            allele_blocks: list of list
            var_pos_blocks: list of list
            stats_vcf: [homozygous0_num, homozygous1_num, hetrozygous_nonphased, hetrozygous_phased, genomic_length_blocks, n50]
    
    
    """
    
    
    
    
    allele_truth_dic = {}
    
    vcf_file = open(vcf_file_address,'r')

    lines_list=[]                  #  lines  of phased vcf  needed for reporting improved VCF 
    var_pos_het_list=[]                # position of phased hetrozygous variants all blocks consequently. 

    # The followings are for phased hetrozygous variants.
    id_blocks = []                 # list of list. Outer list corresponds to phase block. Inner list contains block_id
    
    allele_blocks = []             # list of list. Outer list corresponds to phase block. Inner list contains alleles of hetro variants 
    var_pos_blocks = []            # list of list. Outer list corresponds to phase block. Inner list contains genomic positions of hetro variants   

    
    line_number_het_list = []            # line number of phased hetrozygous variant. We need it for reporting improved version
    lines_list = []
    
    first_het_variant = True
    line_number = 0
    
    homozygous0_num = 0
    homozygous1_num = 0
    hetrozygous_nonphased = 0
    hetrozygous_phased = 0
    
    
    for line in vcf_file:
        
        line_number += 1
        
        line_strip = line.strip() 
        
        lines_list.append(line_strip)
        
        if line_strip.startswith('#'):
            pass 
            #header_lines_list.append(line_strip)
            #sample_names = line_strip.split('\t')[9:11]            # last line of header contains sample name
        else:

            line_parts=line_strip.split('\t') 

            chrom = line_parts[0]

#             if str(chrom_output)!=chrom:
#                 print(chrom_output,chrom)

            var_pos = int(line_parts[1])                           # genomic position of variants
            

            
            format_genotype, values_genotype = line_parts[8:10]    # 'GT:GQ:DP:AF:GL:PS', '0|1:255:.:.:.,0,.:60780'
            
            values_genotype_splitted = values_genotype.split(':')
            format_genotype_splitted = format_genotype.split(':')
            
            gt_index = format_genotype_splitted.index("GT")           #  index of allele in  values_genotype 
            
            allele = values_genotype_splitted[gt_index]
            
            
            # how should we handle '2' in allele ?  
            
#             if './.' in allele:
#                 print("There is a vriant with genomic position "+str(var_pos)+" that is not genotyped. Remove it first.")
#                 exit(1)
                
            
            # if '/' in allele: print("There is a vriant with genomic position "+str(var_pos)+" that is not phased. Remove it first.")            

            if allele == '0/1' or allele == '1/0':
                hetrozygous_nonphased += 1

            if allele == '0|0' or allele == '0/0':
                homozygous0_num += 1
                
            if allele == '1|1' or allele == '1/1':
                homozygous1_num += 1
                    
            if (allele == '0|1' or allele == '1|0'):
            
                hetrozygous_phased += 1
                
                var_pos_het_list.append(var_pos)
                
                line_number_het_list.append(line_number)
                
                ps_index = format_genotype_splitted.index("PS")           #  index of phase set in values_genotype 
                id_block = values_genotype_splitted[ps_index]
                
                
                
                if first_het_variant:           # for the first het variant
                    first_het_variant = False
                    
                    allele_block = [allele]
                    var_pos_block = [int(var_pos)]
                    id_blocks.append(id_block)
                    
                    
                else:                              # for the rest of het variants
                    if id_block in id_blocks:
                        allele_block.append(allele)
                        var_pos_block.append(int(var_pos)) 

                    else:

                        # add previous block to the list of all blocks
                        allele_blocks.append(allele_block)
                        var_pos_blocks.append(var_pos_block)

                        # creat new phase block
                        allele_block = [allele]
                        var_pos_block = [int(var_pos)]
                        id_blocks.append(id_block)

                values_genotype_truth = line_parts[10]
                values_genotype_truth_splitted = values_genotype_truth.split(':')
                allele_truth = values_genotype_truth_splitted[gt_index]
                id_block_truth = values_genotype_truth_splitted[ps_index]
                #if allele_truth != './.':
                if allele_truth == '0|1' or allele_truth == '1|0':
                    allele_truth_dic[var_pos]= str(id_block_truth)+':'+str(allele_truth)
    
    
    # # for the last het variant, we  finish the last block.
    allele_blocks.append(allele_block)
    var_pos_blocks.append(var_pos_block)
    
    
          
    genomic_length_blocks = [] 
    
    for var_pos_block in var_pos_blocks:
        genomic_length_blocks.append(var_pos_block[-1]-var_pos_block[0])
        




    return genomic_length_blocks








if __name__ == "__main__":

    tech=argv[1]
    vcf_file_address_pre = "results_now/"#argv[1] # results_now/ccs/11/11_ccs_true_improved_.90.vcf 
    
    block_all=[]
    for chr in range(1,23):
        vcf_file_address=vcf_file_address_pre+tech+"/"+str(chr)+"/"+str(chr)+"_"+tech+"_true.vcf"
        print(vcf_file_address)
        genomic_length_blocks = read_vcf_file_with_truth(vcf_file_address)
        block_all= block_all+genomic_length_blocks
    


    values_sorted = sorted(block_all, reverse=True)
    csum = np.cumsum(values_sorted)
    
    n2 = int(sum(values_sorted)/2)
    csumn2 = min(csum[csum >= n2])
    ind = np.where(csum == csumn2)
    n50 = values_sorted[int(ind[0])]

    print(n50)




