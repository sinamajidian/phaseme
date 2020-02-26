#!/usr/bin/env python3




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
        
    values_sorted = sorted(genomic_length_blocks, reverse=True)
    csum = np.cumsum(values_sorted)
    
    n2 = int(sum(values_sorted)/2)
    csumn2 = min(csum[csum >= n2])
    ind = np.where(csum == csumn2)
    n50 = values_sorted[int(ind[0])]

        
    
    stats_vcf = [homozygous0_num, homozygous1_num, hetrozygous_nonphased, hetrozygous_phased, genomic_length_blocks, n50]



    return lines_list, var_pos_het_list, line_number_het_list, id_blocks, allele_blocks, var_pos_blocks, stats_vcf, allele_truth_dic, chrom


def read_file_pairs_forward(file_pairs_address):
    
    """
    
    Reading the pairs.txt file and save it in a dictinary only once.

    input: pairs.txt
    Each line of input file has three part 
        position of var1, position of  var2, realtion between phasing of var1 and that of var2
        10 20 identical
        
        Each pair reports only one in output dictionary. We report the example for the key   20:[[10],[]]
        Becuase, in phase block we go variant by vairant. So the phasing of previous variants are important not later's
        
    output: a dictinary, key: genomic position of variant, value: [[],[]] 
    
            first list: the genomic position of variant with identical phasing 
            second list: the genomic position of variant with opposite phasing
    
    """
    
    
    file_pairs = open(file_pairs_address,'r'); 
    
    pop_inf_dic = {}                              # key variant index, value two lists 
    for line in file_pairs:
        
        line_parts = line.strip().split('\t')     # ['42081', '0', '42096', '0']

        
        
        snv1_pos = int(line_parts[0]) 
        snv2_pos = int(line_parts[1])                

        relation_phasing = line_parts[2]
        
        snvs_pos = [snv1_pos, snv2_pos]
        host_pos = max(snvs_pos)
        guest_pos = min(snvs_pos)
        
        if host_pos not in pop_inf_dic.keys():
            
            # first list identical phase
            # second list opposite phase
            pop_inf_dic[host_pos] = [[], []]
            
        if relation_phasing == 'identical':
            pop_inf_dic[host_pos][0].append(guest_pos)  
        elif relation_phasing == 'opposite':
            pop_inf_dic[host_pos][1].append(guest_pos) 
            
    return pop_inf_dic  # [vars_identical_phase, vars_opposite_phase]=pop_inf_dic[var_pos]





def read_file_pairs_forward_backward(file_pairs_address):

    """
    see read_file_pairs_forward
    
    The only difference is that a pair is reported twic. 
    
    """
    
    file_pairs = open(file_pairs_address,'r'); 
    
    pop_inf_dic = {}                              # key variant index, value two lists 
    for line in file_pairs:
        
        line_parts = line.strip().split('\t')     # ['42081', '0', '42096', '0']

         
        snv1_pos = int(line_parts[0]) 
        snv2_pos = int(line_parts[1])                

        relation_phasing = line_parts[2]
        
        host_pos = snv1_pos
        guest_pos = snv2_pos
        
        if host_pos not in pop_inf_dic.keys():
            
            # first list identical phase
            # second list opposite phase
            pop_inf_dic[host_pos] = [[], []]
            
        if relation_phasing == 'identical':
            pop_inf_dic[host_pos][0].append(guest_pos)  
        elif relation_phasing == 'opposite':
            pop_inf_dic[host_pos][1].append(guest_pos) 
   
            
        host_pos = snv2_pos
        guest_pos = snv1_pos
        
        if host_pos not in pop_inf_dic.keys():
            
            # first list identical phase
            # second list opposite phase
            pop_inf_dic[host_pos] = [[], []]
            
        if relation_phasing == 'identical':
            pop_inf_dic[host_pos][0].append(guest_pos)  
        elif relation_phasing == 'opposite':
            pop_inf_dic[host_pos][1].append(guest_pos) 
            
    return pop_inf_dic  # [vars_identical_phase, vars_opposite_phase]=pop_inf_dic[var_pos]

      



def compare_phase_block_pop(allele_block, var_pos_block, pop_inf_dic, lower_bound, upper_bound):

    
    """
    Compare alleles of phased VCF with pairs (from the population information)
    
    input:  allele_block, var_pos_block    (phased VCF)
            pop_inf_dic 
    
    output: a list of [[],[]] 
    
                        Inner list correspond to a variant (host). 
                        first inner list: positons of variants that are matched with the host variant (based on population information)
                                         the status of match can be either identical or opposite phasing.                
                        second inner list mismatched
    
    
    Exmpale 1:  pair from population information:  10 20 identical
                if phased vcf   10 0|1    20 0|1, the 10 and 20 are matched.
                if phased vcf   10 0|1    20 1|0, the 10 and 20 are mismatched.

    Exmpale 2:  pair from population information:  10 20 opposite
                if phased vcf   10 0|1    20 0|1, the 10 and 20 are mismatched.
                if phased vcf   10 0|1    20 1|0, the 10 and 20 are matched.


    """
    
    
    # the results of comparison between 
    
    comparison_result_block = []
    for var_i in range(len(allele_block)):

        #var_i is the index of variant within the block of phased vcf

        var_pos = var_pos_block[var_i]
        allele = allele_block[var_i]

        if var_pos >= lower_bound and  var_pos <= upper_bound:

            #alleles=hap_block[var_i]
            #var_idx=idc_block[var_i]
            #var_i is the index of variant within block
            #var_idx is the index of vriant globally in the VCF file 
            if var_pos in pop_inf_dic.keys():
    
                [vars_identical_phase, vars_opposite_phase] = pop_inf_dic[var_pos]  
        
                
                matched_identical_phase_list = []
                mismatched_identical_phase_list = []
                matched_opposite_phase_list = []
                mismatched_opposite_phase_list = []
                #for sim_idx in vars_identical_phase: # sim shows the relation between two elements of a pop pair
                
                
                # The differene between the two following for is the comparing the allele_var_guest and  allele (flliped)
                
                for var_pos_identical_phase in vars_identical_phase: # Those variant that have the same phasing with var_pos

                    if var_pos_identical_phase >= lower_bound and var_pos_identical_phase <= upper_bound:
                        try:         # if the SNP  of pairs is in this block of phased VCF
                            
                            allele_var_guest= allele_block[var_pos_block.index(var_pos_identical_phase)]
                                                   
                            if allele_var_guest == allele: 
                                matched_identical_phase_list.append(var_pos_identical_phase)
                            else:
                                mismatched_identical_phase_list.append(var_pos_identical_phase) 

                        except:
                            pass

                    
                for var_pos_opposite_phase in vars_opposite_phase: # Those variant that have the same phasing with var_pos
                    
                    if var_pos_opposite_phase >= lower_bound and  var_pos_opposite_phase <= upper_bound:
                        try:         # if the SNP  of pairs is in this block of phased VCF        
                            allele_var_guest = allele_block[var_pos_block.index(var_pos_opposite_phase)]
                                                
                            if allele_var_guest == str(1-int(allele[0]))+'|'+str(1-int(allele[2])):                   
                                matched_opposite_phase_list.append(var_pos_opposite_phase)
                            else:
                                mismatched_opposite_phase_list.append(var_pos_opposite_phase) 
                        except:
                            pass                
                    
                matched_list = matched_identical_phase_list + matched_opposite_phase_list
                mismatched_list = mismatched_identical_phase_list + mismatched_opposite_phase_list
                
                comparison_result = [matched_list, mismatched_list]
            else:
                comparison_result = [[],[]]
        else:
            comparison_result = [[],[]]
            
        comparison_result_block.append(comparison_result)
                
        #print('comparison is done')
    return  comparison_result_block




def report_comparison_withtruth(report_out_address, comparison_result_blocks, chrom, allele_truth_dic):

    file_report= open(report_out_address,'w'); 
    file_report.write("#chr\t var_pos\t blockid:alleles_ont \t mismatched_pair \t  matched_pair  \n")
    
    qual_blocks=[]
    
    for block_i, block_id  in enumerate(id_blocks):
        
        allele_block = allele_blocks[block_i]
        var_pos_block = var_pos_blocks[block_i]
        
                
    
        file_report.write('\t'.join(["#Phase Block with ID  ",str(block_id)])+"\n")

        comparison_result_block = comparison_result_blocks[block_i]
        
        qual_block=[]
        
        
        for var_i in range(len(allele_block)): # var_i internal index of this phase block

            allele=allele_block[var_i]
            var_pos=var_pos_block[var_i]


            list_write=[]

            
            if var_pos in allele_truth_dic:
                list_write.append(str(allele_truth_dic[var_pos]))
                
                allele_truth= allele_truth_dic[var_pos][-3:]
                allele_truth_revert = str(1-int(allele_truth[0]))+'|'+str(1-int(allele_truth[2]))
                
                if allele == allele_truth_dic[var_pos][-3:]:
                    sample_vs_truth = 'same'
                elif allele == allele_truth_revert:
                    sample_vs_truth = 'flipped'
                    
                list_write.append(sample_vs_truth)
                
                
            else:
                list_write.append('.')
                list_write.append('.')
                
                
                
            
            list_write.append(str(chrom))
            list_write.append(str(var_pos))
            
            list_write.append(str(block_id)+':'+str(allele))

            [matched_list, mismatched_list]=comparison_result_block[var_i]
            
            if len(matched_list)+len(mismatched_list):                # if we have population information
                if len(matched_list):
                    qual = len(matched_list)/(len(mismatched_list)+len(matched_list))
                else:
                    qual = 0
                qual_block.append(qual)
                
                
            if len(mismatched_list):
                list_write.append(','.join([str(i) for i in mismatched_list]))
            else: 
                list_write.append('.')

            if len(matched_list):
                list_write.append(','.join([str(i) for i in matched_list]))
            else: 
                list_write.append('.')

            file_report.write('\t'.join(list_write)+"\n")
        
        
        qual_blocks.append(qual_block)
        
        #file_report.write("\n") # new  block

    file_report.close()

    return qual_blocks




def report_qc(report_qc_address, id_blocks, qual_blocks, allele_blocks, stats_vcf):
    
    [homozygous0_num, homozygous1_num, hetrozygous_nonphased, hetrozygous_phased, genomic_length_blocks, n50] = stats_vcf
    
    file_report_qc= open(report_qc_address, 'w'); 
    
    #file_report_qc.write("##per chr:homozygous variants  .\n")
    file_report_qc.write("##Block_i:\tStart_pos\tblock_length\tnum_phased_SNV\t\tblock_quality \n \n")
    file_report_qc.write("Quality report for chromosome "+str(chrom)+":\n \n")
    
    
    #file_report_qc.write("Number of homozygous variants is "+str(homozygous0_num+homozygous1_num)+", including "+str(homozygous0_num)+" variants with 0|0 and "+str(homozygous1_num)+" variants with 1|1:.\n")
    #file_report_qc.write("Number of non-phased heterozygous variants is "+str(hetrozygous_nonphased)+".\n")
    # file_report_qc.write("Number of heterozygous variants in all blocks is "+str(hetrozygous_phased)+".\n\n")
    
    # file_report_qc.write("N50 of phase blocks is "+str(n50)+"bp.\n\n")
    
    for block_i, block_id  in enumerate(id_blocks):
        
        qual_block=qual_blocks[block_i]
        allele_block=allele_blocks[block_i]
        genomic_length_block = genomic_length_blocks[block_i]
        
        if len(qual_block):
            q = round(np.mean(qual_block),5)
        else:
            q= 'NaN'
            
            
            
        file_report_qc.write(str(block_i)+":\t"+block_id+"\t"+str(genomic_length_block)+"\t"+str(len(allele_block))+"\t"+str(q))

        file_report_qc.write("\n") # new  block

    file_report_qc.close()
    return 1















if __name__ == "__main__":

    

    """
    
    Input: phased vcf
           pairs.txt
           
    
    Output: report qc
            report match and mismatches
    
    """

    #type_data ='ont'
    chrom = argv[3]

    #vcf_file_address = type_data+'_son/'+str(chrom)+'/'+str(chrom)+'_'+type_data+'_true.vcf'  # truth
    vcf_file_address = argv[1]
    
    lines_list, var_pos_het_list, line_number_het_list, id_blocks, allele_blocks, var_pos_blocks, stats_vcf, allele_truth_dic, chrom = read_vcf_file_with_truth(vcf_file_address)
    

    #file_pairs_address= type_data+'_son/'+str(chrom)+'/'+str(chrom)+'_'+type_data+'_pairs_500_0.9.txt'  #'_pairs_500_0.98.txt'
    file_pairs_address= argv[2] 
    #file_pairs_address=  argv[1][:-4]+'_pairs_500_0.9.txt'  #'_pairs_500_0.98.txt'
    
    
    
    for Forward_or_both in  ['both', 'forward']:

        if Forward_or_both == 'forward':
            pop_inf_dic = read_file_pairs_forward(file_pairs_address)         

        if Forward_or_both == 'both':                       # reproting pairs both forward and backward
            pop_inf_dic = read_file_pairs_forward_backward(file_pairs_address)      




        # comparing input phased vcf (hap_blocks_sample1_dic) with pop (short_blocks_dic_pos)
        comparison_result_blocks=[]
        for block_i in range(len(id_blocks)):

            allele_block = allele_blocks[block_i]
            var_pos_block = var_pos_blocks[block_i]

            # lower_bound and upper_bound is used for limit the span of comparison.
            # Here, there is no limit. So, we set it from the begining to end of genomic position of block.
            # The limit is used for updating the comparison result after cutting the block (in improver.py).
            lower_bound = var_pos_block[0]         
            upper_bound = var_pos_block[-1]

            comparison_result_block = compare_phase_block_pop(allele_block, var_pos_block, pop_inf_dic, lower_bound, upper_bound)
            comparison_result_blocks.append(comparison_result_block)



        report_out_address=argv[1][:-4]+'_report_mismatches_'+Forward_or_both+'.txt'
#type_data+'_son/'+str(chrom)+'/'+str(chrom)+'_report_mismatches_'+Forward_or_both+'.txt'
        qual_blocks = report_comparison_withtruth(report_out_address, comparison_result_blocks, chrom, allele_truth_dic)

        report_qc_address=argv[1][:-4]+'_report_qc_'+Forward_or_both+'.txt' # type_data+'_son/'+str(chrom)+'/'+str(chrom)+'_report_qc_'+Forward_or_both+'.txt'
        report_qc(report_qc_address, id_blocks, qual_blocks, allele_blocks, stats_vcf)

  

