#!/usr/bin/env python3


import numpy as np
from sys import argv



tresh_match_mismatch= 1.1
    
    

def read_vcf_file_with10x(vcf_file_address):    
    
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
    
    
    
    
    allele_10x_dic = {}
    
    vcf_file = open(vcf_file_address,'r')

    lines_list=[]                      #  lines  of phased vcf.  needed for reporting improved VCF 
    var_pos_het_list=[]                # position of phased hetrozygous variants all blocks consequently. 

    # The followings are for phased hetrozygous variants.
    id_blocks = []                 # list of list. Outer list corresponds to phase block. Inner list contains block_id
    
    allele_blocks = []             # list of list. Outer list corresponds to phase block. Inner list contains alleles of hetro variants 
    var_pos_blocks = []            # list of list. Outer list corresponds to phase block. Inner list contains genomic positions of hetro variants   

    
    line_number_het_list = []            # line number of phased hetrozygous variant. We need it for reporting improved version
    
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

            line_parts = line_strip.split('\t') 

            chrom = line_parts[0]

            if str(chrom_output)!=chrom:
                print(chrom_output,chrom)

            var_pos = int(line_parts[1])                           # genomic position of variants
            

            
            format_genotype, values_genotype = line_parts[8:10]    # 'GT:GQ:DP:AF:GL:PS', '0|1:255:.:.:.,0,.:60780'
            
            values_genotype_splitted = values_genotype.split(':')
            format_genotype_splitted = format_genotype.split(':')
            
            gt_index = format_genotype_splitted.index("GT")           #  index of allele in  values_genotype 
            
            allele = values_genotype_splitted[gt_index]
            
            
            # how should we handle '2' in allele ?  
            
            if './.' in allele:
                print("There is a vriant with genomic position "+str(var_pos)+" that is not genotyped. Remove it first.")
                exit(1)
                
            
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
    
    
    

                values_genotype_10x = line_parts[10]
                values_genotype_10x_splitted = values_genotype_10x.split(':')
                allele_10x = values_genotype_10x_splitted[gt_index]
                id_block_10x = values_genotype_10x_splitted[ps_index]
                if allele_10x != './.':
                    allele_10x_dic[var_pos]= str(id_block_10x)+':'+str(allele_10x)
    
    
    
    
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



    return lines_list, var_pos_het_list, line_number_het_list, id_blocks, allele_blocks, var_pos_blocks, stats_vcf, allele_10x_dic




#def read_file_mismatches(file_mismatches_address, id_blocks):
def read_file_mismatches_with10x(file_mismatches_address, id_blocks):

    file = open(file_mismatches_address,'r')
        
    
    first_variant = True
    block_ids = []       # this is only for checking consitency between phased VCF and the report file
    
    
    for line in file:
        
        if line.strip().startswith('#c'): # new chromosome
            comparison_result_blocks = []
        
        elif line.strip().startswith('#P'):  # new block
            
            if not first_variant:
                comparison_result_blocks.append(comparison_result_block)
                
            comparison_result_block = []
        
        else:
            
            first_variant=False
            
            line_parts_with10x = line.strip().split('\t')     # ['42081', '0', '42096', '0']
            
            line_parts = line_parts_with10x[2:]
            
            chrom = line_parts[0]
            var_pos = int(line_parts[1])
                        
            line_part_2 = line_parts[2].split(':')  # 60780:1|0

            block_id = line_part_2[0]
            
            if block_id not in block_ids:
                block_ids.append(block_id)
            
            allele = line_part_2[1]
            
            if line_parts[3] == '.':
                mismatched_list = []
            else:
                mismatched_list = line_parts[3].split(',')

            if line_parts[4] == '.':
                matched_list = []
            else:
                matched_list = line_parts[4].split(',') 
                
            comparison_result_block.append([matched_list, mismatched_list])
                        
            block_id_previous = block_id
            
    # add last block
    comparison_result_blocks.append(comparison_result_block)
    
    if block_ids != id_blocks:
        print("The block ids of phased VCF are not consitent with that of report file. The phased VCF file should be the same in all steps of package.\n")
                
            
                
            
    return comparison_result_blocks


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





def decide_flip_cut(id_blocks, allele_blocks, var_pos_blocks, comparison_result_blocks):
    
    flip_list = []
    cut_list_blocks = []

    for block_i, block_id  in enumerate(id_blocks):

        allele_block = allele_blocks[block_i]       # hap_block
        var_pos_block = var_pos_blocks[block_i] 
        
        
        # after deciding the cut, it is not applied, so the number of block is not changed.
        # but after deciding the match/mismatch information is updated for those varaint afterwards and removed for previose since new block is started practically. 
        
        comparison_result_block = comparison_result_blocks[block_i] # ont_pop_block

        
        cut_list_block=[]

        for var_i, var_pos in enumerate(var_pos_block):  # var_pos_block idc_block
            
            # var_i is the index of variant within block
            # var_pos is the genomic position of variant

            
            
            comparison_result = comparison_result_block[var_i]  # rightafter we decide a cut, we use an updated version of comparison_result_block

            if comparison_result == [[],[]]:  # there is no population information for this variant
                out_pairs = '.'
                
            else:
                [matched_list, mismatched_list] = comparison_result    
                
                num_mismatched = len(mismatched_list)
                num_matched = len(matched_list)
                
                condition_flip_cut = 0
                condition_flip = 0
                condition_cut = 0
                
                if num_mismatched >= 2:

                    [matched_next1, mismatched_next1] = comparison_result_block[var_i+1]
                    [matched_next2, mismatched_next2] = comparison_result_block[var_i+2]

                    num_matched_next1 = len(matched_next1)
                    num_mismatched_next1 = len(mismatched_next1)
                    
                    
                    if num_matched == 0 or num_matched/num_mismatched < tresh_match_mismatch:
                        if num_mismatched_next1 >=1 and (num_matched_next1 == 0 or num_matched_next1/num_mismatched_next1 < tresh_match_mismatch):
                            condition_cut = 1
                    
                    if var_pos in mismatched_next1 and var_pos in mismatched_next2:
                        condition_flip = 1

                        
                if condition_flip:
                    print('flip candidate', var_pos)
                    flip_list.append(var_pos)
                    
                if condition_cut and not condition_flip:
                    cut_pos = var_pos                 # cut_pos the starting position of new block
                    print('cut candidate', cut_pos) 
                    
                    cut_list_block.append(var_pos)
                    
                    lower_bound = var_pos
                    upper_bound = var_pos_block[-1]
                    
                    # after deciding the cut, it is not applied, so the number of block is not changed.
                    # but after deciding the match/mismatch information is updated for those varaint afterwards and removed for previose since new block is started practically. 
        
                    comparison_result_block = compare_phase_block_pop(allele_block, var_pos_block, pop_inf_dic, lower_bound, upper_bound)
        cut_list_blocks.append(cut_list_block)
    return flip_list, cut_list_blocks             





def improve_vcf_flip(lines_list, line_number_het_list, flip_list):
    
    
    # becareful about the concpet of deep copy. the change will affect the lines_list
    lines_list_improved_flipping = lines_list
    
    for var_pos_flip in flip_list:
        
       
        var_index_flip = var_pos_het_list.index(var_pos_flip)
        
        var_line_number_flip = line_number_het_list[var_index_flip]    # 1-based line number 
                
        line = lines_list[var_line_number_flip-1]                      # 0-based  list

        line_parts = line.split('\t')

        format_genotype, values_genotype = line_parts[8:10]    # 'GT:GQ:DP:AF:GL:PS', '0|1:255:.:.:.,0,.:60780'

        values_genotype_splitted = values_genotype.split(':')
        format_genotype_splitted = format_genotype.split(':')

        gt_index = format_genotype_splitted.index("GT")           #  index of allele in  values_genotype 
        allele = values_genotype_splitted[gt_index]
        
        allele_flipped = str(1-int(allele[0]))+'|'+str(1-int(allele[2]))
        values_genotype_splitted[gt_index] = allele_flipped
        
        line_parts[9]=':'.join(values_genotype_splitted)
        
        line = '\t'.join(line_parts)
        
        lines_list_improved_flipping[var_line_number_flip-1] = line


        
    return lines_list_improved_flipping  
   
    


def improve_vcf_cut(lines_list_improved_flipping, id_blocks, cut_list_blocks, var_pos_blocks):


    var_blockid_dic_updated = {}   # key: var_pos (genomic position of variant) value: block_id   after enforcing cuts!
    
    
    for block_i, block_id in enumerate(id_blocks):
        
        cut_list_block = cut_list_blocks[block_i]
        var_pos_block = var_pos_blocks[block_i]
        
        boundries_list = [var_pos_block[0]]+cut_list_block+[var_pos_block[-1]]   # including start, end, and the cut list of  each block
        
        for i in range(1,len(boundries_list)):
            
            prev_cut_pos = boundries_list[i-1]
            next_cut_pos = boundries_list[i]
            
            blockid = prev_cut_pos
            
            for var_pos in var_pos_block:
                
                if var_pos >= prev_cut_pos and var_pos < next_cut_pos:
                    var_blockid_dic_updated[var_pos] = str(prev_cut_pos)
        # for last variant in the block
        var_blockid_dic_updated[boundries_list[i]] = str(boundries_list[i-1])
        
    
    lines_list_improved_cut = lines_list_improved_flipping

    for var_i in range(len(var_pos_het_list)): 
        
        var_pos = var_pos_het_list[var_i]
        line_number = line_number_het_list[var_i]   # 1-based line number 
         
        line = lines_list_improved_cut[line_number-1]            # 0-based  list
         
        line_parts = line.split('\t')

        format_genotype, values_genotype = line_parts[8:10]    # 'GT:GQ:DP:AF:GL:PS', '0|1:255:.:.:.,0,.:60780'

        values_genotype_splitted = values_genotype.split(':')
        format_genotype_splitted = format_genotype.split(':')

        gt_index = format_genotype_splitted.index("GT")           #  index of allele in  values_genotype
        ps_index = format_genotype_splitted.index("PS")

        allele = values_genotype_splitted[gt_index]

        block_id_updated = var_blockid_dic_updated[var_pos]

        values_genotype_splitted[ps_index] = block_id_updated

        line_parts[9] = ':'.join(values_genotype_splitted)

        lines_list_improved_cut[line_number-1]  = '\t'.join(line_parts)


    return lines_list_improved_cut


    

    
    
def write_out_vcf(lines_list, lines_list_improved_cut):

    vcf_file_improved = open(vcf_file_improved_address,'w');  # phased_vcf_dic

    for header_line in lines_list_improved_cut:
        vcf_file_improved.write(header_line+'\n')
        

    vcf_file_improved.close()
    
    return 1 





if __name__ == "__main__":

    

    """
    
    Input:
    
    Output:
    
    """


    chrom_output = 22
    vcf_file_address = 'uploaded_dropbox/input/'+str(chrom_output)+'_ont_10x.vcf' #tst.vcf' #
        
    lines_list, var_pos_het_list, line_number_het_list, id_blocks, allele_blocks, var_pos_blocks, stats_vcf, allele_10x_dic  =  read_vcf_file_with10x(vcf_file_address)


    
    #file_mismatches_address = 'uploaded_dropbox/'+str(chrom_output)+'_report_mismatches.txt'
    file_mismatches_address = 'uploaded_dropbox/22_report_mismatches_onlyforward.txt'
    comparison_result_blocks = read_file_mismatches_with10x(file_mismatches_address, id_blocks)
    
    file_pairs_address = 'uploaded_dropbox/input/'+str(chrom_output)+'_pairs.txt'
    pop_inf_dic = read_file_pairs_forward(file_pairs_address, var_pos_het_list) #_backward for reporting both
    

    flip_list, cut_list_blocks = decide_flip_cut(id_blocks, allele_blocks, var_pos_blocks, comparison_result_blocks)

    
    
    # becareful about the concpet of deep copy.  Changes will affect the lines_list
    lines_list_improved_flipping = improve_vcf_flip(lines_list, line_number_het_list, flip_list)  # it may contain homo vars
    
    
    
    lines_list_improved_cut= improve_vcf_cut(lines_list_improved_flipping, id_blocks, cut_list_blocks, var_pos_blocks)

    vcf_file_improved_address = 'uploaded_dropbox/input/'+str(chrom_output)+'_ont_10x_improved.vcf' 

    write_out_vcf(vcf_file_improved_address, lines_list_improved_cut)
    



