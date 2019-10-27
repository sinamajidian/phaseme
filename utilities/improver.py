#!/usr/bin/env python3


import numpy as np




tresh_match_mismatch= 1.1
    
    
    
def read_vcf_file(vcf_file_address):
    
    
    """
    Reading the vcf file
    
    input: vcf file
    outputs: 
            header_lines_list: list of string. each string is a line that starts with #.
            var_pos_all: genomic position of all variants in the vcf file. 
            line_number_het_list: list of line numbers in the vcf file that are phased hetrozygous variant (needed in phasing)
            id_blocks: list of ids of phas blocks
            allele_blocks: list of list
            var_pos_blocks: list of list
            stats_vcf: [homozygous0_num, homozygous1_num, hetrozygous_nonphased, hetrozygous_phased, genomic_length_blocks, n50]
    
    
    """
    
    
    
    vcf_file = open(vcf_file_address,'r')

    header_lines_list=[]           # header lines  
    var_lines_list=[]              # needed for reporting improved VCF containing all lines except header.  (both homo and hetro) 
    
    var_pos_all=[]                # position of variants all blocks consequently. It can be used for converting variant index to genomic position (in pairs.txt)

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
        
        line_strip=line.strip() 
        if line_strip.startswith('#'):
            header_lines_list.append(line_strip)
            #sample_names=line_strip.split('\t')[9:11]            # last line of header contains sample name
        else:

            line_parts=line_strip.split('\t') 
            var_lines_list.append(line_parts)

            chrom = line_parts[0]

            if str(chrom_output)!=chrom:
                print(chrom_output,chrom)

            var_pos = int(line_parts[1])                           # genomic position of variants
            var_pos_all.append(var_pos)

            
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



    return header_lines_list, var_pos_all, line_number_het_list, id_blocks, allele_blocks, var_pos_blocks, stats_vcf




def read_file_mismatches(file_mismatches_address, id_blocks):

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
            
            line_parts = line.strip().split('\t')     # ['42081', '0', '42096', '0']
            
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


def read_file_pairs_forward(file_pairs_address, var_pos_all):
    
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





def read_file_pairs_forward_backward(file_pairs_address, var_pos_all):

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


    flip_list=[]
    cut_list_blocks=[]

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

            
            
            comparison_result = comparison_result_block[var_i]  # after a cut, it may use a newer version of this vsrisble
            if comparison_result == [[],[]]:
                out_pairs = '.'
            else:
                [matched_list, mismatched_list] = comparison_result    
                
                num_mismatched = len(mismatched_list)
                num_matched = len(matched_list)
                
                condition = 0
                condition_flip = 0
                condition_cut = 0
                
                if num_mismatched >= 2:
                    if num_matched == 0:
                        condition = 1
                        
                    elif num_matched/num_mismatched < tresh_match_mismatch :
                        condition=1
                        
                if condition:
                    
                    #print(len_mismatched_list, var_pos)
#                     ont_pop_prev = ont_pop_block[var_idx-1] 
#                     ont_pop_nex = ont_pop_block[var_idx+1] 
                    comparison_result_next_var = comparison_result_block[var_i+1]

                    num_matched_next = len(comparison_result_next_var[0])
                    num_mismatched_next = len(comparison_result_next_var[1])
#                     len_match_list_prev = len(ont_pop_prev[0])
#                     len_mismatch_list_prev = len(ont_pop_prev[1])

                    if num_matched_next>1.3*num_mismatched_next: #and len_match_list_prev>1.3*len_mismatch_list_prev: #  ?? if no pop information in last var
                        condition_flip=1
#                     if len_match_list_prev==0:
                    if  num_mismatched_next>=2:
                        condition_cut=1


                if condition_flip:
                    print('flip candidate', var_pos)
                    flip_list.append(var_pos)
                    
                if condition_cut and not condition_flip:
                    cut_pos = var_pos # cut_pos the starting position of new block
                    print('cut candidate', cut_pos) 
                    
                    cut_list_block.append(var_pos)
                    
                    lower_bound = var_pos
                    upper_bound = var_pos_block[-1]
                    
                    # after deciding the cut, it is not applied, so the number of block is not changed.
                    # but after deciding the match/mismatch information is updated for those varaint afterwards and removed for previose since new block is started practically. 
        
                    comparison_result_block = compare_phase_block_pop(allele_block, var_pos_block, pop_inf_dic, lower_bound, upper_bound)
        cut_list_blocks.append(cut_list_block)
    return flip_list, cut_list_blocks             



def improve_vcf_flip(lines_vcf_dic, flip_list):


    lines_vcf_dic_improved_flipping = lines_vcf_dic
    
    for var_pos in flip_list:
        
        line_parts=lines_vcf_dic[var_pos]

        format_genotype, values_genotype = line_parts[8:10]    # 'GT:GQ:DP:AF:GL:PS', '0|1:255:.:.:.,0,.:60780'

        values_genotype_splitted = values_genotype.split(':')
        format_genotype_splitted = format_genotype.split(':')

        gt_index = format_genotype_splitted.index("GT")           #  index of allele in  values_genotype 
        allele = values_genotype_splitted[gt_index]
        
        
        allele_flipped=str(1-int(allele[0]))+'|'+str(1-int(allele[2]))
        
        values_genotype_splitted=values_genotype.split(':')
        values_genotype_splitted[gt_index]=allele_flipped
        
        line_parts_flipped=line_parts
        
        line_parts_flipped[9]=':'.join(values_genotype_splitted)
        
        
        lines_vcf_dic_improved_flipping[var_pos]=line_parts_flipped
    
    return lines_vcf_dic_improved_flipping  
    
    

def improve_vcf_cut(lines_vcf_dic_improved_flipping, id_blocks, cut_list_blocks, var_pos_blocks):


    # allele_blocks
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
        

    lines_vcf_dic_improved_cut = lines_vcf_dic_improved_flipping

    sorted_var_pos = sorted(lines_vcf_dic_improved_cut.keys())
    for var_pos in sorted_var_pos: 
        
        line_parts = lines_vcf_dic_improved_cut[var_pos]

        format_genotype, values_genotype = line_parts[8:10]    # 'GT:GQ:DP:AF:GL:PS', '0|1:255:.:.:.,0,.:60780'
        values_genotype_splitted = values_genotype.split(':')
        format_genotype_splitted = format_genotype.split(':')

        gt_index = format_genotype_splitted.index("GT")           #  index of allele in  values_genotype 
        ps_index = format_genotype_splitted.index("PS")

        
        allele = values_genotype_splitted[gt_index]
        
        if (allele == '0|1' or allele == '1|0'):

            id_block = values_genotype_splitted[ps_index]

            block_id_updated = var_blockid_dic_updated[var_pos]

            values_genotype_splitted_updated = values_genotype.split(':')
            values_genotype_splitted_updated[ps_index] = block_id_updated

            line_parts_updated = line_parts
            line_parts_updated[9] = ':'.join(values_genotype_splitted_updated)

            lines_vcf_dic_improved_cut[var_pos] = line_parts_updated


    return lines_vcf_dic_improved_cut




def write_out_vcf(vcf_file_improved_address, lines_vcf_dic_improved_cut, header_lines_list ):


    # reading phased VCF file,  use grep before that, only 0|1 or 1|0. bi -allelic
    vcf_file_improved = open(vcf_file_improved_address,'w');  # phased_vcf_dic

    sorted_var_pos=sorted(lines_vcf_dic_improved_cut.keys())

    for header_line in header_lines_list:
        vcf_file_improved.write(header_line+'\n')
        
    for var_pos in sorted_var_pos: 
        
        line_parts=lines_vcf_dic_improved_cut[var_pos]
        vcf_file_improved.write('\t'.join(line_parts)+'\n')
        
    vcf_file_improved.close()
    
    return 1 




    



if __name__ == "__main__":

    

    """
    
    Input:
    
    Output:
    
    """


    chrom_output = 22
    vcf_file_address = 'data/'+str(chrom_output)+'/a22_10k.vcf' #tst.vcf' #
    
    # header_lines_list, lines_vcf_dic, var_pos_all, id_blocks, allele_blocks, var_pos_blocks  = read_vcf_file(vcf_file_address)
    
    header_lines_list, var_pos_all, line_number_het_list, id_blocks, allele_blocks, var_pos_blocks, stats_vcf  =  read_vcf_file(vcf_file_address)


    
    file_mismatches_address = 'data/'+str(chrom_output)+'/'+str(chrom_output)+'_report_mismatches.txt'
    comparison_result_blocks = read_file_mismatches(file_mismatches_address, id_blocks)
    
    file_pairs_address='data/'+str(chrom_output)+'/pairs_500_10k.txt'
    pop_inf_dic=read_file_pairs_forward(file_pairs_address, var_pos_all) #_backward for reporting both
    

    flip_list, cut_list_blocks = decide_flip_cut(id_blocks, allele_blocks, var_pos_blocks, comparison_result_blocks)

    lines_vcf_dic_improved_flipping = improve_vcf_flip(lines_vcf_dic, flip_list)  # it may contain homo vars
    
    lines_vcf_dic_improved_cut = improve_vcf_cut(lines_vcf_dic_improved_flipping, id_blocks, cut_list_blocks, var_pos_blocks)
    vcf_file_improved_address = 'data/'+str(chrom_output)+'/10k_improved_'+str(tresh_match_mismatch)+'.vcf' #tst.vcf' #

    write_out_vcf(vcf_file_improved_address, lines_vcf_dic_improved_cut, header_lines_list )
    
    phased_vcf_dic_improved2 = improve_vcf_cut(phased_vcf_dic_improved, hap_blocks, cut_dic)
    write_out_vcf(phased_vcf_dic_improved2, vcf_file_improved_address)

