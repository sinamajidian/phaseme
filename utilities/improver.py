#!/usr/bin/env python3

#!/usr/bin/env python3


import numpy as np




tresh_match_mismatch= 1.1
    
def read_vcf_file(vcf_file_address):
    vcf_file = open(vcf_file_address,'r')

    header_lines_list=[]           # header lines  
    #var_lines_list=[]              # needed for reporting improved VCF containing all lines except header.  (both homo and hetro) 
    
    
    var_pos_all=[]                # position of variants all blocks consequently. It can be used for converting variant index to genomic position (in pairs.txt)

    # The followings are for hetrozygous variants.
    id_blocks = []                 # list of list. Outer list corresponds to phase block. Inner list contains block_id
    
    allele_blocks = []             # list of list. Outer list corresponds to phase block. Inner list contains alleles of hetro variants 
    var_pos_blocks = []            # list of list. Outer list corresponds to phase block. Inner list contains genomic positions of hetro variants   
    
    
    lines_vcf_dic = {}            # key: var_pos  value: line (both het and homo) 
    #line_number_het_list = []            # line number of hetrozygous variant. We need it for reporting improved version

    first_het_variant = True
    #line_number=0
    for line in vcf_file:
        
        #line_number += 1
        
        line_strip=line.strip() 
        if line_strip.startswith('#'):
            header_lines_list.append(line_strip)
            #sample_names=line_strip.split('\t')[9:11]            # last line of header contains sample name
        else:

            line_parts=line_strip.split('\t') 
            
            
            
            
            #var_lines_list.append(line_parts)

            chrom = line_parts[0]

            if str(chrom_output)!=chrom:
                print(chrom_output,chrom)

            var_pos = int(line_parts[1])                           # genomic position of variants
            var_pos_all.append(var_pos)
            
            lines_vcf_dic[var_pos]=line_parts

            
            format_genotype, values_genotype = line_parts[8:10]    # 'GT:GQ:DP:AF:GL:PS', '0|1:255:.:.:.,0,.:60780'
            
            values_genotype_splitted = values_genotype.split(':')
            format_genotype_splitted = format_genotype.split(':')
            
            gt_index = format_genotype_splitted.index("GT")           #  index of allele in  values_genotype 
            ps_index = format_genotype_splitted.index("PS")           #  index of phase set in values_genotype 
            
            
            allele = values_genotype_splitted[gt_index]
            
            
            # how should we handle '2' in allele ?  
            
            if './.' in allele:
                print("There is a vriant with genomic position "+str(var_pos)+" that is not genotyped. Remove it first.")
                exit(1)
                
            elif '/' in allele:
                print("There is a vriant with genomic position "+str(var_pos)+" that is not phased. Remove it first.")            
                exit(1)
                
            elif (allele == '0|0' or allele == '1|1'):
                # homozygous variant
                pass
            
            elif (allele == '0|1' or allele == '1|0'):
                #line_number_het_list.append(line_number)
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
                    
            



    #return header_lines_list, var_lines_list, var_pos_all, line_number_het_list, id_blocks, allele_blocks, var_pos_blocks
    return header_lines_list, lines_vcf_dic, var_pos_all, id_blocks, allele_blocks, var_pos_blocks





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
        print("The block ids of phased VCF are not consitent with that of report file. The phased VCF file should be the same for all steps of code.\n")
                
            
                
            
    return comparison_result_blocks


def read_file_pairs_forward(file_pairs_address, var_pos_all):

    # identical opposite
    
    file_pairs = open(file_pairs_address,'r'); 
    
    pop_inf_dic = {}                              # key variant index, value two lists 
    for line in file_pairs:
        
        line_parts = line.strip().split('\t')     # ['42081', '0', '42096', '0']

        
        
        snv1_index = int(line_parts[0])
        snv1_pos = var_pos_all[snv1_index-1]
        allele1_snv1 = int(line_parts[1])
        
        snv2_index = int(line_parts[2])
        snv2_pos = var_pos_all[snv2_index-1]
        allele1_snv2 = int(line_parts[3])
        

        snvs_pos = [snv1_pos, snv2_pos]
        host_pos = max(snvs_pos)
        guest_pos = min(snvs_pos)
        
        if host_pos not in pop_inf_dic.keys():
            
            # first list identical phase
            # second list opposite phase
            pop_inf_dic[host_pos] = [[], []]
            
        if allele1_snv1 == allele1_snv2:
            pop_inf_dic[host_pos][0].append(guest_pos)  
        else:
            pop_inf_dic[host_pos][1].append(guest_pos) 
            
    return pop_inf_dic






def read_file_pairs_forward_backward(file_pairs_address, var_pos_all):

    # identical opposite
    
    file_pairs = open(file_pairs_address,'r'); 
    
    pop_inf_dic = {}                              # key variant index, value two lists 
    for line in file_pairs:
        
        line_parts = line.strip().split('\t')     # ['42081', '0', '42096', '0']

        
        
        snv1_index = int(line_parts[0])
        snv1_pos = var_pos_all[snv1_index-1]
        allele1_snv1 = int(line_parts[1])
        
        snv2_index = int(line_parts[2])
        snv2_pos = var_pos_all[snv2_index-1]
        allele1_snv2 = int(line_parts[3])
          
        
        host_pos = snv1_pos
        guest_pos = snv2_pos
        
        if host_pos not in pop_inf_dic.keys():
            
            # first list identical phasing 
            # second list opposite phasing 
            pop_inf_dic[host_pos] = [[], []]
            
        if allele1_snv1 == allele1_snv2:
            pop_inf_dic[host_pos][0].append(guest_pos)  
        else:
            pop_inf_dic[host_pos][1].append(guest_pos) 
            
            
            
        host_pos = snv2_pos
        guest_pos = snv1_pos
        
        if host_pos not in pop_inf_dic.keys():
            
            # first list identical phasing 
            # second list opposite phasing 
            pop_inf_dic[host_pos] = [[], []]
            
        if allele1_snv1 == allele1_snv2:
            pop_inf_dic[host_pos][0].append(guest_pos)  
        else:
            pop_inf_dic[host_pos][1].append(guest_pos) 
            
    return pop_inf_dic




def compare_phase_block_pop(allele_block, var_pos_block, pop_inf_dic, lower_bound, upper_bound):

    
    # the results of comparison between sample (phased VCF) and pairs (the population information, pop_inf_dic)
    
    comparison_result_block = []
    for var_i in range(len(allele_block)):

        #var_i is the index of variant within the block of phased vcf

        var_pos = var_pos_block[var_i]
        allele=allele_block[var_i]



        if var_pos >= lower_bound and  var_pos <= upper_bound:

            #alleles=hap_block[var_i]
            #var_idx=idc_block[var_i]
            #var_i is the index of variant within block
            #var_idx is the index of vriant globally in the VCF file 
            if var_pos in pop_inf_dic.keys():
    
                [vars_identical_phase, vars_opposite_phase]=pop_inf_dic[var_pos]  
                matched_identical_phase_list=[]
                mismatched_identical_phase_list=[]
                matched_opposite_phase_list=[]
                mismatched_opposite_phase_list=[]
                #for sim_idx in vars_identical_phase: # sim shows the relation between two elements of a pop pair
                
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
                    
                matched_list = matched_identical_phase_list+matched_opposite_phase_list
                mismatched_list = mismatched_identical_phase_list + mismatched_identical_phase_list
                
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


    chrom_output=19
    vcf_file_address = 'data/'+str(chrom_output)+'/10k_out.vcf' #tst.vcf' #

    #hap_blocks, idc_blocks, var_pos_list, phased_vcf_dic, header_lines_list = read_vcf_file(vcf_file_address)
    #header_lines_list, var_lines_list, var_pos_all, line_number_het_list, id_blocks, allele_blocks, var_pos_blocks = read_vcf_file(vcf_file_address)
    header_lines_list, lines_vcf_dic, var_pos_all, id_blocks, allele_blocks, var_pos_blocks  = read_vcf_file(vcf_file_address)
    
    
    file_mismatches_address='data/'+str(chrom_output)+'/10k_'+str(chrom_output)+'_report_mismatches.txt'

    comparison_result_blocks=read_file_mismatches(file_mismatches_address, id_blocks)
    
    file_pairs_address='data/'+str(chrom_output)+'/10k_pairs.txt'
    pop_inf_dic=read_file_pairs_forward(file_pairs_address, var_pos_all) #_backward for reporting both
    

    flip_list, cut_list_blocks = decide_flip_cut(id_blocks, allele_blocks, var_pos_blocks, comparison_result_blocks)

    lines_vcf_dic_improved_flipping = improve_vcf_flip(lines_vcf_dic, flip_list)  # it may contain homo vars
    
    lines_vcf_dic_improved_cut = improve_vcf_cut(lines_vcf_dic_improved_flipping, id_blocks, cut_list_blocks, var_pos_blocks)
    vcf_file_improved_address = 'data/'+str(chrom_output)+'/10k_improved_'+str(tresh_match_mismatch)+'.vcf' #tst.vcf' #

    write_out_vcf(vcf_file_improved_address, lines_vcf_dic_improved_cut, header_lines_list )
    
#     phased_vcf_dic_improved2 = improve_vcf_cut(phased_vcf_dic_improved, hap_blocks, cut_dic)
#     write_out_vcf(phased_vcf_dic_improved2, vcf_file_improved_address)

