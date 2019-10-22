#!/usr/bin/env python3


import numpy as np




tresh_match_mismatch= 1.1
    
    
    
def read_vcf_file(vcf_file_address):
    vcf_file = open(vcf_file_address,'r')

    hap_blocks={} # key: id of phase block (i.e. phase st), value=allele_1 
    idc_blocks={} # key: id of phase block (i.e. phase st), value=var_idx   

    phased_vcf_dic={}

    var_pos_list=[] # all blocks consequently 
    header_lines_list=[] 

    var_idx=0 # for the first variant, it's 1. 1-based indexing 
    for line in vcf_file:
        line_strip=line.strip() 
        if line_strip.startswith('#'):
            header_lines_list.append(line_strip)
            sample_names=line_strip.split('\t')[9:11]# last line of header contains sample name
        else:
            var_idx+=1
            line_parts=line_strip.split('\t') 

            gt_flags, sample1=line_parts[8:10] # sample of interest that is merged        
            chrom=line_parts[0]

            if str(chrom_output)!=chrom:
                print(chrom_output,chrom)

            var_pos=int(line_parts[1])
            var_pos_list.append(var_pos)
            #if var_pos < 7398830: #60780

            phased_vcf_dic[var_idx]=line_parts


            sample1_split=sample1.split(":")
            #gt_flags_split=gt_flags.split(":")
            gt_id_flags=gt_flags.split(":").index("GT")
            ps_id_flags=gt_flags.split(":").index("PS")
            allele_sample1=sample1_split[gt_id_flags]        
            if '.' not in allele_sample1:
                block_id=int(sample1_split[ps_id_flags])
                if block_id in hap_blocks:
                    hap_blocks[block_id].append(allele_sample1) # append new variants to the existing phase block
                    idc_blocks[block_id].append(var_idx)
                else:
                    hap_blocks[block_id]=[allele_sample1] # creat new phase block
                    idc_blocks[block_id]=[var_idx]

            else:
                print('There is a vriant that ONT is not phased.. removed it first')

    return hap_blocks, idc_blocks, var_pos_list, phased_vcf_dic, header_lines_list








def read_file_mismatches(file_mismatches_address, var_pos_list):

    file = open(file_mismatches_address,'r')
    
    ont_pop_blocks={} # key phase block
    
    # key variant index, value two lists 
    for line in file:
        if not line.strip().startswith('#'):
            
            line_parts=line.strip().split('\t')     # ['42081', '0', '42096', '0']
            
            chrom=line_parts[0]
            var_pos=int(line_parts[1])
            
            var_idx=var_pos_list.index(int(var_pos))+1  # 1 basesd indexing
            
            line_part_2=line_parts[2].split(':')  # 60780:1|0

            block_id=int(line_part_2[0])
            allele=line_part_2[1]
            
            
            if line_parts[3] == '.':
                mismatch_list_pos=[]
            else:
                mismatch_list_pos=line_parts[3].split(',')

            if line_parts[4] == '.':
                match_list_pos=[]
            else:
                match_list_pos=line_parts[4].split(',')
                
            
            #match_list_len=len(match_list)
            
            match_list=[]
            for pos in match_list_pos:
                var_idx_match=var_pos_list.index(int(pos))+1 # 1 basesd indexing
                match_list.append(var_idx_match) 
            
            mismatch_list=[]
            for pos in mismatch_list_pos:
                var_idx_mismatch=var_pos_list.index(int(pos))+1 # 1 basesd indexing
                mismatch_list.append(var_idx_mismatch)  

            ont_pop=[match_list, mismatch_list]
            
            
            if block_id in ont_pop_blocks:
                ont_pop_blocks[block_id][var_idx]=ont_pop
            else:
                ont_pop_blocks[block_id]={var_idx:ont_pop}
                
            
    return ont_pop_blocks




def decide_flip_cut(hap_blocks, idc_blocks, var_pos_list, ont_pop_blocks):



    flip_list=[]
    cut_dic={}

    for block_id, hap_block in hap_blocks.items():
        cut_dic[block_id]=[]

        idc_block=idc_blocks[block_id]
        ont_pop_block=ont_pop_blocks[block_id]

        for var_i, var_idx in enumerate(idc_block): 
            #var_i is the index of variant within block
            #var_idx is the index of vriant globally in the VCF file  1-based

            ont_pop = ont_pop_block[var_idx]  # after a cut, it may use a newer version of this vsrisble
            if ont_pop==[[],[]]:
                out_pairs='.'
            else:
                [match_list,mismatch_list]=ont_pop    
                len_mismatch_list=len(mismatch_list)
                len_match_list=len(match_list)
                condition=0
                condition_flip=0
                condition_cut=0
                if len_mismatch_list>=2:
                    if len_match_list==0:
                        condition=1
                    elif len_match_list/len_mismatch_list < tresh_match_mismatch :
                        condition=1
                if condition:
                    #print(len_mismatch_list,var_idx,var_pos_list[var_idx-1])
                    ont_pop_prev = ont_pop_block[var_idx-1] 
                    ont_pop_nex = ont_pop_block[var_idx+1] 

                    len_match_list_nex=len(ont_pop_nex[0])
                    len_mismatch_list_nex=len(ont_pop_nex[1])
                    len_match_list_prev=len(ont_pop_prev[0])
                    len_mismatch_list_prev=len(ont_pop_prev[1])

                    if len_match_list_nex>1.3*len_mismatch_list_nex and len_match_list_prev>1.3*len_mismatch_list_prev: #  ?? if no pop information in last var
                        condition_flip=1
                    if len_match_list_prev==0:
                        if  len_mismatch_list_nex>=2:
                            condition_cut=1
                    elif len_match_list/len_mismatch_list < tresh_match_mismatch:
                        condition_cut=1

                if condition_flip:
                    print('flip candidate', var_idx, var_pos_list[var_idx-1])
                    flip_list.append(var_idx)
                if condition_cut and not condition_flip:
                    cut_pos=var_pos_list[var_idx-1] # cut_pos the starting position of new block
                    print('cut candidate', cut_pos, var_i, var_idx)   
                    cut_dic[block_id].append(var_idx)
                    interval=range(var_i,len(hap_block))




                    ont_pop_block=compare_phase_block_pop(hap_block,idc_block,pop_inf_dic,interval)
    #                 if cut_pos==43529026:
    #                     print('after')
    #                     print(ont_pop_block[49889])
    #                     print(ont_pop_block[49900])
    #                     break

    return flip_list, cut_dic             



  
def improve_vcf_flip(phased_vcf_dic, flip_list):

    phased_vcf_dic_improved=phased_vcf_dic
    for var_idx in flip_list:
        line_parts=phased_vcf_dic[var_idx]
        gt_flags, sample1=line_parts[8:10]
        sample1_split=sample1.split(":")
        gt_id_flags=gt_flags.split(":").index("GT")        
        allele_sample1=sample1_split[gt_id_flags]   

        allele_sample1_flp=str(1-int(allele_sample1[0]))+'|'+str(1-int(allele_sample1[2]))
        sample1_split_flp=sample1_split
        sample1_split_flp[gt_id_flags]=allele_sample1_flp
        line_parts_flp=line_parts
        line_parts_flp[9]=':'.join(sample1_split_flp)
        phased_vcf_dic_improved[var_idx]=line_parts_flp
    
    return phased_vcf_dic_improved  
    
    
    
    

def improve_vcf_cut(phased_vcf_dic_improved, hap_blocks, cut_dic):

# creating new block id dic

    new_blockid_dic={}

    for block_id in hap_blocks.keys():

        cut_block=cut_dic[block_id]
        idc_block=idc_blocks[block_id]

        cut_block_complete=[idc_block[0]]+cut_block+[idc_block[-1]] # including start and end var_idx of block

        for i in range(1,len(cut_block_complete)):
            pre_cut_idx=cut_block_complete[i-1]
            next_cut_idx=cut_block_complete[i]

            blockid=var_pos_list[pre_cut_idx-1]
            for k in range(pre_cut_idx,next_cut_idx+1):
                new_blockid_dic[k]=blockid


    phased_vcf_dic_improved2=phased_vcf_dic_improved

    sorted_var_idx=sorted(phased_vcf_dic_improved.keys())

    cut_dic_remain=cut_dic

    for var_idx in sorted_var_idx: # imprvd_phased_vcf_dic

        #if var_idx in dic1: # those variant that we need to change their block id 
        line_parts=phased_vcf_dic_improved[var_idx]
        gt_flags, sample1=line_parts[8:10]
        sample1_split=sample1.split(":")
        ps_id_flags=gt_flags.split(":").index("PS")
        block_id=int(sample1_split[ps_id_flags])

        new_blockid=str(new_blockid_dic[var_idx])


        sample1_split_changed=sample1.split(":")
        sample1_split_changed[ps_id_flags]=new_blockid
        line_parts_changed=line_parts
        line_parts_changed[9]=':'.join(sample1_split_changed)
        phased_vcf_dic_improved2[var_idx]=line_parts_changed

    return phased_vcf_dic_improved2






def write_out_vcf(phased_vcf_dic_improved2, vcf_file_improved_address):


    # reading phased VCF file,  use grep before that, only 0|1 or 1|0. bi -allelic
    vcf_file_improved = open(vcf_file_improved_address,'w');  # phased_vcf_dic

    sorted_var_idx=sorted(phased_vcf_dic_improved2.keys())

    for header_line in header_lines_list:
        vcf_file_improved.write(header_line+'\n')
    for var_idx in sorted_var_idx: # imprvd_phased_vcf_dic
        line_parts=phased_vcf_dic_improved2[var_idx]
        vcf_file_improved.write('\t'.join(line_parts)+'\n')
    vcf_file_improved.close()
    
    return 1 













if __name__ == "__main__":

    

    """
    
    Input:
    
    Output:
    
    """


    chrom_output=19
    vcf_file_address = 'data/'+str(chrom_output)+'/out10k.vcf' #tst.vcf' #

    hap_blocks, idc_blocks, var_pos_list, phased_vcf_dic, header_lines_list = read_vcf_file(vcf_file_address)

    
    

    
    file_mismatches_address='data/'+str(chrom_output)+'/'+str(chrom_output)+'_report_mismatches.txt'

    ont_pop_blocks=read_file_mismatches(file_mismatches_address, var_pos_list)

    flip_list, cut_dic = decide_flip_cut(hap_blocks, idc_blocks, var_pos_list, ont_pop_blocks)
#     phased_vcf_dic_improved = improve_vcf_flip(phased_vcf_dic, flip_list)
#     phased_vcf_dic_improved2 = improve_vcf_cut(phased_vcf_dic_improved, hap_blocks, cut_dic)
#     vcf_file_improved_address = 'data/'+str(chrom_output)+'/out_improved_'+str(tresh_match_mismatch)+'.vcf' #tst.vcf' #
#     write_out_vcf(phased_vcf_dic_improved2, vcf_file_improved_address)


