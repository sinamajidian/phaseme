#!/usr/bin/env python3


import numpy as np

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
                block_id1=int(sample1_split[ps_id_flags])
                if block_id1 in hap_blocks:
                    hap_blocks[block_id1].append(allele_sample1) # append new variants to the existing phase block
                    idc_blocks[block_id1].append(var_idx)
                else:
                    hap_blocks[block_id1]=[allele_sample1] # creat new phase block
                    idc_blocks[block_id1]=[var_idx]

            else:
                print('There is a vriant that ONT is not phased.. removed it first')

    return hap_blocks, idc_blocks, var_pos_list, phased_vcf_dic, header_lines_list





def read_file_hap(file_hap_address):

    file_hap = open(file_hap_address,'r'); 
    pop_inf_dic={} # key variant index, value two lists 
    for line in file_hap:
        line_parts=line.strip().split('\t')     # ['42081', '0', '42096', '0']
        snp1_idx=int(line_parts[0])
        allele1_snp1=int(line_parts[1])
        snp2_idx=int(line_parts[2])
        allele1_snp2=int(line_parts[3])
        #pairs_list.append([snp1_idx,allele1_snp1,snp2_idx,allele1_snp2])

        snps_idx=[snp1_idx,snp2_idx]
        snps_allele=[allele1_snp1,allele1_snp2]
        host_idx=max(snps_idx)
        guest_idx=min(snps_idx)
        if host_idx not in pop_inf_dic.keys():
            pop_inf_dic[host_idx]={'sim':list(),'dis':list()}
        if allele1_snp1==allele1_snp2:
            pop_inf_dic[host_idx]['sim'].append(guest_idx)  # sim shows the relation between two elements of a pop pair
        else:
            pop_inf_dic[host_idx]['dis'].append(guest_idx) 
            
    return pop_inf_dic




# comparing one phase block with all pop information

# block_id1=8889839 # test for chr 19
# hap_block=hap_blocks[block_id1]
# idc_block=idc_blocks[block_id1]

# if 1:
def compare_phase_block_pop(hap_block, idc_block, pop_inf_dic, interval):
    
    # input: interval is based on var_i 
    # output: is ont_pop_block[var_idx]=ont_pop 
    
    #var_i is the index of variant within block
    #var_idx is the index of vriant globally in the VCF file     
    
    ont_pop_block={}
    min_interval=min(list(interval))
    max_interval=max(list(interval))
    interval_idx=set(idc_block[min_interval:max_interval]) # set is faster than list for check if element is in it
    for var_i in interval: # the ONT hap block
        alleles=hap_block[var_i]
        var_idx=idc_block[var_i]
        
        #var_i is the index of variant within block
        #var_idx is the index of vriant globally in the VCF file 
        
        
        if var_idx in pop_inf_dic.keys():
            list_pairs=pop_inf_dic[var_idx]
            match_sim_list=[]
            match_dis_list=[]
            mismatch_sim_list=[]
            mismatch_dis_list=[]
            for sim_idx in pop_inf_dic[var_idx]['sim']: # sim shows the relation between two elements of a pop pair
                
                #if sim_idx in idc_block: # if the SNP  of pairs is in this ONT block
                try:
                    alleles_pop=hap_block[idc_block.index(sim_idx)] # extract the allele of the snp 'sim' which is the other snp in pair
                    if sim_idx in interval_idx:
                        if alleles_pop==alleles: 
                            match_sim_list.append(str(sim_idx))  
                        else:
                            mismatch_sim_list.append(str(sim_idx))     
                except:
                    pass
                        
            for sim_idx in pop_inf_dic[var_idx]['dis']: # sim shows the relation between two elements of a pop pair
                #if sim_idx in idc_block: # if the SNP  of pairs is in this ONT block
                try:
                    alleles_pop=hap_block[idc_block.index(sim_idx)] # extract the allele of the snp 'sim' which is the other snp in pair
                    if sim_idx in interval_idx:
                        if alleles_pop==str(1-int(alleles[0]))+'|'+str(1-int(alleles[2])):  
                            match_sim_list.append(str(sim_idx))  
                        else:
                            mismatch_sim_list.append(str(sim_idx))            
                except:
                    pass
        
            match_list=match_sim_list+match_dis_list
            mismatch_list=mismatch_sim_list+mismatch_dis_list
            ont_pop=[match_list,mismatch_list]
        else:
            ont_pop=[[],[]]

        ont_pop_block[var_idx]=ont_pop
    #print('comparison is done')
    return  ont_pop_block
    

    
    
#aa=0


#block_id1=27731803 # 60780  # 20588536 27731803 8889839
#hap_block=hap_blocks[block_id1]

#if 1:

def decide_flip_cut(hap_blocks, idc_blocks, var_pos_list):



    flip_list=[]
    cut_dic={}



    tresh_match_mismatch= 0.5 


    for block_id1, hap_block in hap_blocks.items():
        cut_dic[block_id1]=[]

        idc_block=idc_blocks[block_id1]
        ont_pop_block=ont_pop_blocks[block_id1]

        for var_i, var_idx in enumerate(idc_block): 
            #var_i is the index of variant within block
            #var_idx is the index of vriant globally in the VCF file 

            ont_pop = ont_pop_block[var_idx]  # after a cut, it may use a newr version of this vsrisble
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
                    print('cut candidate', cut_pos,var_i,var_idx)   
                    cut_dic[block_id1].append(var_idx)
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

    for block_id1 in hap_blocks.keys():

        cut_block=cut_dic[block_id1]
        idc_block=idc_blocks[block_id1]

        cut_block_complete=[idc_block[0]]+cut_block+[idc_block[-1]] # including start and end var_idx of block

        for i in range(1,len(cut_block_complete)):
            pre_cut_idx=cut_block_complete[i-1]
            next_cut_idx=cut_block_complete[i]

            blockid=var_pos_list[pre_cut_idx-1]
            for k in range(pre_cut_idx,next_cut_idx+1):
                new_blockid_dic[k]=blockid


    phased_vcf_dic_improved2=phased_vcf_dic_improved

    sorted_var_indx=sorted(phased_vcf_dic_improved.keys())

    cut_dic_remain=cut_dic

    for var_indx in sorted_var_indx: # imprvd_phased_vcf_dic

        #if var_indx in dic1: # those variant that we need to change their block id 
        line_parts=phased_vcf_dic_improved[var_indx]
        gt_flags, sample1=line_parts[8:10]
        sample1_split=sample1.split(":")
        ps_id_flags=gt_flags.split(":").index("PS")
        block_id1=int(sample1_split[ps_id_flags])

        new_blockid=str(new_blockid_dic[var_indx])


        #cut_dic[block_id1]
        sample1_split_changed=sample1.split(":")
        sample1_split_changed[ps_id_flags]=new_blockid
        line_parts_changed=line_parts
        line_parts_changed[9]=':'.join(sample1_split_changed)
        phased_vcf_dic_improved2[var_indx]=line_parts_changed

    return phased_vcf_dic_improved2






def write_out_vcf(phased_vcf_dic_improved2, vcf_file_improved_address):


    # reading phased VCF file,  use grep before that, only 0|1 or 1|0. bi -allelic
    vcf_file_improved = open(vcf_file_improved_address,'w');  # phased_vcf_dic

    sorted_var_indx=sorted(phased_vcf_dic_improved2.keys())

    for header_line in header_lines_list:
        vcf_file_improved.write(header_line+'\n')
    for var_indx in sorted_var_indx: # imprvd_phased_vcf_dic
        line_parts=phased_vcf_dic_improved2[var_indx]
        vcf_file_improved.write('\t'.join(line_parts)+'\n')
    vcf_file_improved.close()
    
    return 1 















if __name__ == "__main__":

    

    """
    
    Input:
    
    output:
    
    """


    chrom_output=19
    vcf_file_address = 'data/'+str(chrom_output)+'/out.vcf' #tst.vcf' #

    hap_blocks, idc_blocks, var_pos_list, phased_vcf_dic, header_lines_list = read_vcf_file(vcf_file_address)


    #### Parsing pusod read
    # handling pop information reporting only once

    # reading short haplotype 1kblock file
    file_hap_address='data/'+str(chrom_output)+'/pairs.txt'


    pop_inf_dic=read_file_hap(file_hap_address)


    # line_parts
    # print(line_parts,'\n', pairs_list[0])
    # print(len(pairs_list))
    # print(pop_inf_dic[host_idx])
    print(len(pop_inf_dic))


    print('number of blocks in ont',len(hap_blocks))

    # print('ont')
    # for block_id, pos_block in pos_blocks.items():
    #     print(block_id,len(pos_block))
    # print('ont')
    # for block_id1, hap_block  in hap_blocks.items():
    #     print(block_id1,len(hap_block))
    # hap_block[:2]



    # comparing input phased vcf (hap_blocks_sample1_dic)with pop (short_blocks_dic_pos)

    ont_pop_blocks={}
    for block_id1, hap_block in hap_blocks.items(): # all ONT hap blocks  
        idc_block=idc_blocks[block_id1]

        #if block_id1==27731803: #27731803: 60780 20588536 8889839
        ont_pop_block=compare_phase_block_pop(hap_block,idc_block,pop_inf_dic,range(len(hap_block)))
        ont_pop_blocks[block_id1]=ont_pop_block

    #hap_blocks.keys()
    #var_pos_list[31566-1] #grep -v "#" out.vcf | sed -n 31566,31566p
    #ont_pop_blocks[60780][31576]
    # print(ont_pop_blocks[60780][292])
    # print('for variant',var_pos_list[292-1],'we have',var_pos_list[268-1],var_pos_list[269-1],var_pos_list[261-1])

    flip_list, cut_dic = decide_flip_cut(hap_blocks, idc_blocks, var_pos_list)

    phased_vcf_dic_improved = improve_vcf_flip(phased_vcf_dic, flip_list)

    phased_vcf_dic_improved2 = improve_vcf_cut(phased_vcf_dic_improved, hap_blocks, cut_dic)


    vcf_file_improved_address = 'data/'+str(chrom_output)+'/out_improved.vcf' #tst.vcf' #



    write_out_vcf(phased_vcf_dic_improved2, vcf_file_improved_address)

