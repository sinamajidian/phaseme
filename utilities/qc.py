




if __name__ == "__main__":

    

    """
    
    Input:
    
    Output:
    
    """


    chrom_output=19
    vcf_file_address = 'data/'+str(chrom_output)+'/out1k.vcf' #tst.vcf' #

    hap_blocks, idc_blocks, var_pos_list, phased_vcf_dic, header_lines_list = read_vcf_file(vcf_file_address)

    
    # reading short haplotype 1kblock file
    file_hap_address='data/'+str(chrom_output)+'/pairs.txt'


    pop_inf_dic=read_file_hap(file_hap_address)



    print(len(pop_inf_dic))


    print('number of blocks in ont',len(hap_blocks))



    # comparing input phased vcf (hap_blocks_sample1_dic)with pop (short_blocks_dic_pos)

    ont_pop_blocks={}
    for block_id1, hap_block in hap_blocks.items(): # all ONT hap blocks  
        idc_block=idc_blocks[block_id1]

        #if block_id1==27731803: #27731803: 60780 20588536 8889839
        ont_pop_block=compare_phase_block_pop(hap_block,idc_block,pop_inf_dic,range(len(hap_block)))
        ont_pop_blocks[block_id1]=ont_pop_block

        
        
    report_out_address='data/'+str(chrom_output)+'/'+str(chrom_output)+'_report_mismatches.txt' #
    qual_blocks=report_match(report_out_address)


    report_qc_address='data/'+str(chrom_output)+'/'+str(chrom_output)+'_report_qc.txt' #
    report_qc(report_qc_address, qual_blocks, hap_blocks)

  