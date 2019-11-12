#!/usr/bin/env python3


from sys import argv



if __name__ == "__main__":

    
    vcf_file_address = argv[1]
    # "/Volumes/uni/working_phasme/data_bcm/grand_truth/new_trio/1_trio.vcf"


    vcf_file = open(vcf_file_address,'r')

    lines_new =[]
    header_lines_list = []

    num=0
    for line in vcf_file:
        line_strip = line.strip() 

        if line_strip.startswith('#'):
            pass 
            header_lines_list.append(line_strip)

        else:

            line_parts = line_strip.split('\t') 
            # chrom = line_parts[0]
            var_pos = int(line_parts[1])                               # genomic position of variants
            ref_allele = line_parts[3]
            alt_allele = line_parts[4]


            format_genotype = line_parts[8]                          #  'GT:GQ:DP:AF:GL:PS', 
            format_genotype_splitted = format_genotype.split(':')


            values_genotype_sample1 = line_parts[9]                         # first sample '0|1:255:.:.:.,0,.:60780'

            values_genotype_sample1_splitted = values_genotype_sample1.split(':')


            gt_index = format_genotype_splitted.index("GT")            #  index of allele in  values_genotype 
            allele = values_genotype_sample1_splitted[gt_index]

            values_genotype_hg3 = line_parts[10]
            values_genotype_hg3_splitted = values_genotype_hg3.split(':')
            allele_hg3 = values_genotype_hg3_splitted[gt_index]

            values_genotype_hg4 = line_parts[11]
            values_genotype_hg4_splitted = values_genotype_hg4.split(':')
            allele_hg4 = values_genotype_hg4_splitted[gt_index]


            phase_son = './.'
            #if '.' not in allele and len(ref_allele)==1 and len(alt_allele)==1 : #(allele == '0|1' or allele == '1|0'):
            if (allele == '0|1' or allele == '1|0' or allele == '1/0' or allele == '0/1') and len(ref_allele)==1 and len(alt_allele)==1:


                #ps_index = format_genotype_splitted.index("PS")           #  index of phase set in values_genotype 


                if  allele_hg3 != './.' and (not allele_hg3 == '0/0') and (not allele_hg3 == '0|0') and allele_hg4 == './.' : 
                    phase_son = '0|1'

                if  allele_hg4 != './.' and (not allele_hg4 == '0/0') and (not allele_hg4 == '0|0') and allele_hg3 == './.' :
                    phase_son = '1|0'


            if phase_son == '0|1' or phase_son == '1|0':
                values_genotype_sample1_splitted[gt_index] = phase_son
                num +=1

                line_parts[2] = '.'
                line_parts[5] = '.'
                line_parts[6] = '.'
                line_parts[7] = '.'
                line_parts[8] = "GT:PS"#values_genotype
                line_parts[9] =  phase_son+":ps" # values_genotype

                line_new = '\t'.join(line_parts[:10]) 
                lines_new.append(line_new)

    vcf_file.close()

    print("Number of variant in grand truth ",num)

    vcf_file_write = open(vcf_file_address[:-4]+"_true.vcf",'w')

    header_lines_list[-1]= header_lines_list[-1][:51]
    for header_line in header_lines_list[:-1]:
        vcf_file_write.write(header_line+'\n')

    vcf_file_write.write("##FORMAT=<ID=PS,Number=1,Type=String,Description=\"Phase set in which this variant falls\">\n")
    vcf_file_write.write(header_lines_list[-1]+'\n')


    for line_new in lines_new:
        vcf_file_write.write(line_new+'\n')

    vcf_file_write.close()

