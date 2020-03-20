#### Parsing phased VCF file containing two sample for comparison


# vcf_file_address = "/Volumes/uni/working_phasme/from_ssm/19_nei10/19_ont_true.vcf" # 19_ont_true.vcf" 19_ont_true_improved_.90
from sys import argv
import numpy as np


#### Parsing phased VCF file containing two sample for comparison


vcf_file_address =  argv[1] #"/home/ssm/Documents/phaseme/ont/1/1_ont_true.vcf"
tr = int(argv[2])

# "/Volumes/uni/myjupyter/jupyter_phaseme/ont_son/1/1_ont_true.vcf" # son_improved_.90_tr1


vcf_file = open(vcf_file_address,'r');

hap_blocks_sample1_dic={} # key: id of phase block (i.e. phase st), value= {genomic_position:allele_1}
hap_blocks_sample2_dic={}

hap_len_sample1_dic={}

header_lines_list=[]
for line in vcf_file:
    line_strip=line.strip()
    if line_strip.startswith('#'):
        header_lines_list.append(line_strip)
        sample_names=line_strip.split('\t')[9:11]# last line of header contains sample name
    else:
        #var_indx+=1

        # '22\t51244182\t.\tG\tA\t.\t.\t.\tGT:GQ:DP:AF:GL:PS\t
        #  0|1:127:.:.:-58.12751217387621,-8.336610948353968e-14,-12.716893473415803:50415657'
        line_parts=line_strip.split('\t')

        gt_flags, sample1, sample2 = line_parts[8:11]
        if "|" in sample1:

            var_pos=int(line_parts[1])

            sample1_split=sample1.split(":")

            gt_flags_split=gt_flags.split(":")

            block_id1 = sample1_split[gt_flags.split(":").index("PS")] #int(
            allele_sample1=sample1_split[gt_flags.split(":").index("GT")]
            if allele_sample1 == '0|1' or '1|0':
                if block_id1 in hap_len_sample1_dic:
                    hap_len_sample1_dic[block_id1].append(var_pos) # append new variants to the existing phase block
                else:
                    hap_len_sample1_dic[block_id1]= [var_pos] # creat new phase block


            if "|" in sample2: # both sample are phased for this variant




                sample2_split=sample2.split(":")


                block_id2 = sample2_split[gt_flags.split(":").index("PS")] #int(

                allele_sample2=sample2_split[gt_flags.split(":").index("GT")]

                if (allele_sample1 == '0|1' or allele_sample1 == '1|0') and (allele_sample2 == '0|1' or allele_sample2 == '1|0'):
                    if block_id1 in hap_blocks_sample1_dic:
                        hap_blocks_sample1_dic[block_id1][var_pos]=allele_sample1 # append new variants to the existing phase block
                    else:
                        hap_blocks_sample1_dic[block_id1]={var_pos:allele_sample1} # creat new phase block

                    if block_id2 in hap_blocks_sample2_dic:
                        hap_blocks_sample2_dic[block_id2][var_pos]=allele_sample2 # append new variants to the existing phase block
                    else:
                        hap_blocks_sample2_dic[block_id2]={var_pos:allele_sample2} # creat new phase block


print('Number of blocks in sample:',len(hap_blocks_sample1_dic))



parental_origin_blocks = {}
varpos_blocks={}



for block_id_sample1, hap_block_sample1 in hap_blocks_sample1_dic.items():

    parental_origin_block ={}

    varpos_block={}


    for block_id_sample2, hap_block_sample2 in hap_blocks_sample2_dic.items():
        parental_origin_block_shared = []
        varpos_block_shared=[]

        var_pos_list = sorted(list(hap_block_sample1.keys()))
        for var_i, var_pos in enumerate(var_pos_list):

            allele_sample1 = hap_block_sample1[var_pos]
            if var_pos in hap_block_sample2.keys():

                allele_sample2 = hap_block_sample2[var_pos]
                allele_sample2_revert = str(1-int(allele_sample2[0]))+'|'+str(1-int(allele_sample2[2]))


                if allele_sample1 == allele_sample2:
                    parental_origin = 1 #zero_one  # considering the sample1 block constant, for each grand truth, we need to change this value

                if allele_sample1 == allele_sample2_revert :
                    parental_origin = 0 #1- zero_one


                parental_origin_block_shared.append(parental_origin)
                varpos_block_shared.append(var_pos)
        if len(parental_origin_block_shared)>= 20:
            parental_origin_block[block_id_sample2]=parental_origin_block_shared
            varpos_block[block_id_sample2]=varpos_block_shared


    if len(parental_origin_block)> 0:
        parental_origin_blocks[block_id_sample1] = parental_origin_block
        varpos_blocks[block_id_sample1] = varpos_block





list_all_consecutive= []
list_all_blocks ={}
all_switch_list=[]
for block_id_sample1, parental_origin_block in parental_origin_blocks.items():

    list_all_block ={}
    varpos_block=varpos_blocks[block_id_sample1]
    for block_id_sample2, parental_origin_block_shared in parental_origin_block.items():
        varpos_block_shared=varpos_block[block_id_sample2]


        new_list=[0]

        for i in range(len(parental_origin_block_shared)):

            if i >=1:
                if parental_origin_block_shared[i] != parental_origin_block_shared[i-1]:
                    new_list.append(i)  #  consists of the starting position of a var that its parental origin (comapred to true) is different than previous ar
                    all_switch_list.append(varpos_block_shared[i])
        new_list.append(len(parental_origin_block_shared))
        list_all_consecutive.append(new_list)
        list_all_block[block_id_sample2]=new_list
    list_all_blocks[block_id_sample1] = list_all_block






list_length_minor = [] # flipped or same

for list1 in list_all_consecutive: #  list1  consists of the starting position of a var that its parental origin (comapred to true) is different than previous ar

    list_length = []
    for i in range(1,len(list1)):
        length1= list1[i]-list1[i-1]
        list_length.append(length1)
        #list_length_all_consecutive.append(length1)
        #if i!=1: list_length_all_consecutive_except_first.append(length1)


    odd_sum = sum(list_length[0::2])
    even_sum = sum(list_length[1::2])
    min_sum= min([odd_sum, even_sum] )
    val = [odd_sum, even_sum].index(min_sum)
    length_minor=list_length[val::2]
    list_length_minor +=length_minor

list_length_minor= np.array(list_length_minor)

# tr=50
lower_than_21=[length1<tr+1 for length1 in list_length_minor]
greater_than_20=[length1>tr for length1 in list_length_minor]

print("smaller than",tr+1)
print("smaller ",sum(lower_than_21))
print("greater than",tr)
print("greater ",sum(greater_than_20))





# length_minor_blocks ={}
# firstIdx_minor_blocks ={}
# for block_id_sample1, list_all_block in list_all_blocks.items():

#     length_minor_block = {}
#     firstIdx_minor_block = {}

#     for block_id_sample2, list1 in list_all_block.items():

#         list_length = []
#         for i in range(1,len(list1)):
#             length1= list1[i]-list1[i-1]
#             list_length.append(length1)

#         odd_sum = sum(list_length[0::2])
#         even_sum = sum(list_length[1::2])
#         min_sum= min([odd_sum, even_sum] )
#         val = [odd_sum, even_sum].index(min_sum)
#         length_minor=list_length[val::2]


#         firstIdx_minor = list1[:-1][val::2]       # switches
#         if len(length_minor):
#             length_minor_block[block_id_sample2]=length_minor
#             firstIdx_minor_block[block_id_sample2]=firstIdx_minor

#     if len(length_minor_block):
#         length_minor_blocks[block_id_sample1]=length_minor_block
#         firstIdx_minor_blocks[block_id_sample1]=firstIdx_minor_block



# first_het_variant = True

# vcf_file = open(vcf_file_address,'r')

# var_pos_blocks = [] # position of those variants that the truth is known

# id_blocks=[]
# for line in vcf_file:

#     line_strip = line.strip()

#     if line_strip.startswith('#'):
#         pass
#         #header_lines_list.append(line_strip)
#         #sample_names = line_strip.split('\t')[9:11]            # last line of header contains sample name
#     else:

#         line_parts = line_strip.split('\t')

#         chrom = line_parts[0]
#         var_pos = int(line_parts[1])                           # genomic position of variants

#         format_genotype, values_genotype = line_parts[8:10]    # 'GT:GQ:DP:AF:GL:PS', '0|1:255:.:.:.,0,.:60780'

#         values_genotype_splitted = values_genotype.split(':')
#         format_genotype_splitted = format_genotype.split(':')

#         gt_index = format_genotype_splitted.index("GT")           #  index of allele in  values_genotype

#         allele = values_genotype_splitted[gt_index]

#         if (allele == '0|1' or allele == '1|0'):

#             ps_index = format_genotype_splitted.index("PS")           #  index of phase set in values_genotype
#             id_block = values_genotype_splitted[ps_index]

#             values_genotype_truth = line_parts[10]
#             values_genotype_truth_splitted = values_genotype_truth.split(':')
#             allele_truth = values_genotype_truth_splitted[gt_index]
#             id_block_truth = values_genotype_truth_splitted[ps_index]

#             #if allele_truth != './.':
#             if allele_truth == '0|1'  or allele_truth == '1|0':

#                 if first_het_variant:           # for the first het variant
#                     first_het_variant = False

#                     var_pos_block = [int(var_pos)]
#                     id_blocks.append(id_block)

#                 else:                              # for the rest of het variants
#                     if id_block in id_blocks:
#                         var_pos_block.append(int(var_pos))

#                     else:

#                         # add previous block to the list of all blocks
#                         var_pos_blocks.append(var_pos_block)

#                         # creat new phase block
#                         var_pos_block = [int(var_pos)]
#                         id_blocks.append(id_block)

# # # for the last het variant, we  finish the last block.
# var_pos_blocks.append(var_pos_block)




# for block_id_sample1, length_minor_block in length_minor_blocks.items():
#     firstIdx_minor_block = firstIdx_minor_blocks[block_id_sample1]


#     for block_id_sample2, length_minor in length_minor_block.items():

#         firstIdx_minor = np.array(firstIdx_minor_block[block_id_sample2])

#         list_length_minor = np.array(length_minor)

#         lower_than_21_logic=[length1<tr for length1 in list_length_minor]

#         greater_than_21_logic=[not i for i in lower_than_21_logic] # greater than and equal

#         print(block_id_sample1,sum(lower_than_21_logic),sum(greater_than_21_logic))
#         if sum(greater_than_21_logic):
#             var_pos_list = var_pos_blocks[id_blocks.index(block_id_sample1)]
#             var_pos_switch_list =  [ var_pos_list[i] for i in firstIdx_minor[greater_than_21_logic]]

#             print("The switch starts at",var_pos_switch_list)
#             # grep "20610824" 22_report_mismatches_both.txt | grep -v "^\." | sed -n 10783,10785p
#             print("with the length", list_length_minor[greater_than_21_logic])

list_len=[len(pos_list) for block_id, pos_list in hap_len_sample1_dic.items()]
list_len_kb=[pos_list[-1]-pos_list[0] for block_id, pos_list in hap_len_sample1_dic.items() if len(pos_list)>1]

# print(list_len)
print("length" ,round(np.mean(list_len_kb)/1000,2))





print(all_switch_list)
