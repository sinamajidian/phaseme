#!/usr/bin/python3

def split_vcf(vcf_file, out_prefix):

    """
    Input: a tab delimited vcf file

    output: several vcf files in out_prefix folder. Every of them is for only one chromosome


    In this version, the first column of the VCF file should be in this format  22.


    """


    subprocess.call("mkdir "+out_prefix , shell=True)
    subprocess.call("cp "+vcf_file+" "+out_prefix+"/input.vcf", shell=True)



    chrs_list = []
    for i in range(1,23):
        i=str(i)

        extract_chri_bash = "grep -P \"^"+i+"\t\" "+out_prefix+"/input.vcf > "+out_prefix+"/temp"+i
        subprocess.call(extract_chri_bash, shell=True)

        length_raw = check_output("wc -l "+out_prefix+"/temp"+i,shell=True) #","wc","-l""

        length_raw2 = str(length_raw)
        length_str = length_raw2[2:].strip().split(' ')[0]

        if length_str:
            length =int(length_str)

            if length>1:
                chrs_list.append(i)
                chr_based_vcf_bash=["mkdir "+out_prefix+"/"+i+";",
                            "grep \"#\" "+out_prefix+"/input.vcf > "+out_prefix+"/"+i+"/chr"+i+".vcf;",
                            "cat "+out_prefix+"/temp"+i+" >> "+out_prefix+"/"+i+"/chr"+i+".vcf;"]

                subprocess.run("".join(chr_based_vcf_bash), shell=True)

        subprocess.call("rm -rf "+out_prefix+"/temp"+i, shell=True)

    return chrs_list





def run_shapeit_graph(shapeit_address, data_1000G_address,out_prefix ,chrom): # run per chromosom

    d1kg=data_1000G_address
    shapeit_check= [shapeit_address,"/shapeit -check",
                    " --input-vcf "+out_prefix+"/"+chrom+"/chr",chrom,".vcf",
                    " -R ",d1kg,"1000GP_Phase3_chr",chrom,".hap.gz ",d1kg,"1000GP_Phase3_chr",chrom,".legend.gz ",d1kg,"1000GP_Phase3.sample",
                    " --output-log "+out_prefix+"/",chrom,"/shapeit_check",
                    " >> "+out_prefix+"/",chrom,"/shapeit_check_out.log"]

    shapeit_check_bash="".join(shapeit_check)
    #print(shapeit_check_bash)
    subprocess.call(shapeit_check_bash, shell=True)

    shapeit_graph= [shapeit_address,"/shapeit ",
                    " --input-vcf  "+out_prefix+"/"+chrom+"/chr",chrom,".vcf",
                    " -R "+d1kg+"1000GP_Phase3_chr"+chrom+".hap.gz "+d1kg+"1000GP_Phase3_chr"+chrom+".legend.gz "+d1kg+"1000GP_Phase3.sample ",
                    " -M "+d1kg+"genetic_map_chr"+chrom+"_combined_b37.txt",
                    " --output-log "+out_prefix+"/",chrom+"/shapeit_graph",
                    " --output-graph "+out_prefix+"/",chrom,"/chr",chrom,".graph"]

    exclude_file= out_prefix+"/"+chrom+"/shapeit_check.snp.strand.exclude"
    if  path.exists(exclude_file):

        shapeit_graph=shapeit_graph+[" --exclude-snp ",exclude_file,
                    " >> "+out_prefix+"/",chrom,"/shapeit_graph_strand_out.log" ]

    else:
        shapeit_graph=shapeit_graph+[" >> "+out_prefix+"/",chrom,"/shapeit_graph_out.log" ]

    shapeit_graph_bash="".join(shapeit_graph)
    #print(shapeit_graph_bash)

    subprocess.call(shapeit_graph_bash, shell=True)

    return 1


def sample_haplotype_graph(input_graph, num_samples,chrom):

    subprocess.call("mkdir "+out_prefix+"/"+chrom+"/samples", shell=True)

    for sample_i in range(num_samples):

        sample_i=str(sample_i)
        seed = str(random.randint(1,10000000))


        sample_graph= [shapeit_address,"/shapeit  -convert --seed ",seed,
                       " --input-graph ",input_graph,
                       " --output-sample "+out_prefix+"/",chrom,"/samples/sample_"+sample_i,
                       " -L "+out_prefix+"/",chrom,"/samples/sample_",sample_i,
                       " >> "+out_prefix+"/",chrom,"/samples.log"]

        sample_graph_bash = "".join(sample_graph)
        #print(sample_graph_bash)
        subprocess.call(sample_graph_bash, shell=True)

    return 1


def read_haplotype_samples(haplotype_sample_address, num_samples):
    haplotype1_samples = []

    for sample_i in range(num_samples):

        haplotype_sample_address_file =  haplotype_sample_address+"/sample_"+str(sample_i)+".haps" #"files_"+str(chrom)+"/samples/sample_"+str(sample_i)+".haps"

        var_pos_list, haplotype1  = read_haplotype_sample(haplotype_sample_address_file)

        haplotype1_samples.append(haplotype1)

        if sample_i>0:
            if var_pos_list_pre!=var_pos_list:
                print("inconsistency in shapeit output, sample index",sample_i,". Please run it again")
                exit(1)

        var_pos_list_pre=var_pos_list


    return haplotype1_samples, var_pos_list



def read_haplotype_sample(haplotype_sample_address_file):

    sample_file = open(haplotype_sample_address_file,'r')

    var_pos_list = []   # genomic position of variants
    haplotype1 = []  # hetrozygous variant only

    for line in sample_file:

        line_splitted = line.strip().split();

        chrom = line_splitted[0]
        var_pos = line_splitted[2]

        allele1 = line_splitted[5]
        allele2 = line_splitted[6]

        if allele1 != allele2:
            haplotype1.append(allele1)
            var_pos_list.append(var_pos)

    return var_pos_list, haplotype1




def pairwise(haplotype1_samples, nsamples, i, j):

    identical_phasing = 0
    opposite_phasing = 0

    for sample_i in range(nsamples):

        if haplotype1_samples[sample_i][i] == '0' and haplotype1_samples[sample_i][j] == '0': identical_phasing +=1
        elif haplotype1_samples[sample_i][i] == '1' and haplotype1_samples[sample_i][j] == '1': identical_phasing +=1

        elif haplotype1_samples[sample_i][i] == '0' and haplotype1_samples[sample_i][j] == '1': opposite_phasing +=1
        elif haplotype1_samples[sample_i][i] == '1' and haplotype1_samples[sample_i][j] == '0': opposite_phasing +=1

    return identical_phasing, opposite_phasing



def extract_pairs(haplotype1_samples):

    num_samples =  len(haplotype1_samples)
    num_variants = len(haplotype1_samples[0])

    pairs =[];

    for i in range(num_variants):

        for j in range(i+1, min(num_variants, i+NEIGHBOURS)):

            identical_phasing, opposite_phasing= pairwise(haplotype1_samples, num_samples, i, j)

            f = float(identical_phasing+0.5)/(identical_phasing+opposite_phasing+1)

            if f > THRESH:
                pairs.append([i,j,'identical']);

            elif 1.0-f > THRESH:
                pairs.append([i,j,'opposite']);

    #print(i)
    pairs_sorted = sorted(pairs, key=lambda item: item[0])

    return pairs_sorted



def report_pairs(file_pair_address, pairs, var_pos_list):

    file_pairs= open(file_pair_address, 'w');  # output

    for pair in pairs:

        line_out = str(var_pos_list[pair[0]])+'\t'+str(var_pos_list[pair[1]])+'\t'+pair[2]+'\n'
        file_pairs.write(line_out)

    file_pairs.close()
    return 1






def pair_linkage(chrs_list, shapeit_address, data_1000G_address, num_samples, ):

    """

    Extracting indiviudual specific linakge information

        Previously population.py


    """

    for chrom in chrs_list:

        run_shapeit_graph(shapeit_address, data_1000G_address,out_prefix ,chrom)
        haplotype_graph = out_prefix+"/"+chrom+"/chr"+chrom+".graph"
        print("Haplotype graph is generated in "+haplotype_graph)

        sample_haplotype_graph(haplotype_graph, num_samples,chrom)


        print(str(num_samples)+" haplotype samples are generated from haplotype graph.")

        haplotype_sample_address= out_prefix+"/"+chrom+"/samples"  # /sample_"+str(sample_i)+".haps"

        haplotype1_samples, var_pos_list = read_haplotype_samples(haplotype_sample_address, num_samples)
        print("samples are read.")


        pairs = extract_pairs(haplotype1_samples)

        file_pair_address = out_prefix+"/"+chrom+"/"+chrom+"_pairs.txt" # _"+str(num_samples)+"_"+str(THRESH)+".txt"
        report_pairs(file_pair_address, pairs, var_pos_list)
        print(str(len(pairs))+" pairs are reported in "+file_pair_address)

    return 1





def read_vcf_file(vcf_file_address):

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
            stats_vcf = [homozygous0_num, homozygous1_num, hetrozygous_nonphased, hetrozygous_phased, genomic_length_blocks, n50,phase_rate]

    """



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
    first_first= True


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

            var_pos = int(line_parts[1])                           # genomic position of variants
            if first_first==True:
                var_pos_first=var_pos
                first_first=False

            format_genotype, values_genotype = line_parts[8:10]    # 'GT:GQ:DP:AF:GL:PS', '0|1:255:.:.:.,0,.:60780'

            values_genotype_splitted = values_genotype.split(':')
            format_genotype_splitted = format_genotype.split(':')

            gt_index = format_genotype_splitted.index("GT")           #  index of allele in  values_genotype

            allele = values_genotype_splitted[gt_index]


#             ## how should we handle '2' in allele ?
#              if './.' in allele:
#                  print("There is a vriant with genomic position "+str(var_pos)+" that is not genotyped. Remove it first.")
#                  exit(1)

#             ## if '/' in allele: print("There is a vriant with genomic position "+str(var_pos)+" that is not phased. Remove it first.")

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


    var_pos_last=var_pos
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

    phase_rate=np.sum(genomic_length_blocks)/ (var_pos_last-var_pos_first)


    stats_vcf = [homozygous0_num, homozygous1_num, hetrozygous_nonphased, hetrozygous_phased, genomic_length_blocks, n50,phase_rate]


    return lines_list, var_pos_het_list, line_number_het_list, id_blocks, allele_blocks, var_pos_blocks, stats_vcf, chrom


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




def report_comparison(report_out_address, comparison_result_blocks, id_blocks, chrom):

    file_report= open(report_out_address,'w');
    file_report.write("#chr\t var_pos\t blockid:alleles_ont \t mismatched_pair \t  matched_pair  \n")

    qual_blocks=[]

    for block_i, block_id  in enumerate(id_blocks):
        #print(block_i,block_id)



        allele_block = allele_blocks[block_i]
        var_pos_block = var_pos_blocks[block_i]


        file_report.write('\t'.join(["#Phase Block with ID  ",str(block_id)])+"\n")

        comparison_result_block = comparison_result_blocks[block_i]

        qual_block=[]


        for var_i in range(len(allele_block)): # var_i internal index of this phase block

            allele=allele_block[var_i]
            var_pos=var_pos_block[var_i]


            list_write=[]

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





def report_qc(report_qc_address, id_blocks, qual_blocks, allele_blocks, stats_vcf, chrom):

    [homozygous0_num, homozygous1_num, hetrozygous_nonphased, hetrozygous_phased, genomic_length_blocks, n50, phase_rate] = stats_vcf


    q_list=[]
    for block_i, block_id  in enumerate(id_blocks):
        qual_block=qual_blocks[block_i]

        if len(qual_block):
            q_list.append(round(np.mean(qual_block),5))



    file_report_qc= open(report_qc_address, 'w');


    first_line_list=["Chromosome","N50(Kb)","Avg. phase block quality","Number of phased heterozygous variants",
                     "Number of non-phased heterozygous variants","Number of homozygous variants","Phase rate"]


    file_report_qc.write("##"+",\t".join(first_line_list)+"\n\n")
    second_line_list=[str(chrom),str(n50),str(round(np.mean(q_list),5)),
                      str(hetrozygous_phased),str(hetrozygous_nonphased),
                      str(homozygous0_num+homozygous1_num),str(phase_rate)]
    file_report_qc.write(",\t".join(second_line_list)+"\n\n")

    file_report_qc.write("Quality report per block for chromosome "+str(chrom)+":\n")
    file_report_qc.write("##Block_i,\tStart_pos,\tblock_length,\tnum_phased_SNV,\t\tblock_quality,\n")




    for block_i, block_id  in enumerate(id_blocks):

        qual_block=qual_blocks[block_i]
        allele_block=allele_blocks[block_i]
        genomic_length_block = genomic_length_blocks[block_i]

        if len(qual_block):
            q = round(np.mean(qual_block),5)
        else:
            q= 'NaN'



        file_report_qc.write(str(block_i)+",\t"+block_id+",\t"+str(genomic_length_block)+",\t"+str(len(allele_block))+",\t"+str(q))

        file_report_qc.write("\n") # new  block

    file_report_qc.close()
    return 1








def decide_cut(id_blocks, allele_blocks, var_pos_blocks, comparison_result_blocks):

    cut_list_blocks = []

    for block_i, block_id  in enumerate(id_blocks):


        # print('working on block :',block_id)
        allele_block = allele_blocks[block_i]
        var_pos_block = var_pos_blocks[block_i]

        # after deciding the cut, it is not applied, so the number of block is not changed.
        # but after deciding the match/mismatch information is updated for those varaint afterwards and removed for previose since new block is started practically.

        comparison_result_block = comparison_result_blocks[block_i] # ont_pop_block


        cut_list_block=[]


        lower_bound = var_pos_block[0]
        upper_bound = var_pos_block[-1]

        for var_i, var_pos in enumerate(var_pos_block):  # var_pos_block idc_block

            # var_i is the index of variant within block
            # var_pos is the genomic position of variant



            comparison_result = comparison_result_block[var_i]  # rightafter we decide a cut, we use an updated version of comparison_result_block

            if comparison_result == [[],[]]:  # there is no population information for this variant
                out_pairs = '.'

            else:
                [matched_list_raw, mismatched_list_raw] = comparison_result
                matched_list=   [i for i in matched_list_raw if int(i)>=lower_bound and int(i)<=upper_bound]
                mismatched_list=[i for i in mismatched_list_raw if int(i)>=lower_bound and int(i)<=upper_bound]


                num_mismatched = len(mismatched_list)
                num_matched = len(matched_list)

                condition_cut = 0

                if num_mismatched >= 2:


                    try:
                        [matched_next1_raw, mismatched_next1_raw] = comparison_result_block[var_i+1]

                        matched_next1=   [i for i in matched_next1_raw if int(i)>=lower_bound and int(i)<=upper_bound]
                        mismatched_next1=[i for i in mismatched_next1_raw if int(i)>=lower_bound and int(i)<=upper_bound]

                    except IndexError: # Next variants are not in the block
                        matched_next1 = []
                        mismatched_next1 = []

                    num_matched_next1 = len(matched_next1)
                    num_mismatched_next1 = len(mismatched_next1)


                    for var_neighbour_j in range(var_i+1, min(len(var_pos_block),var_i+NEIGHBOURS+1)):

                        [matched_neighbour_j_raw, mismatched_neighbour_j_raw] = comparison_result_block[var_neighbour_j]
                        matched_neighbour_j=   [i for i in matched_neighbour_j_raw if int(i)>=lower_bound and int(i)<=upper_bound]
                        mismatched_neighbour_j=[i for i in mismatched_neighbour_j_raw if int(i)>=lower_bound and int(i)<=upper_bound]

                        num_matched_j = len(matched_neighbour_j)
                        num_mismatched_j = len(mismatched_neighbour_j)

                        if num_mismatched_j >=1 and (num_matched_j == 0 or num_matched_j/num_mismatched_j < tresh_match_mismatch_cut):
                            condition_cut = 1

                        # ideal case, mismatchin_in_next_ten should be before the  var_pos to be a good cut


                if condition_cut :#and not condition_flip:

                    cut_pos = var_pos                 # cut_pos the starting position of new block
                    #print('cut candidate', cut_pos)

                    cut_list_block.append(var_pos)

                    lower_bound = var_pos # var_pos_block[0] #
                    upper_bound = var_pos_block[-1]

        cut_list_blocks.append(cut_list_block)

    return cut_list_blocks




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





























from sys import argv
import numpy as np

import subprocess

from subprocess import check_output
import random


from os import path












if __name__ == "__main__":



    """


    # precomputed mode- only QC
        python phaseme.py my.vcf qc precomputed

    # precomputed mode- imrpoving
        python phaseme.py my.vcf improver precomputed

    # Complete usage- imrpoving
        python phaseme.py  my.vcf improver individual path_shapeit path_1000g

        please make sure by running the following in terminal
           ls path_shapeit/shapeit
           ./path_shapeit/shapeit
           ls path_1000g/1000GP_Phase3_chr1.hap.gz


    """

#    shapeit_address =  "/home/ssm/Documents/phaseme/"
#    data_1000G_address = "/home/ssm/Documents/phaseme/data1/1000g/"
#    mode_phasme = "precomputed"  # "precomputed" "individual"
#    mode_phasme_qc_improver = "improver"  # "improver"  "qc"


    help_note=[ "\n   python phaseme.py my.vcf out qc        precomputed",
            "   python phaseme.py my.vcf  out improver  precomputed",
            "   python phaseme.py my.vcf  out improver  individual path_shapeit path_1000g",
            "If you choose individual, please make sure that the followings work well in bash:  ls path_1000g/1000GP_Phase3_chr1.hap.gz; ls path_shapeit/shapeit; ./path_shapeit/shapeit;",
            "  "]


    if len(argv)==1:


        print("\n".join(help_note))
        exit()


    elif len(argv)<5:
        print("\nPlease provide enough argumnets \n \n")
        print("\n".join(help_note))
        exit()

    elif len(argv)==5:
        mode_phasme = "precomputed"

    elif len(argv)==6:
        mode_phasme = "individual"
        shapeit_address = argv[5]
        data_1000G_address = argv[6]

    mode_phasme_qc_improver = argv[3]
    out_prefix = argv[2]
    vcf_file_address = argv[1]




    chrs_list = split_vcf(vcf_file_address, out_prefix)        # # "In this version, the first column of the VCF file should be in this format  22  ."
    if not len(chrs_list):
        print("The input VCF file at "+vcf_file_address+"is empty or does not exist or the first column (showing the chromosomes) is not a sole number (It shouldn't be chr1!)")
        exit()
    print("Input VCF file contains "+str(len(chrs_list))+" chromosome(s) and is split into chr-based VCFs in "+out_prefix+".")


    NEIGHBOURS =  20     # number of neighbour variants to be checked
    tresh_match_mismatch_cut = 1

    if  mode_phasme == "individual":

        THRESH=0.90          # the extent of  between samples
        num_samples =  500   # the number of times that we sample the haplotype graph (output of shapeit)

        print("In case of error, please check the log files as well.")
        pair_linkage(chrs_list, shapeit_address, data_1000G_address, num_samples)



    elif mode_phasme == "precomputed":

        file_pairs_address_pre= "precomputed/"

        for chrom in chrs_list:

            copy_pre_pair_bash= "cp "+file_pairs_address_pre+"/pair_"+chrom+".txt "+out_prefix+"/"+chrom+"/"+chrom+"_pairs.txt"
            subprocess.call(copy_pre_pair_bash, shell=True)



    else:

        print("Please specify the mode of phaseME correctly: either individual or precomputed")

        # 22_ont_pairs_500_0.9.txt"

    if mode_phasme_qc_improver == "qc":


        #### QC part ######

        for chrom in chrs_list:

            print("working on chromosome ", chrom)

            vcf_file_address = out_prefix+"/"+chrom+"/chr"+chrom+".vcf"
            lines_list, var_pos_het_list, line_number_het_list, id_blocks, allele_blocks, var_pos_blocks, stats_vcf, chrom = read_vcf_file(vcf_file_address)
            #print(len(var_pos_het_list))

            file_pairs_address = out_prefix+"/"+chrom+"/"+chrom+"_pairs.txt"

            pop_inf_dic = read_file_pairs_forward(file_pairs_address)

            out_name_prefix= vcf_file_address[:-4] #argv[3]  #

            # comparing input phased vcf (hap_blocks_sample1_dic) with pop (short_blocks_dic_pos)
            comparison_result_blocks=[]

            for block_i in range(len(id_blocks)):

                allele_block = allele_blocks[block_i]
                var_pos_block = var_pos_blocks[block_i]

                # lower_bound and upper_bound is used for limit the span of comparison.
                # Here, there is no limit. So, we set it from the begining to end of genomic position of block.
                # The limit was used in previous version of the code for updating the comparison result after cutting the block (in improver.py).
                lower_bound = var_pos_block[0]
                upper_bound = var_pos_block[-1]

                comparison_result_block = compare_phase_block_pop(allele_block, var_pos_block, pop_inf_dic, lower_bound, upper_bound)
                comparison_result_blocks.append(comparison_result_block)

            report_out_address=out_name_prefix+'_comparison.txt'

            qual_blocks = report_comparison(report_out_address, comparison_result_blocks, id_blocks, chrom)

            report_qc_address=out_name_prefix+'_qc.txt'
            report_qc(report_qc_address, id_blocks, qual_blocks, allele_blocks, stats_vcf,chrom)

            if chrom == chrs_list[0]:
                subprocess.call("head -n 2 "+report_qc_address+" > "+out_prefix+"/QC.csv", shell=True)

            else:
                subprocess.call("sed -n 3,3p "+report_qc_address+" >> "+out_prefix+"/QC.csv", shell=True)



    elif mode_phasme_qc_improver == "improver":


        for chrom in chrs_list:

            print("working on chromosome ", chrom)

            vcf_file_address = out_prefix+"/"+chrom+"/chr"+chrom+".vcf"
            lines_list, var_pos_het_list, line_number_het_list, id_blocks, allele_blocks, var_pos_blocks, stats_vcf, chrom = read_vcf_file(vcf_file_address)
            #print(len(var_pos_het_list))

            file_pairs_address = out_prefix+"/"+chrom+"/"+chrom+"_pairs.txt"

            pop_inf_dic = read_file_pairs_forward(file_pairs_address)

            out_name_prefix= vcf_file_address[:-4] #argv[3]  #

            # comparing input phased vcf (hap_blocks_sample1_dic) with pop (short_blocks_dic_pos)
            comparison_result_blocks=[]

            for block_i in range(len(id_blocks)):

                allele_block = allele_blocks[block_i]
                var_pos_block = var_pos_blocks[block_i]

                # lower_bound and upper_bound is used for limit the span of comparison.
                # Here, there is no limit. So, we set it from the begining to end of genomic position of block.
                # The limit was used in previous version of the code for updating the comparison result after cutting the block (in improver.py).
                lower_bound = var_pos_block[0]
                upper_bound = var_pos_block[-1]

                comparison_result_block = compare_phase_block_pop(allele_block, var_pos_block, pop_inf_dic, lower_bound, upper_bound)
                comparison_result_blocks.append(comparison_result_block)


            report_out_address=out_name_prefix+'_comparison.txt'

            qual_blocks = report_comparison(report_out_address, comparison_result_blocks, id_blocks, chrom)

            report_qc_address=out_name_prefix+'_qc.txt'
            report_qc(report_qc_address, id_blocks, qual_blocks, allele_blocks, stats_vcf,chrom)

            if chrom == chrs_list[0]:
                subprocess.call("head -n 3 "+report_qc_address+" > "+out_prefix+"/QC.csv", shell=True)

            else:
                subprocess.call("sed -n 3,3p "+report_qc_address+" >> "+out_prefix+"/QC.csv", shell=True)



            cut_list_blocks = decide_cut(id_blocks, allele_blocks, var_pos_blocks, comparison_result_blocks)
            #print(cut_list_blocks) # a list of list corresponding to the phase blocks

            lines_list_improved_cut= improve_vcf_cut(lines_list, id_blocks, cut_list_blocks, var_pos_blocks)

            vcf_file_improved_address = vcf_file_address[:-4]+'_improved.vcf'
            write_out_vcf(vcf_file_improved_address, lines_list_improved_cut)

            if chrom == chrs_list[0]:
                subprocess.call("cat "+vcf_file_improved_address+" > "+out_prefix+"/improved.vcf", shell=True)

            else:
                subprocess.call("grep -v \"#\" "+vcf_file_improved_address+" >> "+out_prefix+"/improved.vcf", shell=True)

        print("The QC report is ready at "+out_prefix+"/QC.vcf")
        print("The improved VCF file is ready at "+out_prefix+"/improved.vcf")



    else:

        print("Please specify the sub mode of phaseME correctly: either qc or improver")
