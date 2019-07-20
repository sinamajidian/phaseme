#include <stdio.h>
#include <string.h>
#include "htslib/vcf.h"

#define maxsamples 100 // max population size
#define maxsamplename 10 // max sample name and also pop code

int sample_pop_code(char popfilename[],char sample_names[maxsamplename][maxsamples],char population_codes[maxsamplename][maxsamples]);


int main(int argc, char *argv[]) {
/*
argv[1] phased vcf file names
argv[2] a short string of code of the population of interest
argv[3] name of the sample-membership txt file
    1st_column: sample name, 2nd_column: population code the project

*/

    char popfilename[10]="pop.txt"; //argv[3]
    char sample_names[10][maxsamples], population_codes[maxsamplename][maxsamples];
    int number_samples;
    number_samples= sample_pop_code(popfilename, sample_names, population_codes);




    char popname[10]= "GBR"; // argv[2]
    char buf[256];
    for(int i = 1; i<number_samples; i++) // number_samples in popfilename
    {

      if (strcmp(population_codes[i],popname)==0){
      //fprintf(stderr,"  Element %dth is %s and code is %s \n",i, sample_names[i], population_codes[i]);
        fprintf(stderr,"  Element %dth is %s and code is %s \n",i, sample_names[i], population_codes[i]);
        snprintf(buf, sizeof buf, "%s%s%s", sample_names[i],",", buf);
        fprintf(stderr,"  list of samples is  %s \n \n ", buf);
       }
    }
    // fprintf(stderr," ");

    char sample_list[maxsamples*3+1]="A1,A2";


    htsFile *inf;
    char vcfname[10]="my.vcf" ;  // argv[1]
    inf = bcf_open(vcfname, "r");    // opening vcf file
    bcf_hdr_t *hdr = bcf_hdr_read(inf); // reading the header
    // the list in the next line should be without space no  "A1, A2"
    int subsmaple_results =bcf_hdr_set_samples(hdr,sample_list, 0); // third elemtn  is a file (1) or a comma-separated list (0)
    //Returns 0 on success, -1 on error or a positive integer if the list contains samples not present in the VCF header.

    if (subsmaple_results==0) {fprintf(stderr, "subsmaple is perfect! \n");}
    if (subsmaple_results>0) {fprintf(stderr, " the list contains samples not present in the VCF header\n");}

    // #CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	NA00001	NA00002	NA00003	NA00004 A00003
    // should be tab delimited (no space)
    int nsmpl = bcf_hdr_nsamples(hdr);
    fprintf(stderr, "For the mentioned population, %s contains %i sample(s) \n", vcfname,nsmpl );


    return 0;
}





    //
    // // struc for storing each record
    // bcf1_t *rec = bcf_init();
    // int32_t *gt_arr = NULL, ngt_arr = 0;
    // while (bcf_read(inf, hdr, rec) == 0) {
    //   n++;
    //   if (bcf_is_snp(rec)) {
    //     nsnp++;
    //   } else {
    //     continue;
    //   }
    //   //ngt = bcf_get_format_int32(hdr, rec, "GT", &gt, &ngt_arr);
    //   //fprintf(stderr, "gt[0] is %i \n",  bcf_gt_allele(gt) ); //bcf_gt_allele(
    //   ngt = bcf_get_genotypes(hdr, rec, &gt_arr, &ngt_arr); // ngt is a saclar ploidy*num_samples
    //   int max_ploidy = ngt/nsmpl;
    //
    //   }
    //
    // }
    // free(gt_arr);



int sample_pop_code(char popfilename[], char sample_names[maxsamplename][maxsamples], char population_codes[maxsamplename][maxsamples])
{
/*
This function is for extracting sample names and the corresponding population names
  from a text file containing two cloumns.

Output: The sample and population names are saved in the second and third argument array.
        The number of samples is also reported as the ouput.

*/
    FILE* fp = fopen(popfilename,"r");
    int counter_sample = 1;
    char sample_name[maxsamplename], population_code[maxsamplename];

    while (fscanf(fp,"%s %s",sample_name, population_code)==2){
        strcpy(sample_names[counter_sample],sample_name);
        strcpy(population_codes[counter_sample],population_code);
        counter_sample++;
    }
  return counter_sample-1;
}
