#include <stdio.h>  //puts and printf
#include "htslib/vcf.h"


int main(int argc, char *argv[])
{
    htsFile *inf;
    // counters
    int n    = 0;  // total number of records in file
    int nsnp = 0;  // number of SNP records in file
    // genotype data for each call

    int ngt     = 0;
    int *gt     = NULL;


    // opening vcf file
    inf = bcf_open(argv[1], "r");
    // read header
    bcf_hdr_t *hdr = bcf_hdr_read(inf);
    // #CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	NA00001	NA00002	NA00003	NA00004 A00003
    // should be tab delimited
    int nsmpl = bcf_hdr_nsamples(hdr);
    fprintf(stderr, "File %s contains %i samples\n", argv[1],nsmpl );


    // struc for storing each record
    bcf1_t *rec = bcf_init();
    int32_t *gt_arr = NULL, ngt_arr = 0;
    while (bcf_read(inf, hdr, rec) == 0) {
      n++;
      if (bcf_is_snp(rec)) {
        nsnp++;
      } else {
        continue;
      }
      //ngt = bcf_get_format_int32(hdr, rec, "GT", &gt, &ngt_arr);
      //fprintf(stderr, "gt[0] is %i \n",  bcf_gt_allele(gt) ); //bcf_gt_allele(
      ngt = bcf_get_genotypes(hdr, rec, &gt_arr, &ngt_arr); // ngt is a saclar ploidy*num_samples
      int max_ploidy = ngt/nsmpl;
      //if ( ngt<=0 ) return; // GT not present
      for (int i=0; i<nsmpl; i++)
      {
        int32_t *ptr = gt + i*max_ploidy;
        for (int j=0; j<max_ploidy; j++)
        {
          int is_phased = bcf_gt_is_phased(ptr[j]);


        }
      }

    }
    free(gt_arr);



    fprintf(stderr, "File %s contains %i SNPs\n", argv[1], n );



    bcf_close(inf);
    return 0;
}
