#!/bin/sh


# bcftools merge arises an error of different ADALL lengths. It happens when there are duplicated position. 

#bcftools annotate -x INFO,^FORMAT/GT HG003_GRCh37_GIAB_highconf_CG-IllFB-IllGATKHC-Ion-10X_CHROM1-22_v.3.3.2_all.vcf > hg3.vcf
#bcftools annotate -x INFO,^FORMAT/GT HG004_GRCh37_GIAB_highconf_CG-IllFB-IllGATKHC-Ion-10X_CHROM1-22_v.3.3.2_all.vcf > hg4.vcf
#bcftools annotate -x INFO,^FORMAT/GT NA24385_LongRanger_snpindel.vcf  > hg2_10x.vcf


for chr in $(seq 1 21); do
echo $chr
grep "#" hg3.vcf | sed 's/INTEGRATION/HG003/' > ${chr}_hg3.vcf 
grep "^${chr}\b" hg3.vcf >> ${chr}_hg3.vcf

grep "#" hg4.vcf | sed 's/INTEGRATION/HG004/' > ${chr}_hg4.vcf 
grep "^${chr}\b" hg4.vcf >> ${chr}_hg4.vcf



grep "#" hg2_10x.vcf | sed 's/FORMAT    61199/FORMAT    HG002_10x/' | sed 's/Type=Integer,Description="ID of Phase Set for Variant">/Type=String,Description="Phase set in which this variant falls">/' > chr${chr}_10x.vcf

grep "^chr${chr}\b" hg2_10x.vcf >> chr${chr}_10x.vcf
sed "s/chr${chr}/${chr}/g" chr${chr}_10x.vcf   > ${chr}_10x.vcf


for n in ${chr}_10x ${chr}_hg3 ${chr}_hg4; do
    bgzip -c ${n}.vcf >${n}.vcf.gz
    tabix ${n}.vcf.gz
    bcftools index ${n}.vcf.gz
done


bcftools merge ${chr}_10x.vcf.gz ${chr}_hg3.vcf.gz ${chr}_hg4.vcf.gz > ${chr}_trio.vcf

python trio_grand_truth.py ${chr}_trio.vcf



done
