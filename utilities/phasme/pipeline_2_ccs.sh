#!/bin/sh

ty='ccs' # ont
#echo ${ty}

for chr in  $(seq 1 22 ) ; do #
echo $chr

cd /home/ssm/Documents/phaseme/${ty}/${chr}

# the folder chr and chr.vcf in it,  are created by population.py 

cp  /home/ssm/Documents/phaseme/data1/grand_truth/${chr}_trio_true.vcf .

for n in ${chr} ${chr}_trio_true ; do
    bgzip -c ${n}.vcf >${n}.vcf.gz
    tabix ${n}.vcf.gz
    bcftools index ${n}.vcf.gz
done

bcftools merge ${chr}.vcf.gz ${chr}_trio_true.vcf.gz > ${chr}_${ty}_true.vcf

python /home/ssm/Documents/phaseme/qc_withtruth.py ${chr}_${ty}_true.vcf ${chr}_pairs_500_0.9.txt ${chr}
python /home/ssm/Documents/phaseme/improver1.py ${chr}_${ty}_true.vcf ${chr}_pairs_500_0.9.txt >${chr}_cutflip.txt

done



