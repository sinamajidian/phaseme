#!/bin/sh





for chr in  $(seq 1 22)  ; do 
echo $chr


python /home/ssm/Documents/phaseme/swer.py /home/ssm/Documents/phaseme/ccs/${chr}/${chr}_ccs_true_improved_.90.vcf > results/ccs_${chr}_improved_90.txt


python /home/ssm/Documents/phaseme/swer.py /home/ssm/Documents/phaseme/ccs/${chr}/${chr}_ccs_true.vcf > results/ccs_${chr}.txt

done
