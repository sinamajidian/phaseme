#!/bin/sh





for chr in   $(seq 1 22)  ; do 
echo $chr


python /home/ssm/Documents/phaseme/swer.py /home/ssm/Documents/phaseme/clr/${chr}/${chr}_clr_true_improved_.90.vcf > results/clr_${chr}_improved_90.txt


python /home/ssm/Documents/phaseme/swer.py /home/ssm/Documents/phaseme/clr/${chr}/${chr}_clr_true.vcf > results/clr_${chr}.txt

done
