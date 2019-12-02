#!/bin/sh





for chr in  8; do  # $(seq 19 20)  ; do 
echo $chr


python /home/ssm/Documents/phaseme/swer.py /home/ssm/Documents/phaseme/ont/${chr}/${chr}_ont_true_improved_.90.vcf > results/${chr}_ont_improved_90.txt


python /home/ssm/Documents/phaseme/swer.py /home/ssm/Documents/phaseme/ont/${chr}/${chr}_ont_true.vcf > results/${chr}_ont.txt

done
