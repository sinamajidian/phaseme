#!/bin/sh



tr=20
ty='ont'
#for tr in 5 10 20 50 100; do
echo ${tr}
#mkdir results_${tr}
for chr in $(seq 1 3)  ; do
echo $chr

#for ty in 'ont' 'ccs' 'clr'; do 

python /home/ssm/Documents/phaseme/swer.py /home/ssm/Documents/phaseme/${ty}/${chr}/${chr}_${ty}_true_improved_.90.vcf ${tr} > results_${tr}/${ty}_${chr}_improved_90.txt
python /home/ssm/Documents/phaseme/swer.py /home/ssm/Documents/phaseme/${ty}/${chr}/${chr}_${ty}_true.vcf ${tr} > results_${tr}/${ty}_${chr}.txt

done
#done
#done





