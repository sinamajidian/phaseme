#!/bin/sh



ty=$1
chr=$2
tr=20
add=$3
#for tr in 5 10 20 50 100; do
#echo ${tr}
#mkdir results_${tr}
#for chr in $(seq 1 3)  ; do
echo $chr

#for ty in 'ont' 'ccs' 'clr'; do 

python /home/ssm/Documents/phaseme/swer5.py ${add}_true_improved_.90.vcf ${tr} > result/${ty}_${chr}_improved_90.txt
python /home/ssm/Documents/phaseme/swer5.py  ${add}_true.vcf ${tr} > result/${ty}_${chr}.txt

#done
#done
#done





