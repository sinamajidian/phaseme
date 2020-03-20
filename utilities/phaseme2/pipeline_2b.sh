#!/bin/sh

ty=$1
chr=$2
pair_address=$3  # "/home/ssm/Documents/phaseme_v2/vcf/${chr}_fake_pairs_500_0.9.txt"

cd ${ty}/${chr}

python /home/ssm/Documents/phaseme/qc_withtruth.py ${chr}_${ty}_true.vcf ${pair_address}  ${chr}
python /home/ssm/Documents/phaseme/improver2.py ${chr}_${ty}_true.vcf ${pair_address}  >${chr}_cutflip.txt




