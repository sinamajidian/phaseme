Phaseme
======





## Prerequisite  


1- shapeit

2- 1000 Genomes reference panel haplotypes from [here](https://mathgen.stats.ox.ac.uk/impute/1000GP_Phase3.html)








## input file

- phased VCF





## Step 1.


```



inputvcf='hg002_ont_53x_minimap_2019.vcf'


mkdir ${chr}; cd ${chr}

grep "#" ../${inputvcf} > out.vcf
grep -P "^${chr}\b" ../${inputvcf}  | grep -v "/" | grep -v "\.:\.:\.">> out.vcf

shapeit -check --input-vcf  out.vcf -R ../../data/1000GP_Phase3/1000GP_Phase3_chr${chr}.hap.gz  ../../data/1000GP_Phase3/1000GP_Phase3_chr${chr}.legend.gz  ../../data/1000GP_Phase3/1000GP_Phase3.sample   --output-log out

shapeit --input-vcf out.vcf  -R ../../data/1000GP_Phase3/1000GP_Phase3_chr${chr}.hap.gz  ../../data/1000GP_Phase3/1000GP_Phase3_chr${chr}.legend.gz  ../../data/1000GP_Phase3/1000GP_Phase3.sample  -M ../../data/1000GP_Phase3/genetic_map_chr${chr}_combined_b37.txt --output-log out1 --output-graph out.graph --exclude-snp  out.snp.strand.exclude

python2 ../utilities/samplehaps.py out 500 >log_samplehaps
python ../utilities/encoderead3.py out.hapsamples
 
```




## Step 2.


```
python ../utilities/improve.py out.vcf pairs.txt
```








## Citation:

[IntegratedPhasing](https://github.com/vibansal/IntegratedPhasing)

[shapeit](https://mathgen.stats.ox.ac.uk/genetics_software/shapeit/shapeit.html)





## future



The package contains two parts.


1- Reporting quality of phased VCF (including estimated haplotypes).

2- Improving the phased VCF 




