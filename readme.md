Phaseme
======





## Prerequisite  


1- shapeit

2- 1000 Genomes reference panel haplotypes from [here](https://mathgen.stats.ox.ac.uk/impute/1000GP_Phase3.html)








## input file

- phased VCF





## Step 1: Optain population information over Shapeit


```

shapeit -check --input-vcf  your_out.vcf -R data/1000GP_Phase3/1000GP_Phase3_chr${chr}.hap.gz  data/1000GP_Phase3/1000GP_Phase3_chr${chr}.legend.gz  data/1000GP_Phase3/1000GP_Phase3.sample   --output-log out

shapeit --input-vcf out.vcf  -R ../../data/1000GP_Phase3/1000GP_Phase3_chr${chr}.hap.gz  data/1000GP_Phase3/1000GP_Phase3_chr${chr}.legend.gz  data/1000GP_Phase3/1000GP_Phase3.sample  -M data/1000GP_Phase3/genetic_map_chr${chr}_combined_b37.txt --output-log out1 --output-graph out.graph --exclude-snp  out.snp.strand.exclude

```



## Step 2: Run PhaseMe to obtain stats and improve the quality of phase blocks

```

python2 utilities/samplehaps.py out 500 >log_samplehaps
python utilities/encoderead3.py out.hapsamples
python utilities/improve.py out.vcf pairs.txt
```

## Future steps

The package contains two parts.


1- Reporting quality of phased VCF (including estimated haplotypes).

2- Improving the phased VCF 


