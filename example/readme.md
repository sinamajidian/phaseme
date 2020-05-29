Examples
======


## Precomputed  Mode 

Quality assessment

```
python phaseme.py quality example/my.vcf example/out_pre    
```

Improver

```
python phaseme.py improver example/my.vcf example/out_pre_imp    
```


## Individual  Mode 

Quality assessment

```
python phaseme.py quality example/my.vcf example/out_quality /home/ssm/Documents/phaseme/shapeit_folder /home/ssm/Documents/phaseme/data1/1000g_folder
```

Improver

```
python phaseme.py improver example/my.vcf example/out_imp /home/ssm/Documents/phaseme/shapeit_folder /home/ssm/Documents/phaseme/data1/1000g_folder
```

## Parental mode

Quality assessment

```
python phaseme.py quality example/trio.vcf example/out_trio_quality trio
```

Improver

```
python phaseme.py improver example/trio.vcf example/out_trio_imp trio
```


