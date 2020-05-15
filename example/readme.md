Examples
======


## Precomputed  Mode 

QC

```
python phaseme.py qc example/my.vcf example/out_pre    
```

Improver

```
python phaseme.py improver example/my.vcf example/out_pre_imp    
```


## Individual  Mode 

QC

```
python phaseme.py qc example/my.vcf example/out_qc /home/ssm/Documents/phaseme/shapeit_folder /home/ssm/Documents/phaseme/data1/1000g_folder
```

Improver

```
python phaseme.py improver example/my.vcf example/out_imp /home/ssm/Documents/phaseme/shapeit_folder /home/ssm/Documents/phaseme/data1/1000g_folder
```

## Parental mode

QC

```
python phaseme.py qc example/trio.vcf example/out_trio_qc trio
```

Improver

```
python phaseme.py improver example/trio.vcf example/out_trio_imp trio
```


