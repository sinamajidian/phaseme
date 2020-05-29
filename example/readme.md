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

## Parental mode

Quality assessment

```
python phaseme.py quality example/trio.vcf example/out_trio_quality trio
```

Improver

```
python phaseme.py improver example/trio.vcf example/out_trio_imp trio
```




## Individual Mode 
It needs Shapeit.

Quality assessment

```
python phaseme.py quality example/my.vcf example/out_quality /path/to/shapeit/ /path/to/1000G
```

Improver

```
python phaseme.py improver example/my.vcf example/out_imp /path/to/shapeit/ /path/to/1000G
```

