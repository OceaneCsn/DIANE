# Expected expression file

---

You should feed DIANE the raw counts that were obtained after the bioinformatic pipeline of mapping and quantification of your reads.
It gives, for each gene, the transcript aboundance found in each of your experimental condition and repliactes.
The needed matrix-like file shoud contain a column named "Gene", containing all your gene IDs.
The other columns should be your sample names, noted as follow : conditionName_replicate (see example below).
Expression values should be zeros or positive integers.


In order for your input to be compatible with the proposed organisms, it should contain gene IDs as follow :

+ For Arabidopsis thaliana : TAIR IDs (ex: AT1G01020, or AT1G01020.1) 
+ For Human : ensembl IDs (ex: ENSG00000005513) 

Here is the head of our companion expression file, from 24 Arabidopsis thaliana's samples (8 conditions with triplicates) :

```
Gene,C_1,C_2,C_3,S_1,S_2,S_3,M_1,M_2,M_3,H_1,H_2,H_3,SM_1,SM_2,SM_3,SH_1,SH_2,SH_3,MH_1,MH_2,MH_3,SMH_1,SMH_2,SMH_3
AT1G01010.1,127.0,67.9,65.5,94.0,88.1,95.9,65.1,100.3,126.8,95.4,135.0,117.2,96.7,104.4,98.1,94.7,96.1,101.3,82.8,107.4,97.1,100.1,96.7,121.8
AT1G01020.1,207.9,220.8,186.8,192.5,225.1,197.8,234.2,196.9,179.4,312.9,366.0,318.0,169.0,179.6,186.5,340.8,352.6,345.0,331.2,315.8,327.5,267.7,313.7,319.3
AT1G01030.1,32.7,34.4,55.8,33.6,15.9,31.3,21.4,29.8,33.5,47.7,39.0,51.1,25.2,31.2,34.4,61.0,65.5,61.2,45.0,36.2,55.9,46.0,42.7,56.8
```

