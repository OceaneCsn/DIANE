# Expected expression file

---

You should feed DIANE the raw counts that were obtained after the bioinformatic pipeline of mapping and quantification of your reads.
It gives, for each gene, the transcript aboundance found in each of your experimental condition and repliactes.
The needed matrix-like file shoud contain a column named "Gene", containing all your gene IDs.
The other columns should be your sample names, noted as follow : conditionName_replicate (see example below).
Expression values should be zeros or positive integers.


Here is the head of an example file that would be accepted, from 12 Arabidopsis thaliana's samples :

```
Gene,cNF_3,cNF_2,cNF_1,cnF_2,cnF_1,cnF_3,CNF_1,CnF_2,CnF_1,CnF_3,cNf_1,cnf_2
AT1G01010,1526,1006,1116,1275,967,1018,854,1132,1294,1364,2325,2113
AT1G01020,416,285,289,349,364,370,226,513,502,561,461,407
AT1G01030,31,15,19,29,36,28,12,47,34,47,18,40
```

Please, for now, only comma seperated files are accepted.
