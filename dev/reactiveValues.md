
# Reactive values in DIANE shared accross modules
---

Variable reactive globale **r** :
```
r
| raw_counts (dataframe : (genes * samples))
| normalized_counts (dataframe: (genes * samples))
| normalized_counts_pre_filter (dataframe : (genes * samples))
| tcc (TCC class object)
| conditions (vector of the samples condition names)
| design (dataframe : (conditionNames * factors))
| regulators
| gene_info dataframe 
| splicing_aware (boolean)
| DEGs
|   | ref trt (vector of genes)
| top_tags
|   | ref trt (dataframe : (genes * (logFC, logCPM, FDR))
| fit (dispersion and glmFit for dea)
| clusterings
|   | ref trt 
|   |    | model (coseqResult)
|   |    | membership (named vector)
|   |    | conditions (vector)
| networks
|   | ref trt 
|   |    | nodes
|   |    | edges
|   |    | membership clustering louvain(named vector)
|   |    | conditions
|   |    | graph
``` 

Juste pour savoir ce qui est accessible à tout moment dans l'appli.

Peut être aussi pour être stockée dans une session et reloadé direct.


Organism specific :

+ GO analysis, change id to entrez

+ gene info

+ TF list for network inference