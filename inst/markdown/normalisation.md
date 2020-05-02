
# Normalisation methods

---

The [TCC R package](https://rdrr.io/bioc/TCC) is used for the normalisation step here, to make samples comparable by correcting for their differences in **sequencing depths**. This step is mandatory before further statitical analysis.

You can choose to normalize using the methods implemented in edgeR, referenced as 'tmm', or the one used in DESeq, referenced as 'deseq2'.

Those normalisation methods rely on the hypothesis that a very small proportion of genes are differentially expressed between your samples. If you suspect a lot of genes could be differentially expressed in your data, TCC offers the possibility to proceed to a first detection of potential differentially expressed genes, to remove them, and then provide a final less biased normalisation.

In that case, enable "prior removal of differentially expressed genes".
TCC will perform the following setp, depending on the normalisation method you chose :
+ tmm/deseq2 temporary normalisation
+ potential DEG identification and removal using edgeR test method
+ tmm/deseq2 definitive normalisation   