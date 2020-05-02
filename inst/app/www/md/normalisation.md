
# Normalisation


The **TCC** package is used here, to make samples comparable by correcting for their differences in **sequencing depths**. This step is mandatory before further statitical analysis.

You can choose to normalize using the methods implemented in edge, 'tmm', or the one of 'deseq2'.

Those normalisation methods rely on the hypothesis that few genes are differentially expressed between your samples. If you suspect this is not the case in your data, TCC offers the possibility to proceed to a first detection of differentially expressed genes, to remove them, and then provide a final less biased normalisation.

In that case, enable "prior removal of differentially expressed genes".
TCC will perform the following setp, depending on the normalisation method you chose :
+ tmm/deseq2
+ DEG identification and removal using edgeR test method
+ tmm/deseq2
