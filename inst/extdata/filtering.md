# Removing low count genes


Removing genes with very low aboundance is a common practice in RNA-Seq analysis pipelines for several reasons :


+ The have little biological signifiance, and could be caused either by noise or mapping errors.


+ The statitical modelling we are planning to perform next is not well suited for low counts, as they make the mean-variance relationship harder to estimate.

There is no absolute and commonly accpeted threshold value, but it is recommended to allow only genes with more than 10 counts per sample in average. DIANE thus proposes a threshold at 10*sampleNumber, but feel free to experiment with this value depending on your dataset.
