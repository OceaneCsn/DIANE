# Dashboard for the Inference and Analysis of Networks from Expression data <img src="man/figures/hex-DIANE.png" align="right" alt="" width="120" />

## Application presentation

DIANE is a R-Shiny application for the analysis of high throughput gene expression data (**RNA-Seq**). Its function is to extract important regulatory pathways involved in the response to environmental changes, or any perturbation inducing genomic modifications.

Given the popularity of combinatorial approaches in experimental biology, we designed this tool to process, explore, and perform advanced statistical analysis on **multifactorial expression data** using state of the art methods. It includes :

+ Raw count data pre-processing and normalisation

+ Differential expression analysis and results visualization (Volcano plots, heatmaps, Venn diagrams...)

+ Gene ontology enrichment analysis

+ Expression based clustering in the framework of Mixture Models, and individual characterisation of those clusters (generalized linar models and GO enrichment analysis).

+ Machine learning based Gene regulatory network inference and interactive network analysis, communi

+ Session reporting and results download

+ Demonstration dataset and other ready to be explored datasets on several organisms

All of the features in DIANE are accessible via a signle page Shiny application that can be locally launched, or used online at https://diane.bpmp.inrae.fr.

<img src="man/figures/net.PNG" align="center" alt="" width="900" />

For users more familiar with R programming, all server-side functions in DIANE are exported so they can be called from R scripts. 
Those functions documentation can be found in the [Reference](https://oceanecsn.github.io/DIANE/reference/index.html), and are illustrated in the [corresponding vignette](https://oceanecsn.github.io/DIANE/articles/DIANE_Programming_Interface.html)


## Use DIANE locally

DIANE relies on R >= 4.0.1, available for all OS at https://cloud.r-project.org/.

Download and install DIANE in your R console as follows (you need the remotes package installed ```install.packages("remotes")```) :

```R
remotes::install_github("OceaneCsn/DIANE")
```

You can then launch the application :

```R
library(DIANE)
DIANE::run_app()
```
In case your expression input file exceeds 5MB, you may need to run the command ```options(shiny.maxRequestSize=30*1024^2)``` before calling ```DIANE::run_app()``` to upload up to 30BM.

Once the application is launched, if the resolution poorly fits your screen, you can adjust it with the keyboard shortcuts ```ctrl +``` or  ```ctrl -``` (use ```cmd``` on Mac).
 

Authors : Océane Cassan, Antoine Martin, Sophie Lèbre

PhD Student at BPMP (Plant Biology and Molecular Physiology) research unit, SUPAGRO Montpellier.
