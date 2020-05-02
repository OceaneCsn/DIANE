 
# Dashboard for the Inference and Analysis of Networks from Expression data
![DIANE](www/favicon.ico "DIANE")

DIANE is a shiny application for the analysis of high throughput gene expression data (**RNA-Seq**). The objective is to extract important gene clusters or regulators involved in the response to various perturbations, given the popularity of combinatorial approaches in experimental biology.

This package provides the user with tools to process, explore, and perform advanced statistical analysis on its data using state of the art methods.

+ Count data pre-processing and sample-wise normalisation
+ Differential expression analysis

As several interactive tools already offer those kind of service, we try to go further in the analysis :
+ Expression based clustering in the framework of Poisson Mixture Models
+ Machine learning based Gene regulatory network inference

**DIANE is in a very early stage of development**.

For now, it can be downloaded and installed from this repository with the following R code :

```R
library(remotes)
remotes::install_github("OceaneCsn/DIANE")
```

Packages mentioned as unavailable (because they are not on CRAN) during the installation may have to be manually installed from Bioconductor, usign ```BiocManager::install(...)```.

You can then launch the application with 

```R
library(DIANE)
DIANE::run_app()
```

Author : Océane Cassan

PhD Student at BPMP (Plant Biology and Molecular Physiology) research unit, SUPAGRO Montpellier.

<img src="www/header-logo.png" alt="banner" width="800"/>
