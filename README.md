 
# Dashboard for the Inference and Analysis of Networks from Expression data
![alt text](inst/app/www/favicon.ico "DIANE")

DIANE is a shiny application for the analysis of high throughput gene expression data (RNA-Seq). The objective is to extract important gene clusters or regulators involved in the response to various perturbation factors, given the popularity of combinatorial approaches in experimental biology.

The package provides the user with tools to process, explore, and perform advanced statistical analysis on its data.

It starts with count data processing, normalisation, and then proceeds to differential expression analysis using state of the art methods.

As several interactive tools already offer those kind of data exploration, we try to go further in the analysis by proposing expression based clustering, as well as gene regulatory network inference.

**DIANE is in a very early stage of development**.

For now, it can be downloaded and installed from this repository with the following R code :

```R
library(remotes)

remotes::install_github("OceaneCsn/DIANE")
```


Packages mentioned as unavailable (because they are not on CRAN) during the installation may have to be manually installed from Bioconductor.

You can then launch the application with 

```R
library(DIANE)

DIANE::run_app()
```


