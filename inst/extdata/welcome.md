 

# Dashboard for the Inference and Analysis of Networks from Expression data <img src="www/favicon.ico" align="right" alt="" width="180" />

---

DIANE is a shiny application for the analysis of high throughput gene expression data (**RNA-Seq**). Its function is to extract important regulatory pathways involved in the response to environmental changes, or any perturbation inducing genomic modifications.

We designed this tool to process, explore, and perform advanced statistical analysis on **multifactorial expression data** using state of the art methods. It includes :

+ Raw count data pre-processing and sample-wise normalisation
+ Customizable differential expression analysis

+ Gene ontology enrichment analysis for model organisms

+ Expression based clustering in the framework of Poisson Mixture Models, and characterisation of those clusters with generalized linar models and GO enrichment analysis

+ Machine learning based Gene regulatory network inference

All of the features in DIANE are accessible via a signle page shiny application that can be locally launched.

<img src="www/DIANE.png" alt="banner" width="1200" align="center"/>


For more advanced users, all server-side functions in DIANE are exported so they can be called from R scripts. 

Fore more information, please find full documentation and examples in the github page  https://oceanecsn.github.io/DIANE.

Once the application is launched, if the resolution poorly fits your screen, you can adjust it with the keyboard shortcuts ```ctrl +``` or  ```ctrl -``` (use ```cmd``` on Mac).



**DIANE is in beta version, please report any bug or suggestion via github or at oceane.cassan@supagro.fr**.


To use DIANE locally, download and install DIANE in your R console as follows (you need the remotes package installed) :

```R
remotes::install_github("OceaneCsn/DIANE")
```

DIANE relies on R 4.0.0, available for all OS at https://cloud.r-project.org/.

You can then launch the application :

```R
library(DIANE)
DIANE::run_app()
```

Authors : Océane Cassan, Antoine Martin, Sophie Lèbre

Dev : Océane Cassan, PhD Student at BPMP (Plant Biology and Molecular Physiology) research unit, SUPAGRO Montpellier.
