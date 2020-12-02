# Dashboard for the Inference and Analysis of Networks from Expression data <img src="www/favicon.ico" align="right" width="180"/>

DIANE is a shiny application for the analysis of high throughput gene expression data (**RNA-Seq**). Its function is to extract important regulatory pathways involved in the response to environmental changes, or any perturbation inducing genomic modifications.

We designed this tool to process, explore, and perform advanced statistical analysis on **multifactorial expression data** using state of the art methods. It includes :

-   Raw count data pre-processing and sample-wise normalisation

-   Customizable differential expression analysis

-   Gene ontology enrichment analysis for model organisms

-   Expression based clustering in the framework of Poisson Mixture Models, and characterization of those clusters with generalized linar models and GO enrichment analysis

-   Machine learning based Gene regulatory network inference

All of the features in DIANE are accessible via a single page shiny application that can be locally launched, or accessed online at <https://diane.bpmp.inrae.fr>.

The steps should be performed in the order of the different tabs. For instance, before running clustering or network inference, differential expression should be performed first. The figure above summarizes DIANE's main workflow.

<img src="www/DIANE.png" alt="banner" width="1200" align="center"/>

For more advanced users, all server-side functions in DIANE are exported so they can be called from R scripts.

Fore more information, please find full documentation and examples in the github page <https://oceanecsn.github.io/DIANE>.

Once the application is launched, if the resolution poorly fits your screen, you can adjust it with the keyboard shortcuts `ctrl +` or `ctrl -` (use `cmd` on Mac).

**DIANE is in beta version, please report any bug or suggestion via github or at [oceane.cassan\@supagro.fr](mailto:oceane.cassan@supagro.fr){.email}**.

## Use DIANE locally

To use DIANE locally, download and install DIANE in your R console as follows (you need the remotes package installed) :

``` {.r}
remotes::install_github("OceaneCsn/DIANE")
```

DIANE relies on R 4.0.0, available for all OS at <https://cloud.r-project.org/>.

You can then launch the application :

``` {.r}
library(DIANE)
DIANE::run_app()
```

## Deploy DIANE on your server

We provide a [solution based on Docker and Shiny server](https://hub.docker.com/r/rocker/shiny) to deploy DIANE on any linux server. To do so, see the following command line instructions.

Get DIANE source code via Git :

    git clone https://github.com/OceaneCsn/DIANE.git

Install Docker engine, as described in the [Docker docs](https://docs.docker.com/engine/install/).

Go to DIANE's folder.

First, you can change the default settings for the dockeried shiny-server by editing the file shiny-customized.config (like changing the port, the user to run with, and more)

Now let's build the image, that we'll name diane, from the Dockerfile (superuser rights required).

    cd DIANE
    docker build -t diane .

This might take a while. Then, you can start the container diane :

    docker run -d --cpus 16 --user shiny --rm -p 8086:8086 -v /path/to/app/on/host/:/srv/shiny-server/ -v /path/to/logs/on/host/:/var/log/shiny-server/ diane

In this example, a session of DIANE will be allowed to use 16 CPU cores.

`/path/to/app/on/host/` is the path to DIANE on the host, the location where you cloned it, containing the app.R file. `/path/to/logs/on/host/` is the folder you want to store your app logs.

`-p 8086:8086` is the port to use, change the first 8086 to use another one on the host.

`--user shiny` allows to run as non root, with the shiny user that should have been created before, and granted rights to the folder `/path/to/app/on/host/logs` and `/path/to/logs/on/host/`.

`-d --rm` are options for the detached mode.

You can check that the container is running with `docker ps`.

---

Authors : Océane Cassan, Antoine Martin, Sophie Lèbre

Dev : Océane Cassan, PhD Student at BPMP (Plant Biology and Molecular Physiology) research unit, SUPAGRO Montpellier.
