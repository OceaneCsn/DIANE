# Dashboard for the Inference and Analysis of Networks from Expression data <img src="www/favicon.ico" align="right" width="180"/>

DIANE is a shiny application for the analysis of high throughput gene expression data (**RNA-Seq**). Its function is to extract important regulatory pathways involved in the response to environmental changes, or any perturbation inducing genomic modifications.

We designed this tool to process, explore, and perform advanced statistical analysis on **multifactorial expression data** using state of the art methods. It includes :

-   Raw count data pre-processing and sample-wise normalization

-   Customizable differential expression analysis

-   Gene ontology enrichment analysis for model organisms

-   Expression based clustering in the framework of Poisson Mixture Models, and characterization of those clusters with generalized linear models and GO enrichment analysis

-   Machine learning based Gene regulatory network inference

All of the features in DIANE are accessible via a single page shiny application that can be locally launched, or accessed online at [https://diane.ipsim.inrae.fr](https://diane.bpmp.inrae.fr){.uri}.

The steps should be performed in the order of the different tabs. For instance, before running clustering or network inference, differential expression should be performed first. The figure above summarizes DIANE's main workflow.

<img src="www/DIANE.png" alt="banner" width="1200" align="center"/>

For more advanced users, all server-side functions in DIANE are exported so they can be called from R scripts.

Fore more information, please find full documentation and examples in the github page <https://oceanecsn.github.io/DIANE>.

Once the application is launched, if the resolution poorly fits your screen, you can adjust it with the keyboard shortcuts `ctrl +` or `ctrl -` (use `cmd` on Mac).

**Please report any bug or suggestion via github or at oceane.cassan\@cnrs.fr**.

**To cite DIANE in publications use:**

Cassan, O., Lèbre, S. & Martin, A. Inferring and analyzing gene regulatory networks from multi-factorial expression data: a complete and interactive suite. BMC Genomics 22, 387 (2021). <https://doi.org/10.1186/s12864-021-07659-2>

**A BibTeX entry for LaTeX users is**

    @Article{cassan2021Inferring,
        title = {Inferring and analyzing gene regulatory networks from multi-factorial expression data: a complete and interactive suite},
        author = {Océane Cassan and Sophie Lèbre and Antoine Martin},
        journal = {BMC Genomics},
        year = {2021},
        volume = {22},
        number = {387},
        url = {https://bmcgenomics.biomedcentral.com/articles/10.1186/s12864-021-07659-2}}

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

------------------------------------------------------------------------

## License


Copyright (C) 2020 Oceane Cassan

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.


---

## Deploy DIANE on your server

We provide a [solution based on Docker and Shiny server](https://hub.docker.com/r/rocker/shiny) to deploy DIANE on any linux server, just as it is at <https://diane.bpmp.inrae.fr>. To do so, see the following command line instructions.

Get DIANE source code via Git :

    git clone https://github.com/OceaneCsn/DIANE.git

Install Docker engine, as described in the [Docker docs](https://docs.docker.com/engine/install/).

Go to DIANE's folder.

First, you can change the default settings for the dockerized shiny-server by editing the file shiny-customized.config (like changing the port, the user to run with, and more)

Now let's build the image, that we'll name diane, from the Dockerfile (superuser rights required).

    cd DIANE
    docker build -t diane .

This might take a while. You can check that the container image was built with `docker images`. Then, you can start the container diane, by setting appropriately the following options in the above command:

`/path/to/app/on/host/` is the path to DIANE on the host, that is to say the location where you cloned it, containing the app.R file. `/path/to/logs/on/host/` is the folder you want to store your app logs.

`-p 8086:8086` is the port to use, change the first 8086 to use another one on the host.

`--user shiny` allows to run as non root, with the shiny user that must have been created before, and granted rights to the folder `/path/to/app/on/host/logs` and `/path/to/logs/on/host/`.

`-d --rm` are options for the detached mode.

Plus, in the following example, a session of DIANE will be allowed to use 16 CPU cores :

    docker run -d --cpus 16 --user shiny --rm -p 8086:8086 -v /path/to/app/on/host/:/srv/shiny-server/ -v /path/to/logs/on/host/:/var/log/shiny-server/ diane

You can check that the container is running with `docker ps`.

------------------------------------------------------------------------

Authors : Océane Cassan, Antoine Martin, Sophie Lèbre

Dev : Océane Cassan, PhD Student at IPSIM (Institute for Plant Sciences in Montpellier) research unit, SUPAGRO Montpellier, with contributions from Alexandre Soriano.
