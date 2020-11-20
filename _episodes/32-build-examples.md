---
title: "Breakout room 3: design container images"
teaching: 0
exercises: 20
questions:
objectives:
- Write Dockerfiles for real-world applications
keypoints:
- Picking the appropriate base image is key and can save you lots of work
- The most important instructions in a Dockerfile are 'FROM' to select the base image and 'RUN' to execute commands
- Build a container image with Docker using `docker build`
- Convert a Docker image into Singularity format by using `singularity pull docker-daemon:`
---


### Goal

In this session, you're going to use Docker to build two container images, in particular one RStudio image and one Conda image.  For each of them, you'll pick an appropriate base image, write the Dockerfile and action the build.

The first step in designing a Dockerfile is to choose a base image, that is the starting point for our container.

The best place to find useful base image is the [Docker Hub](https://hub.docker.com) online registry.  Here is a non-comprehensive list of potentially useful base images (the version tags are not necessarily the most recent ones, but they're relevant to this exercise):

* OS images, such as `ubuntu:18.04`, `debian:buster` and `centos:7`;
* R images, in particular the `rocker/` repository by the [Rocker project](https://www.rocker-project.org):
    - Base R: `rocker/r-ver:3.6.1`;
    - RStudio: `rocker/rstudio:3.6.1`;
    - Tidyverse+RStudio: `rocker/tidyverse:3.6.1`;
* Python images, such as `python:3.8` and `python:3.8-slim` (a lightweight version);
* Conda images by [Anaconda](https://www.anaconda.com), such as `continuumio/miniconda3:4.8.2`;
* Jupyter images, in particular the `jupyter/` repository by [Jupyter Docker Stacks](https://jupyter-docker-stacks.readthedocs.io) (unfortunately making extensive use of the `latest` tag), for instance:
    - Base Jupyter: `jupyter/base-notebook:latest`;
    - Jupyter with scientific Python packages: `jupyter/scipy-notebook:latest`;
    - *Data science* Jupyter, including Python, scientific Python packages, R, Tidyverse and more: `jupyter/datascience-notebook:latest`.


### Exercise 1: Write an RStudio Dockerfile

The first exercise in this session is to write a little Dockerfile for the R package `ggtree`.  This package provides functionalities to represent phylogenetic trees using R, by building on top of the [Tidyverse](https://www.tidyverse.org) collection of data science packages, and in particular the plotting library `ggplot2`.  

First, `cd` into the appropriate directory:

```bash
cd /data/abacbs-containers/exercises/build/r-ggtree
```

Use a text editor to create a blank `Dockerfile`; both `nano` and `vi` are available, pick your favourite.


> ## Choose the base image
> 
> Considering the characteristic that we have just stated for `ggtree`, an R package based on Tidyverse, which is the closest base image you would choose from the list of useful images above?
> 
> > ## Solution
> > 
> > We need R and Tidyverse, so let's go with `rocker/tidyverse:3.6.1`.
> {: .solution}
> 
> Now we need to declare this image in the Dockerfile, using an appropriate Docker instruction.
> 
> > ## Solution
> > 
> > ```source
> > FROM rocker/tidyverse:3.6.1
> > ```
> {: .solution}
{: .challenge}


> ## Command to install an R package
> 
> The package `ggtree` is part of the BioConductor project.  So, from inside an R console we could install it by using the command `BiocManager::install("ggtree")`.  However, here we need Docker to execute this command from a bash shell.
> 
> If you are an R user, do you know how you can execute the R command above from the shell?  No worries if you don't, just have a look at the solution.
> 
> > ## Solution
> > 
> > ```source
> > R -e 'BiocManager::install("ggtree")'
> > ```
> {: .solution}
{: .challenge}
> 
> We're almost there with our first Dockerfile... now we just need to embed the shell command above in the Dockerfile, by using the appropriate Docker instruction.
> 
> > ## Solution
> > 
> > ```source
> > FROM rocker/tidyverse:3.6.1
> > 
> > RUN R -e 'BiocManager::install("ggtree")'
> > ```
> > 
> > Very often (albeit not always) preparing Dockerfiles for R images looks as simple as this.  Other times, you will also need to install other packages such as pre-requisites.  And of course things can even get more complicated than this.
> {: .solution}
{: .challenge}


> ## Document your Dockerfile with labels
> 
> It's a good practice to add information in your Dockerfile in the form of "labels", for instance an email contact for the person who developed it.  
> In this case, just add your name, using the appropriate Docker instruction.
> 
> > ## Solution
> > 
> > ```source
> > FROM rocker/tidyverse:3.6.1
> > 
> > LABEL maintainer="Myself" 
> >
> > RUN R -e 'BiocManager::install("ggtree")'
> > ```
> {: .solution}
{: .challenge}


> ## Build the image
> 
> Now it's finally the time for building!  
> Let' use the appropriate `docker` syntax in the shell to build an image called `ggtree:2.0.4`; remember we're running from the directory where the Dockerfile is (`.`). 
> 
> > ## Solution
> > 
> > ```bash
> > sudo docker build -t ggtree:2.0.4 .
> > ```
> > 
> > It will only take a couple of minutes to build, as most required R packages are already provided by the base Tidyverse image.
> {: .solution}
> 
> Note how here we provided the information about the package version; normally we would be able to find it out ourselves after the first build, by inspecting the R installation in the container.
{: .challenge}


### Converting the Docker image into Singularity format

We're not doing it now to save time, but remember that the last step required to use the built image with Singularity is to turn it into a SIF file.  You'll need a Singularity installation in conjunction with Docker:

```bash
singularity pull docker-daemon:ggtree:2.0.4
```


### Exercise 2: Write a Conda Dockerfile

In this second exercise you're going to build an image for the popular bioinformatics tool `samtools`; you'll use the *Conda* package manager to install it.

First, `cd` into the appropriate directory, and then create a blank `Dockerfile` with a text editor:

```bash
cd /data/abacbs-containers/exercises/build/conda-samtools
```


> ## First, pick the base image
> 
> Have a look back at the list of suggested base images above; which one would you pick for a Conda installation?
> 
> > ## Solution
> > 
> > The best match is `continuumio/miniconda3:4.8.2`.
> {: .solution}
> 
> Now embed the base image declaration in the Dockerfile.
> 
> > ## Solution
> > 
> > ```source
> > FROM continuumio/miniconda3:4.8.2
> > ```
> {: .solution}
{: .challenge}


> ## Second, write the command to install a Conda package
> 
> Samtools version 1.9 can be installed with Conda through the channel *bioconda*, so the shell syntax would be (`-y` is to confirm prompts):
> 
> ```source
> conda install -y -c bioconda samtools=1.9
> ```
> 
> With this information, complete the Dockerfile to install samtools using Conda.
> 
> > ## Solution
> > 
> > ```source
> > FROM continuumio/miniconda3:4.8.2
> > 
> > RUN conda install -y -c bioconda samtools=1.9
> > ```
> > 
> > Similar to the case of R, Dockerfiles often turn out to be quite compact when using a `conda` base image.  Things are not always this easy, as for instance package version conflicts are common with Conda;  additional command lines might be required to work around them, or you might even need to install packages entirely manually.
> {: .solution}
> 
> In the interest of time, you're not going to build this Conda-based image.
{: .challenge}


> ## Bonus: example Dockerfiles
> 
> You may have a look at these, to get a taste of what more articulated Dockerfiles look like.
> 
> > ## A large R image
> >
> > ```bash
> > FROM rocker/tidyverse:latest
> > 
> > RUN apt-get update -qq && apt-get -y --no-install-recommends install \
> >       autoconf \
> >       automake \
> >       g++ \
> >       gcc \
> >       gfortran \
> >       make \
> >       && apt-get clean all \
> >       && rm -rf /var/lib/apt/lists/*
> > 
> > RUN mkdir -p $HOME/.R
> > COPY Makevars /root/.R/Makevars
> > 
> > RUN Rscript -e "library('devtools')" \
> >       -e "install_github('Rdatatable/data.table', build_vignettes=FALSE)" \
> >       -e "install.packages('reshape2')" \
> >       -e "install.packages('fields')" \
> >       -e "install.packages('ggbeeswarm')" \
> >       -e "install.packages('gridExtra')" \
> >       -e "install.packages('dynamicTreeCut')" \
> >       -e "install.packages('DEoptimR')" \
> >       -e "install.packages('http://cran.r-project.org/src/contrib/Archive/robustbase/robustbase_0.90-2.tar.gz', repos=NULL, type='source')" \
> >       -e "install.packages('dendextend')" \
> >       -e "install.packages('RColorBrewer')" \
> >       -e "install.packages('locfit')" \
> >       -e "install.packages('KernSmooth')" \
> >       -e "install.packages('BiocManager')" \
> >       -e "source('http://bioconductor.org/biocLite.R')" \
> >       -e "biocLite('Biobase')" \
> >       -e "biocLite('BioGenerics')" \
> >       -e "biocLite('BiocParallel')" \
> >       -e "biocLite('SingleCellExperiment')" \
> >       -e "biocLite('GenomeInfoDb')" \
> >       -e "biocLite('GenomeInfgoDbData')" \
> >       -e "biocLite('DESeq')" \
> >       -e "biocLite('DESeq2')" \
> >       -e "BiocManager::install(c('scater', 'scran'))" \
> >       -e "library('devtools')" \
> >       -e "install_github('IMB-Computational-Genomics-Lab/ascend', ref = 'devel')" \
> >       && rm -rf /tmp/downloaded_packages
> > ```
> {: .solution}
> 
> > ## Samtools compiled in the image
> >
> > ```bash
> > FROM ubuntu:18.04
> > 
> > # Image metadata
> > LABEL maintainer="john.doe@nowhere.com"
> > 
> > # Define version as build variable
> > ARG SAM_VER="1.9"
> > 
> > # Good practice variables
> > ENV DEBIAN_FRONTEND="noninteractive"
> > ENV LANG="C.UTF-8" LC_ALL="C.UTF-8"
> > 
> > # Install apt dependencies
> > RUN apt-get update && \
> >     apt-get -y install \
> >       gcc \
> >       libbz2-dev \
> >       libcurl4-openssl-dev \
> >       liblzma-dev \
> >       libncurses5-dev \
> >       libncursesw5-dev \
> >       make \
> >       perl \
> >       tar \
> >       vim \
> >       wget \
> >       zlib1g-dev \
> >     && apt-get clean all && \
> >     apt-get purge && \
> >     rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*
> > 
> > # Build samtools
> > RUN mkdir /build && \
> >     cd /build && \
> >     wget https://github.com/samtools/samtools/releases/download/${SAM_VER}/samtools-${SAM_VER}.tar.bz2 && \
> >     tar -vxjf samtools-${SAM_VER}.tar.bz2 && \
> >     cd samtools-${SAM_VER} && \
> >     ./configure --prefix=/apps && \
> >     make && \
> >     make install && \
> >     cd htslib-${SAM_VER} && \
> >     make && \
> >     make install && \
> >     cd / && \
> >     rm -rf /build
> > 
> > # Define PATH variable
> > ENV PATH=/apps/bin:$PATH
> > 
> > # Default command to be bash
> > CMD ["/bin/bash"]
> > ```
> {: .solution}
{: .challenge}
