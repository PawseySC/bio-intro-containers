---
title: "Session 3: building containers"
teaching: 10
exercises: 50
questions:
objectives:
---


### Goals

In this session, we're going to build together three container images.

In particular we'll cover:
* one RStudio image;
* one Conda image;
* one image where we're compiling the software ourselves.

As we've seen in the webinars, here the key tool is not Singularity but Docker, due to the higher compatibility and overwhelming popularity of its container image format.  
We're going to touch on relevant aspects of writing a Dockerfile:
* choosing an appropriate base image and declaring it with `FROM`;
* using `RUN` to execute commands;
* using `ENV` to define environment variables;
* other useful commands;
* best practices.

Then, we're going to explore the typical steps that are required in the creation of a container image:
* building an image with `docker build`;
* testing and debugging it with `docker run`;
* assigning a new name to it with `docker tag`;
* sharing it on a public registry with `docker push`;
* converting it offline for usage with Singularity, with `singularity pull docker-daemon:`.

To brush up on building container images with Docker, you can refer to the webinar episode [Building images with Docker](https://pawseysc.github.io/singularity-containers/22-build-docker/index.html).

Writing a Dockerfile for a container is an art, which you can refine over time with practice.  
We don't mean to be exhaustive in this session; instead, we hope to provide you with the basic and most common commands, as well as some good practices.


### GUIDED - Let's write an R Dockerfile

The first exercise in this session is to write a little Dockerfile for the R package `ggtree`.  This package provides functionalities to represent phylogenetic trees using R, by building on top of the [Tidyverse](https://www.tidyverse.org) collection of data science packages, and in particular the plotting library `ggplot2`.  
The end result will actually be the RStudio image that is used in the session on graphical applications.

The first step in designing a Dockerfile is to choose a base image, that is the starting point for our container.

The best place to find useful base image is the [Docker Hub](https://hub.docker.com) online registry.  Here is a non-comprehensive list of potentially useful base images (the version tags are not necessarily the most recent ones):

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


> ## Use at least two terminal tabs
>
> We advise you to open at least two terminal tabs, and connect to the VM from both of them.  In this way, you can use one to edit files, and one to execute commands, thus making your workflow more efficient.
{: .callout}


Now, let's cd into the appropriate directory:

```
$ cd /data/containers-bioinformatics-workshop
$ cd exercises/build/r-ggtree
```
{: .bash}
<!-- $ export WORK=$(pwd) -->

And then use a text editor to create a blank `Dockerfile`.


> ## Pick the base image
> 
> Considering the characteristic that we stated above for the `ggtree` package, which is the simplest base image you would choose from this list?
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
> > ```
> > FROM rocker/tidyverse:3.6.1
> > ```
> > {: .source}
> {: .solution}
{: .challenge}


> ## A non-container question
> 
> The package `ggtree` is part of the BioConductor project.  So, from inside an R console we could install it by using the command `BiocManager::install("ggtree")`.  However, here we need Docker to execute this command from a bash shell.
> 
> If you are an R user, do you know how you can execute the R command above from the shell?  No worries if you don't, just have a look at the solution.
> 
> > ## Solution
> > 
> > ```
> > R -e 'BiocManager::install("ggtree")'
> > ```
> > {: .source}
> {: .solution}
{: .challenge}


> ## Command to install an R package
> 
> We're almost there with our first Dockerfile... now we just need to embed the shell command above in the Dockerfile, by using the appropriate Docker instruction.
> 
> > ## Solution
> > 
> > ```
> > FROM rocker/tidyverse:3.6.1
> > 
> > RUN R -e 'BiocManager::install("ggtree")'
> > ```
> > {: .source}
> > 
> > Very often (albeit not always) preparing Dockerfiles for R images looks as simple as this.  Other times, you will also need to install other packages such as pre-requisites.  And of course things can even get more complicated than this.
> {: .solution}
{: .challenge}


> ## Build the image
> 
> Now it's the time for building.  
> Let' use the appropriate `docker` syntax in the shell to build an image called `ggtree:2.0.4`; remember we're running from the directory where the Dockerfile is (`.`). 
> 
> > ## Solution
> > 
> > ```
> > $ sudo docker build -t ggtree:2.0.4 .
> > ```
> > {: .bash}
> > 
> > It will only take a couple of minutes to build, as most required R packages are already provided by the base Tidyverse image.
> {: .solution}
> 
> Note how here we provided the information about the package version; normally we would be able to find it out ourselves after the first build, by inspecting the R installation in the container.  Let's do it.
{: .challenge}


> ## Test the image
> 
> In order to test that the image we built actually contains a working version of `ggtree`, we're going to execute the following shell command to print the version of this package, from inside the container:
> 
> ```
> $ R -e 'packageVersion("ggtree")'
> ```
> {: .bash}
> 
> Here's something new for you: if you want to run a container with Docker, you need to execute `docker run`; then, similar to Singularity, we need to append the image name and the command.
> 
> > ## Solution
> > 
> > ```
> > $ sudo docker run ggtree:2.0.4 R -e 'packageVersion("ggtree")'
> > ```
> > {: .bash}
> > 
> > ```
> > R version 3.6.1 (2019-07-05) -- "Action of the Toes"
> > Copyright (C) 2019 The R Foundation for Statistical Computing
> > Platform: x86_64-pc-linux-gnu (64-bit)
> > 
> > [..]
> > 
> > > packageVersion("ggtree")
> > [1] ‘2.0.4’
> > > 
> > > 
> > ```
> > {: .output}
> > 
> > It has worked!  The container has got `ggtree` version 2.0.4.
> {: .solution}
{: .challenge}


### ROOM - Write a Dockerfile to compile a tool

This exercise is a bit more articulated.  You're going to build an image for a popular bioinformatics tool, *samtools*.  In the Dockerfile, you're compiling the tool yourself;  this will give the opportunity to touch on a set of important aspects.

Cd into the appropriate directory:

```
$ cd ../compile-samtools
```
{: .bash}

We're going to assume that you already know how to install samtools in a Ubuntu box.  Here is the corresponding bash installation script (there's a copy in your current directory, `install_samtools.sh`):

```
#!/bin/bash

# Install apt dependencies
apt-get update
apt-get -y install \
  gcc \
  libbz2-dev \
  libcurl4-openssl-dev \
  liblzma-dev \
  libncurses5-dev \
  libncursesw5-dev \
  make \
  perl \
  tar \
  vim \
  wget \
  zlib1g-dev

# Build samtools
mkdir /build
cd /build

wget https://github.com/samtools/samtools/releases/download/1.9/samtools-1.9.tar.bz2
tar -vxjf samtools-1.9.tar.bz2

cd samtools-1.9
./configure --prefix=/apps
make
make install

cd htslib-1.9
make
make install
```
{: .bash}

Here's a short description of each block in this script:
* install package dependencies using `apt`, the Ubuntu package manager; note how backslashes `\` followed by newlines are used to improve the readability of the code;
* create a build directory, `/build`;
* download and untar the source code;
* configure and compile samtools (note we're specifying the installation path to be `/apps`);
* compile htslib, a companion library.


> ## Which base image would you pick?
> 
> This installation recipe was tested with Ubuntu.  Which base image would you choose from the list at the beginning of this session?
> 
> > ## Solution
> > 
> > `ubuntu:18.04` seems suitable.
> {: .solution}
{: .challenge}


> ## Dockerfile 1: tool dependencies
> 
> Create a blank `Dockerfile` in the current directory, and insert the following, using the appropriate Docker instructions:
> * a line to define the base image;
> * two lines to execute the `apt` command for the tool prerequisites (see reference bash script above).
> 
> Remember you have a copy of the reference script in the work directory, `install_samtools.sh`, in case you need to copy/paste bits of it.
> 
> > ## Solution
> > 
> > ```
> > FROM ubuntu:18.04
> > 
> > RUN apt-get update
> > RUN apt-get -y install \
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
> >       zlib1g-dev
> > ```
> > {: .source}
> {: .solution}
> 
> In this exercise, you're going to use `docker build` several times as you grow the Dockerfile.  Go ahead now with the first build; you can call the image `sam:1`.
> 
> > ## Solution
> > 
> > ```
> > $ sudo docker build -t sam:1 .
> > ```
> > {: .bash}
> > 
> > This will take a couple of minutes, and generate a lot of output.  
> > 
> > Along the output, do you notice that there are references to *Steps*?  For instance at the very beginning:
> > 
> > ```
> > Sending build context to Docker daemon  13.82kB
> > Step 1/3 : FROM ubuntu:18.04
> >  ---> 775349758637
> > Step 2/3 : RUN apt-get update
> >  ---> Running in 1a102770efa6
> > ```
> > {: .output}
> > 
> > It seems like this build is made up of 3 steps: there's one for each Docker instruction in the Dockerfile!
> > 
> > The last line will look like:
> > 
> > ```
> > Successfully tagged sam:1
> > ```
> > {: .output}
> {: .solution}
{: .challenge}


> ## Dockerfile 2: samtools installation
> 
> Now grab all of the remaining commands from the reference bash script, and append them in the Dockerfile, each with its own Docker instruction.
> 
> > ## Solution
> > 
> > ```
> > 
> > [..]
> > 
> > RUN mkdir /build
> > RUN cd /build
> > RUN wget https://github.com/samtools/samtools/releases/download/1.9/samtools-1.9.tar.bz2
> > RUN tar -vxjf samtools-1.9.tar.bz2
> > RUN cd samtools-1.9
> > RUN ./configure --prefix=/apps
> > RUN make
> > RUN make install
> > RUN cd htslib-1.9
> > RUN make
> > RUN make install
> > ```
> > {: .source}
> {: .solution}
> 
> Now, run the build again, calling the image `sam:2`.  What happens?
> 
> > ## Solution
> > 
> > ```
> > $ sudo docker build -t sam:2 .
> > ```
> > {: .bash}
> > 
> > Well, there are a few things worth noting.
> > 
> > 1. The `apt` commands now execute instantaneously, thanks to Docker layer caching (results from the previous `docker build` are re-used):
> > 
> >     ```
> >     Step 1/14 : FROM ubuntu:18.04
> >      ---> 775349758637
> >     Step 2/14 : RUN apt-get update
> >      ---> Using cache
> >      ---> 14d7f8a970f4
> >     Step 3/14 : RUN apt-get -y install       gcc       libbz2-dev       libcurl4-openssl-dev       liblzma-dev       libncurses5-dev       libncursesw5-dev       make       perl       tar       vim       wget       zlib1g-dev
> >      ---> Using cache
> >      ---> 62984f3be676
> >     ```
> >     {: .output}
> > 
> > 2. The build now comprises 14 steps.
> > 
> > 3. The build stops with an error on step 9:
> > 
> >     ```
> >      ---> 6952729b0c9d
> >     Step 9/14 : RUN ./configure --prefix=/apps
> >      ---> Running in 00f36b1ce58c
> >     /bin/sh: 1: ./configure: not found
> >     The command '/bin/sh -c ./configure --prefix=/apps' returned a non-zero code: 127
> >     ```
> >     {: .output}
> >     
> >     The `configure` executable is supposed to be inside `/build/samtools-1.9`, but instead it is not found.  What's going on?
> {: .solution}
{: .challenge}


> ## Dockerfile 3: fixing the samtools installation
> 
> In order to understand what's going on, you'll now open an interactive shell in the last step of the image build, by using `docker run` and the flags for an interactive shell, `-it`.  This is a great practice to effectively debug the image building process.
> 
> What image shall you use?  Have a look at the error output that was generated.  At the end of each step, an alphanumeric image ID is printed, in lines starting with `--->`.  Retrieve the image ID just before the failing step (*9/14*); in the case of the sample output above the ID would be `6952729b0c9d`.  
> Now use this image ID to execute `bash`:
> 
> ```
> $ sudo docker run -it 6952729b0c9d bash
> ```
> {: .bash}
> 
> ```
> root@c8365ee497d9:/# 
> ```
> {:. output}
> 
> You're inside a container launched from an intermediate image along the build process.  
> Now investigate the content of `/build/` using `ls`.
> 
> > ## Solution
> > 
> > ```
> > # ls /build
> > ```
> > {: .bash}
> > 
> > ```
> > 
> > ```
> > {: .output}
> > 
> > It's empty! Mmh suspicious... 
> {: .solution}
> 
> Now check the root directory `/`.
> 
> > ## Solution
> > 
> > ```
> > # ls /
> > ```
> > {: .bash}
> > 
> > ```
> > bin  boot  build  dev  etc  home  lib  lib64  media  mnt  opt  proc  root  run  samtools-1.9  samtools-1.9.tar.bz2  sbin  srv  sys  tmp  usr  var
> > ```
> > {: .output}
> > 
> > Here you can see both the compressed source code, `samtools-1.9.tar.bz2`, and the uncompressed directory `samtools-1.9`.  
> > If you look back at the reference recipe, and to the latest Dockerfile, you would expect those to be under `/build` instead.
> > 
> > You can now close the container interactive shell by typing `exit`, or hitting `Ctrl-D`.
> > 
> > You're seeing here a behaviour that is characteristic of Docker builds: every time you issue a new `RUN` instruction, commands are executed from the root directory `/` and NOT from the last active directory, as it would be the case in a regular Linux shell.
> {: .solution}
> 
> Hot to fix this?  
> One good way is to pack bash commands, that need to execute from the same directory, within the same single `RUN` instruction.  You can concatenate them by using the Linux syntax `&&`; note that this instructs the shell (or Docker here) to execute the following command only if the previous one has ended successfully (with no error).  For instance:
> 
> ```
> RUN cd /build
> RUN wget https://github.com/samtools/samtools/releases/download/1.9/samtools-1.9.tar.bz2
> ```
> {: .source}
> 
> will become (also using `\` and newlines for clarity):
> 
> ```
> RUN cd /build && \
>     wget https://github.com/samtools/samtools/releases/download/1.9/samtools-1.9.tar.bz2
> ```
> {: .source}
> 
> Now go ahead and pack all of the related commands together in the same `RUN`.
> 
> > ## Solution
> > 
> > ```
> > 
> > [..]
> > 
> > RUN mkdir /build
> > RUN cd /build && \
> >     wget https://github.com/samtools/samtools/releases/download/1.9/samtools-1.9.tar.bz2 && \
> >     tar -vxjf samtools-1.9.tar.bz2 && \
> >     cd samtools-1.9 && \
> >     ./configure --prefix=/apps && \
> >     make && \
> >     make install && \
> >     cd htslib-1.9 && \
> >     make && \
> >     make install
> > ```
> > {: .source}
> {: .solution}
> 
> And now try and build again, using the image name `sam:3`.
> 
> > ## Solution
> > 
> > ```
> > $ sudo docker build -t sam:3 .
> > ```
> > {: .bash}
> > 
> > Depending on the extent of your grouping action, you will now have between 3 and 5 build steps.  
> > It's taking a few minutes, and it looks like it's compiling some code... in the end:
> > 
> > ```
> > Successfully built cbd9a61c748c
> > Successfully tagged sam:3
> > ```
> > {: .output}
> > 
> > It worked!  You successfully compiled samtools in the container image.
> {: .solution}
> 
> Now test that you can actually execute `samtools`, using `docker run`.
> 
> > ## Solution
> > 
> > ```
> > $ sudo docker run sam:3 samtools
> > ```
> > {: .bash}
> > 
> > ```
> > docker: Error response from daemon: OCI runtime create failed: container_linux.go:348: starting container process caused "exec: \"samtools\": executable file not found in $PATH": unknown.
> > ERRO[0000] error waiting for container: context canceled 
> > ```
> > {: .output}
> > 
> > There's an error: *executable file not found in $PATH*.  
> > What's going on now?  Let's go have a look.
> {: .solution}
{: .challenge}


> ## Dockerfile 4: locate samtools executables
> 
> Open again an interactive session with `docker run -it`, to inspect the image `sam:3`.
> 
> > ## Solution
> > 
> > ```
> > $ sudo docker run -it sam:3 bash
> > ```
> > {: .bash}
> {: .solution}
> 
> From the `./configure` line of the recipe and Dockerfile, you see that you asked to install the package under `/apps`.  
> Use `ls` to inspect that directory, *e.g.* to look for a subdirectory `bin`, which is a common location for executables.
> 
> > ## Solution
> > 
> > ```
> > # ls /apps
> > ```
> > {: .bash}
> > 
> > ```
> > bin  include  lib  share
> > ```
> > {: .output}
> > 
> > There's indeed a `bin` directory.
> > 
> > ```
> > # ls /apps/bin
> > ```
> > {: .bash}
> > 
> > ```
> > ace2sam  blast2sam.pl   export2sam.pl  interpolate_sam.pl  maq2sam-short  md5sum-lite  plot-bamstats  sam2vcf.pl  samtools.pl            soap2sam.pl  varfilter.py  wgsim_eval.pl
> > bgzip    bowtie2sam.pl  htsfile        maq2sam-long        md5fa          novo2sam.pl  psl2sam.pl     samtools    seq_cache_populate.pl  tabix        wgsim         zoom2sam.pl
> > ```
> > {: .output}
> > 
> > Bingo!  This directory contains a set of executables, including `samtools`.
> > 
> > You can close the interactive session by typing `quit` or hitting `Ctrl-D`.
> {: .solution}
> 
> So your last edit to the Dockerfile will be to add `/apps/bin` to the `$PATH` environment variable, which tells a Linux shell where to look for executables.  How shall you achieve this?
> 
> In bash, you would use `export PATH=/apps/bin:$PATH`, so you might be tempted to just go with `RUN export PATH=/apps/bin:$PATH`.  
> However, this won't work in a Dockerfile, as such a definition would just hold within the build step for that `RUN`, and then be lost at the next step (a bit like what happened with the current directory being reset to `/`).
> 
> Well, you might remember that Docker has an instruction to declare runtime environment variables, `ENV`.  Now use it to append the `PATH` specification in your Dockerfile.
> 
> > ## Solution
> > 
> > ```
> > 
> > [..]
> > 
> > ENV PATH=/apps/bin:$PATH
> > ```
> > {: .source}
> {: .solution}
> 
> And now for one last (?) build, with name `sam:4`...
> 
> > ## Solution
> > 
> > ```
> > $ sudo docker build -t sam:4 .
> > ```
> > {: .bash}
> > 
> > This is super fast, everything is cached, except for the variable definition.
> {: .solution}
> 
> Finally, try and `docker run` again the `samtools` command.
> 
> > ## Solution
> > 
> > ```
> > $ sudo docker run sam:4 samtools
> > ```
> > {: .bash}
> > 
> > ```
> > Program: samtools (Tools for alignments in the SAM format)
> > Version: 1.9 (using htslib 1.9)
> > 
> > Usage:   samtools <command> [options]
> > 
> > [..]
> > ```
> > {: .output}
> > 
> > Success!
> {: .solution}
{: .challenge}


If you had troubles, you can pick up the final Dockerfile from `solutions/Dockerfile.4`.


### ROOM - Bonus - Write a Conda Dockerfile

Go through this only **if time allows**.  

In this last exercise you're going to build another image for samtools; this time, though, you'll simply use the *Conda* package manager to install it.

First, cd into the appropriate directory, and then create a blank `Dockerfile` with a text editor:

```
$ cd ../conda-samtools
```
{: .bash}


> ## First bit: pick the base image
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
> > ```
> > FROM continuumio/miniconda3:4.8.2
> > ```
> > {: .source}
> {: .solution}
{: .challenge}


> ## Second bit: the command to install a Conda package
> 
> Samtools version 1.9 can be installed with Conda through the channel *bioconda*, so the shell syntax would be (`-y` is to confirm prompts):
> 
> ```
> conda install -y -c bioconda samtools=1.9
> ```
> {: .source}
> 
> With this information, complete the Dockerfile to install samtools using Conda.
> 
> > ## Solution
> > 
> > ```
> > FROM continuumio/miniconda3:4.8.2
> > 
> > RUN conda install -y -c bioconda samtools=1.9
> > ```
> > {: .source}
> > 
> > Similar to the case of R, Dockerfiles often turn out to be quite compact when using a `conda` base image.  Things are not always this easy, as for instance package version conflicts are common with Conda;  additional command lines might be required to work around them, or you might even need to install packages entirely manually.
> {: .solution}
> 
> In the interest of time, you're not going to build and test this Conda-based image.
{: .challenge}


<!-- 
> ## Now build the image...
> 
> Similar to the R exercise, what is the appropriate `docker` shell command to build an image called `samtools-conda:1.9`?
> 
> > ## Solution
> > 
> > ```
> > $ sudo docker build -t samtools-conda:1.9 .
> > ```
> > {: .bash}
> > 
> > The build should only take a couple of minutes.
> {: .solution}
{: .challenge}


> ## ...And then test it!
> 
> To check that the samtools installation in the container is working, use `docker run` to execute the command `samtools` from the image you just built.
> 
> > ## Solution
> > 
> > ```
> > $ sudo docker run samtools-conda:1.9 samtools
> > ```
> > {: .bash}
> > 
> > ```
> > Program: samtools (Tools for alignments in the SAM format)
> > Version: 1.9 (using htslib 1.9)
> > 
> > Usage:   samtools <command> [options]
> > 
> > [..]
> > ```
> > {: .output}
> {: .solution}
{: .challenge}
 -->


> ## ROOM - Questions for Reflection
> 
> These are a set of questions designed to kick off the discussion within your breakout room:
> 
> * Do you need to install R, Python or Conda packages to run your workflows?  Many of them?
> * Do you need to compile tools yourself?
> * Do you often need to reinstall the same tools, for instance on different machines?  Do you and your collaborators use similar tools?
> * Can you see the advantage in making the effort to install a tool in a container image only once, and then just downloading when needed?
> * Look back at what you went through in the previous exercises.  What was the most interesting learning outcome?  Which was the hardest bit to digest?
{: .discussion}


### GUIDED - Share and convert the image

In this session we've built three container images so far.  How can we share them publicly, and run them using our preferred container runtime, *i.e.* Singularity?  
Let's reply to these questions by using the latest image we built for the compiled version of samtools, `sam:4`.


#### Upload an image to Docker Hub

We need a free Docker Hub account to perform this task.  If you don't have one, just relax and watch, it won't be long.

First, we need to give the image a new name that includes our Docker Hub account; we'll call the image `test-samtools:1.9`.  We're using `docker tag` to this end:

```
$ sudo docker tag sam:4 <your-dockerhub-account>/test-samtools:1.9
```
{: .bash}

Then, we're ready to upload the newly renamed image, by means of `docker push`:

```
$ sudo docker push <your-dockerhub-account>/test-samtools:1.9
```
{: .bash}

This will take a few minutes...

```
The push refers to repository [docker.io/marcodelapierre/test-samtools]
[..]
1.9: digest: sha256:48e92fe7287326ca5e582ad562666977be790f884c0ed1b7af5d5284559d2400 size: 1995
```
{: .output}


#### Convert into Singularity Image Format

If we're on the machine where we built the Docker image, we can convert it offline, using `singularity pull docker-daemon:`:

```
$ singularity pull docker-daemon:<your-dockerhub-account>/test-samtools:1.9
```
{: .bash}

```
INFO:    Converting OCI blobs to SIF format
INFO:    Starting build...
Getting image source signatures
[..]
Writing manifest to image destination
Storing signatures
[..]
INFO:    Creating SIF file...
INFO:    Build complete: test-samtools_1.9.sif
```
{: .output}

In alternative, if we pushed the Docker image onto Docker hub, we can `pull` and convert it from there:

```
$ singularity pull docker://<your-dockerhub-account>/test-samtools:1.9
```
{: .bash}

Once the image is available in the Singularity format...

```
$ ls
```
{: .bash}

```
Dockerfile  best_practices  install_samtools.sh  solutions  test-samtools_1.9.sif
```
{: .output}

...we can test it!

```
$ singularity exec test-samtools_1.9.sif samtools
```
{: .bash}

```
Program: samtools (Tools for alignments in the SAM format)
Version: 1.9 (using htslib 1.9)

Usage:   samtools <command> [options]

[..]
```
{: .output}


### GUIDED - Best practices

In the `best_practices` subdirectory of your current work directory, you can find a final version of a `Dockerfile` for samtools, that implements all the best practices described below.  Here's a copy, too:

```
FROM ubuntu:18.04

# Image metadata
LABEL maintainer="john.doe@nowhere.com"

# Define version as build variable
ARG SAM_VER="1.9"

# Good practice variables
ENV DEBIAN_FRONTEND="noninteractive"
ENV LANG="C.UTF-8" LC_ALL="C.UTF-8"

# Install apt dependencies
RUN apt-get update && \
    apt-get -y install \
      gcc \
      libbz2-dev \
      libcurl4-openssl-dev \
      liblzma-dev \
      libncurses5-dev \
      libncursesw5-dev \
      make \
      perl \
      tar \
      vim \
      wget \
      zlib1g-dev \
    && apt-get clean all && \
    apt-get purge && \
    rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

# Build samtools
RUN mkdir /build && \
    cd /build && \
    wget https://github.com/samtools/samtools/releases/download/${SAM_VER}/samtools-${SAM_VER}.tar.bz2 && \
    tar -vxjf samtools-${SAM_VER}.tar.bz2 && \
    cd samtools-${SAM_VER} && \
    ./configure --prefix=/apps && \
    make && \
    make install && \
    cd htslib-${SAM_VER} && \
    make && \
    make install && \
    cd / && \
    rm -rf /build

# Define PATH variable
ENV PATH=/apps/bin:$PATH

# Default command to be bash
CMD ["/bin/bash"]
```
{: .source}


1. **Condensate commands into few `RUN` instructions**

    Grouping commands into single `RUN`s is a good practice when building images, as it reduces the number of image caching *layers*, and thus the total size of the image.  
    
    Why is it so?  
    At one extreme, imagine one `RUN` per command.  Any edit to the same files in the container or any deletion of unnecessary files would create a separate layer.  The final container would be larger, as it would keep a snapshot of each intermediate state, *e.g.* multiple versions of one file or even deleted files.  
    On the other extreme, concentrating every command in one `RUN` would minimise image size, as there would be one single layer for the whole set of files that makes up your containerised application.  As a disadvantage, readability of the Dockerfile would be reduced, and caching during build would become less effective.
    
    In the end, the best practice is to group together commands that relate to the same component of the installation.  Readability and caching are improved in this way, and you would not gain space anyway if you further grouped them, as these groups of commands operate on different files.
    
    In our samtools example, you can for instance have one `RUN` for the `apt` installation of dependencies, and one `RUN` for samtools.


2. **Clean the installation process**

    Installation files that are not required by the application runtime can be deleted, contributing to reduce the final image size (when coupled with practice `1.` above).  
    
    In our example with samtools, you would clean the `apt` installations with the following commands:
    
    ```
    RUN apt-get clean all && \
        apt-get purge && \
        rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*
    ```
    {: .source}
    
    and the samtools build files with:
    
    ```
    RUN cd / && \
        rm -rf /build
    ```
    {: .source}
    
    These two removal operations, when combined with the `RUN` grouping from best practice `1.`, will reduce the final size of your samtools image by 30%!
    
    Note how the following command can be useful in the context of a Conda installation:
    
    ```
    RUN /opt/conda/bin/conda clean -ay
    ```
    {: .source}


3. **Know and set some useful environment variables**

    Dockerfile installations are non-interactive by nature, *i.e.* no installer can ask you questions during the process.  
    You can define a variable prior to running any `apt` command, that informs a shell in Ubuntu or Debian that you are not able to interact, so that no questions will be asked and default values will be picked instead.
    
    This is the variable:
    
    ```
    ENV DEBIAN_FRONTEND="noninteractive"
    ```
    {: .source} 
    
    Another pair of useful variables, again to be put at the beginning of the Dockerfile, are:
    
    ```
    ENV LANG="C.UTF-8" LC_ALL="C.UTF-8"
    ```
    {: .source}
    
    These variables are used to specify the language localisation, or *locale*, to the value *C.UTF-8* in this case.  
    The *locale* specification impacts text rendering, time/date formats, monetary formats, and language-specific characters.  Leaving this undefined can result, for some programs, in warnings or in characters being displayed inappropriately (both at build and run time).  On the other end, setting the *locale* will avoid these issues, and should be considered as a best practice to enforce in any Dockerfile.


4. **Document your Dockerfile with labels**

    This is a way to provide information to third party users of your container image, to enable a more effective use.
    
    The typical label to add is the *maintainer* one, so that people can get in touch with them in case they are in trouble using the image:
    
    ```
    LABEL maintainer="john.doe@nowhere.com"
    ```
    {: .source}


5. **Think about the default command**

    The Docker instruction `CMD` can be used to set the default command that gets executed by `docker run <IMAGE>` (without arguments) or `singularity run <IMAGE>`.  
    If you don't specify it, the default will be the shell.  
    
    You might be tempted to set the application binary as default, *e.g.* samtools, but in the end this is not a very useful setting.  The `CMD` default gets overridden by any command/argument you will add to the `run` commands, so you will need to explicitly state the main command anyway.  Moreover, lots of packages come with multiple executables, in which case some will be left out anyway.
    
    In the end our suggestion for standalone containerised applications (such as samtools) is to set it to `/bin/bash`.  This is the default anyway, but documenting your choices is always good.  
    For Python and R containers, you might set it to the `python` and `R` interpreter, respectively, if you want.  
    Do not set `CMD` for `rocker/` RStudio and `jupyter/` Jupyter images; these come with articulated setups that allow to spawn the web servers, so your choice would be ignored anyway.


6. **Abstract package versions if you can**

    The `ARG` instruction can be used to define variables in the Dockerfile, that only exist within the build process.  This can be especially useful to specify package versions in a more general and flexible way.
    
    For instance, in the samtools Dockerfile you could define a variable for the tool version:
    
    ```
    ARG SAM_VER="1.9"
    ```
    {: .source}
    
    and then, substitute the explicit version `1.9` with `$SAM_VER` all along the commands that install samtools.
    
    A nice feature of `docker build` is that you can dynamically assign these `ARG` at build time using the `--build-arg` flag, *e.g.* if you want to install samtools `1.8`:
    
    ```
    $ sudo docker build -t sam:1.8 --build-arg SAM_VER=1.8 .
    ```
    {: .bash}
