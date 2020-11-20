---
title: "Building images with Docker"
teaching: 15
exercises: 15
questions:
objectives:
- Discuss the pros&cons of building with Singularity vs Docker
- Learn the basic syntax of a Dockerfile
- Build and share an image with Docker
- Convert a Docker image into the Singularity format
keypoints:
- Specify the build starting image with `FROM`
- Execute shell commands with `RUN`
- Declare environment variables with `ENV`
- Build images starting from a `Dockerfile` using `docker build`
- Push images to a web registry using `docker push`
- Convert a Docker image by using `singularity pull docker-daemon:`
- Share Singularity SIF images as any other file
---


### Why using Docker for container builds?

Over the whole of this tutorial we're proposing Singularity as the principal tool to handle and run containerised applications in HPC.  Do we really need to extend our toolkit to include [Docker](https://hub.docker.com/search/?type=edition&offering=community)?

To better inform an answer to this question, here are some of the advantages when building with one or the other tool.

#### Singularity
* Single file image, can be handled as any other file and shared easily;
* Unambiguous container usage modes, via distinct keywords: `exec`, `shell`, `run`, `instance` (see episode on web applications);
* Powerful ways of defining shell environment (see [Sylabs docs on Environment and Metadata](https://sylabs.io/guides/3.5/user-guide/environment_and_metadata.html));
* Ability to sign/verify images for improved security.

#### Docker
* Compatibility: image format can be run by all existing container engines, including Singularity;
* Layered image format allows caching, for reduced build time during prototyping and development (however shipping these images requires more effort).

Note how, at present, both tools require root privileges for building, implying that this step cannot be performed on HPC, and requires a dedicated machine instead.

Although Singularity builds offer some interesting advantages, there's a single item that right now makes using Docker preferred in most situations.  
It's **compatibility**.  Build with Docker, and you'll know the resulting image can be run from every container engine anywhere in the world.


### Building a container image with Docker

This example is adapted from this well crafted [Singularity Tutorial](https://github.com/ArangoGutierrez/Singularity-tutorial).  

Let's `cd` into the relevant demo directory:

```bash
cd /data/abacbs-containers/exercises/lolcow_docker
```

Let us also start building the image with `docker build` straight away.  Meanwhile, we'll bring the discussion on.

```bash
sudo docker build -t lolcow:21Nov20 .
```

First, note that we're typing `sudo`, because Docker needs to be executed with administrative privileges.  If your Linux user belongs to the "docker" group (this is the case in the machines for this tutorial), you can skip typing `sudo`, however the execution will still be privileged.

At the end of the line, `.` is the location of the build context (*i.e.* the directory for the Dockerfile).  

The `-t` flag is used to specify the image name (compulsory) and tag (optional).  
Any lowercase alphanumeric string can be used as image name; here we've used `lolcow`.  The image tag (following the colon) can be optionally used to maintain a set of different image versions on Docker Hub, and is a key feature in enabling reproducibility of your computations through containers; here we've used `21Nov20`.

If you have a Docker Hub account, Adding the prefix `<Your account>/` to the image name allows to push the built image to your Docker Hub registry for sharing (see below).  
Then, the complete format for the image name looks like: `<Your Docker Hub account ^>/<Image name>:<Image tag ^>`. `^`These are optional.

This is the output of our build:

```output
Sending build context to Docker daemon  2.048kB
Step 1/6 : FROM ubuntu:18.04
 ---> 2eb2d388e1a2
Step 2/6 : LABEL maintainer="Pawsey Supercomputing Centre"
 ---> Running in 6b3a2ca5a738
Removing intermediate container 6b3a2ca5a738
 ---> 8b97ddfd477c
Step 3/6 : LABEL version="v0.0.1"
 ---> Running in cc349d1718b6
Removing intermediate container cc349d1718b6
 ---> b7d793a1c012
Step 4/6 : RUN apt-get -y update &&   apt-get -y install fortune cowsay lolcat

[..]

Removing intermediate container 69028dcddc81
 ---> c61504ed1b78
Step 5/6 : ENV PATH=/usr/games:$PATH
 ---> Running in 1108c6ca6db4
Removing intermediate container 1108c6ca6db4
 ---> c8032059066d
Step 6/6 : CMD fortune | cowsay | lolcat
 ---> Running in c19994af1c27
Removing intermediate container c19994af1c27
 ---> 4d20bd8fcf3b
Successfully built 4d20bd8fcf3b
Successfully tagged lolcow:21Nov20
```


### A Dockerfile recipe

Now let's have a look at the `Dockerfile` recipe for this image:

```bash
cat Dockerfile
```

```source
FROM ubuntu:18.04
  
LABEL maintainer="Pawsey Supercomputing Centre"
LABEL version="v0.0.1"

RUN apt-get -y update && \
  apt-get -y install fortune cowsay lolcat

ENV PATH=/usr/games:$PATH

CMD fortune | cowsay | lolcat
```

The directory where the the Dockerfile is stored is the so called the Docker **build context**.  Docker will include files in this directory in the build process and in the final image.  As a by-product, this will make the build process longer and the image larger, so that we want to include only those strictly required for the build, even none if possible.

Let's comment on the Docker instructions that appear in this Dockerfile.

The first line in the Dockerfile is `FROM ubuntu:18.04`.  This tells Docker which base image to use to perform the build.

Lines starting with the instruction `RUN` are the most important.  These are the sequence of commands to be executed to install packages in the image, in essence the same commands you would use for installation in a Linux box.  Here we are ensuring we have an up-to-date list of packages, and then we are installing three Linux utilities.

The instruction `ENV` allows to set up environment variables that need to be defined at runtime rather than at build time.  Here the `PATH` needs to be updated to reflect the location of the three utilities that we installed with the `RUN` instructions.  **DO NOT** use `RUN export <..>` to this end, as the variable will be lost after the `RUN` step is completed.

The `LABEL` instruction is optionally used to add metadata to the container image.  

The final instruction, `CMD`, specifies the default command to be executed with the container, in case no other command is provided.  The most typical choice is the shell, `/bin/bash`, as specific commands are typically provided as arguments to the container execution line; here we're using a different command just to show something fancy later on.

Although we are not using them in this Dockerfile, other useful instructions are `ADD` or `COPY`, like in: 

```source
ADD <src-file> <dst-file>
```

These are used to copy files from the host, *i.e.* <src-file>, inside the container in the destination <dst-file>.

More information on the Dockerfile syntax can be found at the [Dockerfile reference](https://docs.docker.com/engine/reference/builder/).


### List local Docker images

We know that Docker container images are not single files, but rather adopt a multi-layered format.  To keep things tidy, Docker stores images and their layers in a hidden directory, under its own control.  
To get the list of available images, including the ones you built, use `docker images`:

```bash
sudo docker images
```

```output
REPOSITORY                          TAG                                IMAGE ID            CREATED             SIZE
lolcow                              21Nov20                            4d20bd8fcf3b        11 minutes ago      181MB
```


> ### Bonus: What if you need to debug the build process?
> 
> Quite often devising a recipe to install software involves a certain deal of trial and error.  
> 
> If you need to inspect a Docker container during build, you can open an interactive shell session this way:
> 
> ```bash
> sudo docker run --rm -it ubuntu:18.04 bash
> ```
> 
> ```output
> root@dd1ca993f4ad:/#
> ```
> 
> Here, `-it` keeps the container standard input open and allocates a terminal; `--rm` does some clean up when closing the session.  
> 
> Note that Docker containers, unlike Singularity containers, are writable.  So during an interactive sessions you can even trial software installations.  However, edits are ephemeral, *i.e.* you lose them when you close the container.
> 
> When you're done, type `exit`, or hit `Ctrl-D`, to leave the interactive shell.
{: .callout}


### Converting and running the image with Singularity

We've created a container image using Docker.  Now, how can we use it with Singularity?  
As seen in a previous episode, the Singularity `pull` command can take care of converting a Docker image for its own usage.

If your Docker-equipped machine also comes with Singularity, you can grab the image from the local image repo, using the `docker-daemon:` prefix, and convert it into the SIF format.  Note that, due to a current bug in Singularity, you will need to ditch the double slashes `//`:

```bash
singularity pull docker-daemon:lolcow:21Nov20
```

Does it work?  Let's try the command `fortune`:

```bash
singularity exec lolcow_21Nov20.sif fortune
```

```output
Whenever you find that you are on the side of the majority, it is time
to reform.
		-- Mark Twain
```

Now, try and run this pipe of commands: `bash -c 'fortune | cowsay | lolcat'`.

```bash
singularity exec lolcow_21Nov20.sif bash -c 'fortune | cowsay | lolcat'
```

You will get something similar to this, hopefully just more colourful:

```output
 _______________________________________
/ Have a place for everything and keep  \
| the thing somewhere else; this is not |
| advice, it is merely custom.          |
|                                       |
\ -- Mark Twain                         /
 ---------------------------------------
        \   ^__^
         \  (oo)\_______
            (__)\       )\/\
                ||----w |
                ||     ||
```

Hey, we've just containerised a cow that cites novelists!  

And just as a little trick, because we used that weird `CMD` definition in the Dockerfile, we can make the cow talk also via:

```bash
./lolcow_21Nov.sif
```

```output
 ________________________________________
/ Don't go surfing in South Dakota for a \
\ while.                                 /
 ----------------------------------------
        \   ^__^
         \  (oo)\_______
            (__)\       )\/\
                ||----w |
                ||     ||
```


### Sharing the Docker image using Docker Hub

If you have a (free) Docker Hub account you can use it to share your image in Docker format. 

First you must login to Docker:

```bash
sudo docker login
```

Then, create a second tag for the image, that includes your Docker Account.  To this end we need to use `docker tag`:

```bash
sudo docker tag lolcow:21Nov20 <your-dockerhub-account>/lolcow:21Nov20
```

Finally, we can push the image to the Docker Hub web registry:

```bash
sudo docker push <your-dockerhub-account>/lolcow:21Nov20
```

```output
The push refers to repository [docker.io/marcodelapierre/lolcow]
9d2959e72647: Pushed 
317d47a452af: Pushed 
e0b3afb09dc3: Mounted from library/ubuntu 
6c01b5a53aac: Mounted from library/ubuntu 
2c6ac8e5063e: Mounted from library/ubuntu 
cc967c529ced: Mounted from library/ubuntu 
21Nov20: digest: sha256:295c5695e2b05f6123bc2d8669ec7b66e45df5000ab9fc45ce3566ae3c0d839e size: 1571
```

Your image is now publicly available for anyone to pull, including using Singularity.


### Sharing the Singularity SIF file

An image in Singularity SIF format is just a file, so... you can transfer it across systems using Linux command line utilities like `scp` or `rsync`, or even graphical applications such as `Filezilla`.  
Just remember that images can be quite large, typically ranging from tens of MBs up to several GBs.  The *lolcow* image we created is about 70 MB, but for instance a typical *RStudio* image is well above 1 GB.

If you want to keep the images publicly available for the wider community, the most compatible way is to use the Docker image format as described above (`docker push` to Docker Hub).



### Useful base images

At the time of writing, [Docker Hub](https://hub.docker.com) is the most popular web registry for general purpose container images.  Therefore all images mentioned below are hosted in this registry.

#### R
The [Rocker Project](https://www.rocker-project.org) maintains a number of good R base images.  Of particular relevance is [rocker/tidyverse](https://hub.docker.com/r/rocker/tidyverse), which embeds the basic R distribution, an RStudio web-server installation and the [tidyverse](https://www.tidyverse.org) collection of packages for data science.

Other more basic images are [rocker/r-ver](https://hub.docker.com/r/rocker/r-ver) (R only) and [rocker/rstudio](https://hub.docker.com/r/rocker/rstudio) (R + RStudio).

#### Python
[python](https://hub.docker.com/_/python) hosts the official Python images.  Different versions are available for some OS flavours.  Smaller base images have tags ending with `-slim`.

[continuumio/miniconda3](https://hub.docker.com/r/continuumio/miniconda3) are images provided by the maintainers of the [Anaconda](https://anaconda.org) project.  They ship with Python 3, as well as `pip` and `conda` to install and manage packages.

If you need interactive Jupyter Notebooks, [Jupyter Docker Stacks](https://jupyter-docker-stacks.readthedocs.io/en/latest/) maintain a series of dedicated container images.  Among others, there is the base SciPy image [jupyter/scipy-notebook](https://hub.docker.com/r/jupyter/scipy-notebook), the data science image [jupyter/datascience-notebook](https://hub.docker.com/r/jupyter/datascience-notebook), and the machine learning image [jupyter/tensorflow-notebook](https://hub.docker.com/r/jupyter/tensorflow-notebook).

#### CUDA
[nvidia/cuda](https://hub.docker.com/r/nvidia/cuda) has images to build GPU enabled applications.  There are different image types for different needs.  Tags containing `runtime` are suitable for binary applications that are ready to run; if you need to compile GPU code, pick tags containing `devel` instead.  Different OS flavours are available, too.

#### MPI
As you can see in the episode on MPI applications, when containerising this type of software the MPI libraries in the image need to be ABI compatible with the MPI libraries in the host.  The Pawsey Supercomputing Centre maintains some dedicated base images at [pawsey/mpich-base](https://hub.docker.com/r/pawsey/mpich-base), for building images that will run on our HPC systems.

