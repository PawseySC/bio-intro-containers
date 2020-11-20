---
title: "Basics of Singularity"
teaching: 10
exercises: 20
questions:
objectives:
- Download container images
- Run commands from inside a container
- Discuss what are the most popular image registries
keypoints:
- Singularity can run both Singularity and Docker container images
- Execute commands in containers with `singularity exec`
- Open a shell in a container with `singularity shell`
- Download a container image in a selected location with `singularity pull`
- You should not use the `latest` tag, as it may limit workflow reproducibility
- The most commonly used registries are Docker Hub, Quay, Biocontainers and Nvidia GPU Cloud
---


### Get ready for the hands-on

Before we start, let us ensure we have got the required files to run the tutorials.

If you haven't done it already, download the following Github repo, then `cd` into it.  **NOTE**: if the `/data` directory is unavailable, use home `~` instead.

```bash
cd /data
git clone https://github.com/PawseySC/abacbs-containers
cd abacbs-containers
```


### Singularity: a container engine for HPC

[Singularity](https://sylabs.io/singularity/) is developed and maintained by [Sylabs](https://sylabs.io), and was designed from scratch as a container engine for HPC applications, which is clearly reflected in some of its main features:

* *unprivileged* runtime: Singularity containers do not require the user to hold root privileges to run (the Singularity executable needs to be installed and owned by *root*, though);

* *integration*, rather than *isolation*, by default: same user as host, same shell variables inherited by host, current directory bind mounted, communication ports available; as a result, launching a container requires a much simpler syntax than Docker;

* interface with job schedulers, such as *Slurm* or *PBS*;

* ability to run MPI enabled containers using host libraries;

* native execution of GPU enabled containers;

* unfortunately, *root* privileges are required to build container images: users can build images on their personal laptops or workstations, on the cloud, or via a Remote Build service.

Notably, Singularity is able to run containers from images that were built both with Singularity and with Docker-compatible engines.

This tutorial assumes Singularity version 3.0 or higher. Version **3.5.0 or higher** is recommended as it offers a smoother, more bug-free experience.


### Executing a simple command in a Docker container

For these first exercises, we're going to use a plain *Ubuntu* container image from [**Docker Hub**](https://hub.docker.com), *i.e.* the main registry for Docker containers.  This is a small and quick image to download, and will allow us to get to know how containers work by using typical Linux commands.  

From within the tutorial directory, let us cd into `exercises/singularity`:

```bash
cd exercises/singularity
```

Running a command is done by means of `singularity exec`:

```bash
singularity exec docker://ubuntu:16.04 cat /etc/os-release
```

```output
INFO:    Converting OCI blobs to SIF format
INFO:    Starting build...
Getting image source signatures
Copying blob sha256:22e816666fd6516bccd19765947232debc14a5baf2418b2202fd67b3807b6b91
 25.45 MiB / 25.45 MiB [====================================================] 1s
Copying blob sha256:079b6d2a1e53c648abc48222c63809de745146c2ee8322a1b9e93703318290d6
 34.54 KiB / 34.54 KiB [====================================================] 0s
Copying blob sha256:11048ebae90883c19c9b20f003d5dd2f5bbf5b48556dabf06c8ea5c871c8debe
 849 B / 849 B [============================================================] 0s
Copying blob sha256:c58094023a2e61ef9388e283026c5d6a4b6ff6d10d4f626e866d38f061e79bb9
 162 B / 162 B [============================================================] 0s
Copying config sha256:6cd71496ca4e0cb2f834ca21c9b2110b258e9cdf09be47b54172ebbcf8232d3d
 2.42 KiB / 2.42 KiB [======================================================] 0s
Writing manifest to image destination
Storing signatures
INFO:    Creating SIF file...
INFO:    Build complete: /data/singularity/.singularity/cache/oci-tmp/a7b8b7b33e44b123d7f997bd4d3d0a59fafc63e203d17efedf09ff3f6f516152/ubuntu_16.04.sif

NAME="Ubuntu"
VERSION="16.04.6 LTS (Xenial Xerus)"
ID=ubuntu
ID_LIKE=debian
PRETTY_NAME="Ubuntu 16.04.6 LTS"
VERSION_ID="16.04"
HOME_URL="http://www.ubuntu.com/"
SUPPORT_URL="http://help.ubuntu.com/"
BUG_REPORT_URL="http://bugs.launchpad.net/ubuntu/"
VERSION_CODENAME=xenial
UBUNTU_CODENAME=xenial
```

Here is what Singularity has just done:

* downloaded the various layers making up the Docker image;
* assembled them into a single SIF image file;
* stored the SIF file into the default cache directory;
* instantiated a container from that image;
* executed the command `cat /etc/os-release`.

Container images have a **name** and a **tag**, in this case `ubuntu` and `16.04`.  The tag can be omitted, in which case Singularity will default to a tag named `latest`.


> ## Using the *latest* tag
>
> The practice of using the `latest` tag can be handy for quick typing, but is dangerous when it comes to reproducibility of your workflow, as under the hood the *latest* tag could point to different images over time.
{: .callout}


Here Singularity pulled the image from an online image registry, as represented in this example by the prefix `docker://`, that corresponds to Docker Hub.  
Images in there are organised as: `<repository>/<name>:<tag>`.  In the case of the Ubuntu image, the repository is the default `library`, which can be omitted.  
The full specification is used in the next example:

```bash
singularity exec docker://library/ubuntu:16.04 echo "Hello World"
```

```output
Hello World
```

Here we are also experiencing image caching in action: the output has no more mention of the image being downloaded.

Various versions of the same image may exist in an image registry, typically in the form of distinct tags.  For instance, let's now use the latest version of Ubuntu, 20.04:

```bash
singularity exec docker://ubuntu:20.04 cat /etc/os-release
```

```output
[..]
NAME="Ubuntu"
VERSION="20.04 LTS (Focal Fossa)"
[..]
```


### Open up an interactive shell

Sometimes it can be useful to open a shell inside a container, rather than to execute commands, *e.g.* to inspect its contents.

Achieve this by using `singularity shell`:

```bash
singularity shell docker://ubuntu:16.04
```

```output
Singularity ubuntu_16.04.sif:/home/ubuntu/singularity-containers/exercises/singularity>
```

Remember to type `exit`, or hit `Ctrl-D`, when you're done!


### Download and use images via SIF file names

All examples so far have identified container images using their registry name specification, *e.g.* `docker://ubuntu:16.04` or similar.

An alternative option to handle images is to download them to a known location, and then refer to their full directory path and file name.

Let's use `singularity pull` to save the image to a specified path (output might differ depending on the Singularity version you use):

```bash
singularity pull docker://ubuntu:16.04
```

By default, the image is saved in the current directory:

```bash
ls
```

```output
ubuntu_16.04.sif
```

Then you can use this image file by:

```bash
singularity exec ./ubuntu_16.04.sif echo "Hello World"
```

```output
Hello World
```

You can specify the storage location with the `--dir` flag:

```bash
mkdir -p sif_lib
singularity pull --dir ~/path/to/sif/lib docker://library/ubuntu:16.04
```

Being able to specify download locations allows you to keep the local set of images organised and tidy, by making use of a directory tree.  It also allows for easy sharing of images within your team in a shared resource.  In general, you will need to specify the location of the image upon execution, *e.g.* by defining a dedicated variable:

```bash
export image="/data/abacbs-containers/exercises/singularity/ubuntu_16.04.sif"
singularity exec $image echo "Hello Again"
```

```output
Hello Again
```


> ## Contextual help on Singularity commands
>
> Use `singularity help`, optionally followed by a command name, to print help information on features and options.
{: .callout}


### Popular registries (*aka* image libraries)

At the time of writing, **Docker Hub** hosts a much wider selection of container images than **Sylabs Cloud**.  This includes Linux distributions, Python and R deployments, as well as a big variety of applications.

Bioinformaticians should keep in mind another container registry, [Quay](https://quay.io) by Red Hat, that hosts thousands of applications in this domain of science.  These mostly come out of the [Biocontainers](https://biocontainers.pro) project, that aims to provide automated container builds of all of the packages made available through [Bioconda](https://bioconda.github.io).

Nvidia maintains the [Nvidia GPU Cloud (NGC)](https://ngc.nvidia.com), hosting an increasing number of containerised applications optimised to run on GPUs.

Right now, the official [Sylabs Cloud Library](https://cloud.sylabs.io) does not contain a large number of images.

