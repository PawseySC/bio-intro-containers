---
title: "Share files and variables with the host machine"
teaching: 10
exercises: 20
questions:
objectives:
- Mount host directories in a container
- Pass specific variables to the container
keypoints:
- By default Singularity mounts the host current directory, and uses it as the container working directory
- Map additional host directories in the containers with the flag `-B`, or the variable SINGULARITY_BINDPATH
- Avoid mounting the `$HOME` directory, to better protect your sensitive data in the host
- By default Singularity passes all host variables to the container
- If you need to, isolate the container shell environment with the flag `-e`
- Pass specific shell variables to containers by prefixing them with SINGULARITYENV_
---


### Access to directories in the host machine

Let's start and `cd` into the root demo directory:

```bash
cd /data/abacbs-containers/exercises
```

What directories can we access from the container?  
First, let us assess what the content of the root directory `/` looks like from outside *vs* inside the container, to highlight the fact that a container runs on his own distinct filesystem.  From the host:

```bash
$ ls /
```

```output
bin  boot  dev  etc  home  lib  lib64  media  mnt  opt  proc  root  run  sbin  scratch  shared  srv  sys  tmp  usr  var
```

And from the container:

```bash
$ singularity exec docker://ubuntu:18.04 ls /
```

```output
bin  boot  data  dev  environment  etc	home  lib  lib64  media  mnt  opt  proc  root  run  sbin  singularity  srv  sys  tmp  usr  var
```

We can see they differ in content, *i.e.* they're distinct directory structures.

Now, in which directory is the container running?  For reference, this is the host:

```bash
pwd
```

```output
/home/ubuntu/singularity-containers/exercises
```

And now the container:
```bash
singularity exec docker://ubuntu:18.04 pwd
```

```output
/home/ubuntu/singularity-containers/exercises
```

By default, the container runs from the host working directory!

Now, we can see the contents of the current directory from inside the container:

 ```bash
 singularity exec docker://ubuntu:18.04 ls
 ```

 ```output
 blast_1        blast_2        build_examples database_blast lolcow_docker  nextflow       pipeline       singularity
 ```

This is why the BLAST exercise in the previous session worked!  
However, how about other directories in the host?  For instance, let us inspect `../_episodes`:

```bash
singularity exec docker://ubuntu:18.04 ls ../_episodes
```

```output
ls: cannot access '../_episodes': No such file or directory
```

Host directories that are external to the current directory are not visible!  How can we fix this?  Read on...


> ## By the way, can we write inside a container?
> 
> Try and create a file called `example` in the container root directory.  (**Hint**: run `touch /example` inside the container).
> 
> > ## Solution
> > 
> > ```bash
> > $ singularity exec docker://ubuntu:18.04 touch /example
> > ```
> > 
> > ```output
> > touch: cannot touch '/example': Read-only file system
> > ```
> > 
> > We have just learn something more on containers: by default, they are **read-only**.  How can we get a container to write files then?  Read on...
> {: .solution}
{: .challenge}


To summarise what we've learnt in the previous examples, we may say that a container ships an application and its dependencies by encapsulating them in an isolated, read-only filesystem.  
In order for a container to access directories from the host filesystem (and write files), one needs to explicitly bind mount them.  The main exception here is the current work directory, which is bind mounted by default.


### Bind mounting host directories

Singularity has the runtime flag `--bind`, `-B` in short, to mount host directories.

There is a long syntax, which allows to map the host dir onto a container dir with a different name/path, `-B hostdir:containerdir`.  
There is also a short syntax, that just mounts the dir using the same name and path: `-B hostdir`.

Let's use the latter syntax to mount the whole of `/data` into the container and re-run `ls`.

```bash
singularity exec -B /data docker://ubuntu:18.04 ls ../_episodes
```

```output
11-containers-intro.md  13-blast.md             22-pipeline.md          32-build-examples.md
12-singularity-intro.md 21-share-host.md        31-build-docker.md      41-workflow-engines.md
```

Also, we can write files in a host dir which has been bind mounted in the container:

```bash
singularity exec -B /data docker://ubuntu:18.04 touch ../_episodes/example
singularity exec -B /data docker://ubuntu:18.04 ls ../_episodes/example
```

```output
../_episodes/example
```

Now we are talking!

If you need to mount multiple directories, you can either repeat the `-B` flag multiple times, or use a comma-separated list of paths, *i.e.*

```bash
-B dir1,dir2,dir3
```

Also, if you want to keep the runtime command compact, you can equivalently specify directories to be bind mounted using the environment variable `SINGULARITY_BINDPATH`:

```bash
export SINGULARITY_BINDPATH="/data"
```


> ## Default bind mounted directories
> 
> HPC sites will typically default to bind mount their key filesystems at container runtime, so that you won't need to take any action to access files on those locations from inside a container.
> 
> However, knowing the bind mounting syntax can come very handy if you're running in a custom system, such as a virtual machine on a cloud service.  Most likely, no defaults will be set there, and you will have to configure host directory bind mounting by yourself.
{: .callout}


> ## Mounting $HOME
>
> Depending on the site configuration of Singularity, user home directories might or might not be mounted into containers by default.  
> We do recommend that you **avoid mounting home** whenever possible, to avoid sharing potentially sensitive data, such as SSH keys, with the container, especially if exposing it to the public through a web service.  In addition, it's also handy to avoid sharing configuration files located in the home, that may unexpectedly impact the behaviour of applications in the container.
>
> If you need to share data inside the container home, you might just mount that specific file/directory, *e.g.*
>
> ```bash
> -B $HOME/.local
> ```
>
> Or, if you want a full fledged home, you might define an alternative host directory to act as your container home, as in
>
> ```bash
> -B /path/to/fake/home:$HOME
> ```
>
> Finally, you should also **avoid running a container from your host home**, otherwise this will be bind mounted as it is the current working directory.
{: .callout}


### How about sharing environment variables with the host?

By default, shell variables are inherited in the container from the host:

```bash
export HELLO=world
singularity exec docker://ubuntu:18.04 bash -c 'echo $HELLO'
```

```output
world
```

There might be situations where you need to isolate the shell environment of the container.  For instance, this is typically the case for Python containers, as host Python variables may affect the Python runtime in the container.  
To isolate the shell environment, you can use the flag `-e`, or `--cleanenv`:  

```bash
export HELLO=world
singularity exec -e docker://ubuntu:18.04 bash -c 'echo $HELLO'
```

```output

```

If you need to pass only specific variables to the container, that might or might not be defined in the host, you can define variables that start with `SINGULARITYENV_`.  This prefix will be automatically trimmed in the container:

```bash
export SINGULARITYENV_BYE=bye
singularity exec -e docker://ubuntu:18.04 bash -c 'echo $BYE'
```

```output
bye
```

