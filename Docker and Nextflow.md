# Developer environment using Docker, Singularity and Nextflow

Post-GWAS analysis code running environment developing using Docker and Nextflow

https://www.nextflow.io/docs/latest/getstarted.html



# 1. Set up Docker

## Install Docker

### Install Docker or Docker toolbox in windows

Docker Ref: https://docs.docker.com/docker-for-windows/

Docker toolbox, Ref: https://docs.docker.com/toolbox/toolbox_install_windows/

```bash
$ dockerd
$ docker version
```

> Client:
>  Version:           19.03.1
>  API version:       1.40
>  Go version:        go1.12.7
>  Git commit:        74b1e89e8a
>  Built:             Wed Jul 31 15:18:18 2019
>  OS/Arch:           windows/amd64
>  Experimental:      false
>
> Server: Docker Engine - Community
>  Engine:
>   Version:          19.03.5
>   API version:      1.40 (minimum version 1.12)
>   Go version:       go1.12.12
>   Git commit:       633a0ea838
>   Built:            Wed Nov 13 07:28:45 2019
>   OS/Arch:          linux/amd64
>   Experimental:     false
>  containerd:
>   Version:          v1.2.10
>   GitCommit:        b34a5c8af56e510852c35414db4c1f4fa6172339
>  runc:
>   Version:          1.0.0-rc8+dev
>   GitCommit:        3e425f80a8c931f88e6d94a8c831b9d5aa481657
>  docker-init:
>   Version:          0.18.0
>   GitCommit:        fec3683



### Install Docker in Linux

Ref: https://phoenixnap.com/kb/how-to-install-docker-on-ubuntu-18-04

```bash
sudo apt update -y
sudo apt remove docker docker-engine docker.io # remove old version
sudo apt install docker.io

sudo systemctl start docker # Start docker
sudo systemctl enable docker # Automate docker

docker version
```



### (Option) Troubleshoot the Docker daemon

Ref: https://docs.docker.com/config/daemon/



## Set up Dockerfile

Ref: https://velog.io/@lazysoul/Docker-Basic-Usage

Docker 개발환경 팁, Ref: https://jhb.kr/368

Run Docker Quickstart Terminal:



### Build image from Dockerfile

Docker document, Ref: https://docs.docker.com/engine/reference/commandline/image_save/

Build docker image, Ref: https://subicura.com/2017/02/10/docker-guide-for-beginners-create-image-and-deploy.html

Install R packages, Ref: https://stackoverflow.com/questions/45289764/install-r-packages-using-docker-file

Cannot install R packages error, Ref: https://stackoverflow.com/questions/41601496/cannot-install-r-packages-in-docker-image

Avoiding user interaction with tzdata, Ref: https://askubuntu.com/questions/909277/avoiding-user-interaction-with-tzdata-when-installing-certbot-in-a-docker-contai



Generate a `Dockerfile` with this contents:

```txt
FROM ubuntu:latest
MAINTAINER kisudsoe@gmail.com
RUN apt -y update
RUN apt -y upgrade
RUN apt install openjdk-11-jre-headless
...
```



### Link between Github and Docker Hub projects

Automated build by Docker Hub, Ref: https://shine-yeolmae.tistory.com/34

Automated build set up, Ref: https://docs.docker.com/docker-hub/builds/



### (Optional) Local build of image

Ref: https://www.nextflow.io/

* Make sure java 8 (or later, upto 11) is installed which is required for Nextflow.

```bash
sudo apt update
sudo apt -y upgrade
sudo apt install openjdk-11-jre-headless
```

Save image

```bash
cd "C:\Users\kisud\OneDrive\Suh's Lab\Postgwas_v3\docker"
docker build -t kisudsoe/postgwas-env:latest .
```

> ...
>
> Removing intermediate container fbdb5efcc0d4
>  ---> 67febf44105e
> Successfully built 67febf44105e
> Successfully tagged kisudsoe/postgwas-env:latest
> SECURITY WARNING: You are building a Docker image from Windows against a non-Windows Docker host. All files and directories added to build context will have '-rwxr-xr-x' permissions. It is recommended to double check and reset permissions for sensitive files and directories.

```bash
docker images
```

> REPOSITORY              TAG                 IMAGE ID            CREATED             SIZE
> kisudsoe/postgwas-env   latest              67febf44105e        8 minutes ago       1.2GB
> alpine                  latest              a24bb4013296        2 months ago        5.57MB
> ubuntu                  latest              4e5021d210f6        4 months ago        64.2MB
> ubuntu-upstart          latest              b28219773b9b        4 years ago         253MB



### Test Docker Environment image

```bash
docker run -it kisudsoe/postgwas-env:latest
```

```bash
java --version
bedtools
R --version
exit
```



### (Solution) Error: biomaRt is not installed.

biomaRt (Bioconductor)

* Issue for "XML had Non Zero Exit Status" error: https://stackoverflow.com/questions/20671814/non-zero-exit-status-r-3-0-1-xml-and-rcurl
* `$ sudo apt install libcurl4-openssl-dev libxml2-dev`
* Issue for "openssl had non-zero exit status" error: https://github.com/rocker-org/rocker/issues/124
* `$ sudo apt install libssl-dev`
* Issue for "Biobase had non-zero exit status" error: https://stackoverflow.com/questions/56241007/non-zero-exit-status-r-3-6-0-biobase
* `> Sys.setenv(R_INSTALL_STAGED = FALSE)`
* Issue for install XML library: https://stackoverflow.com/questions/26042751/cannot-install-package-xml-to-r 
* `> install.packages("XML", repos = "http://www.omegahat.net/R")`

```R
BiocManager::install("biomaRt")
```



### Upload the Environment image to Docker Hub

```bash
$ docker login
$ docker tag kisudsoe/postgwas-env:latest kisudsoe/postgwas-env:1

$ docker push kisudsoe/postgwas-env:latest
$ docker push kisudsoe/postgwas-env:1
```



# 2. Update Docker image with source codes



## (Optional) Pull latest postgwas-env image

```bash
#docker rmi -f 67febf44105e
$ docker pull kisudsoe/postgwas-env
```

> Using default tag: latest
> latest: Pulling from kisudsoe/postgwas-env
> 3ff22d22a855: Pull complete                                                                                                                e7cb79d19722: Pull complete                                                                                                                323d0d660b6a: Pull complete                                                                                                                b7f616834fd0: Pull complete                                                                                                                b6163fe5c976: Pull complete                                                                                                                63cf5a70f41d: Pull complete                                                                                                                505fb5139b6a: Pull complete                                                                                                                7ebebcf3b64b: Pull complete                                                                                                                
> Digest: sha256:a0d4a148dc8d511afc0efb7f78e78efb6e1070aa7d9b688b31d7b90bfd27a117
> Status: Downloaded newer image for kisudsoe/postgwas-env:latest
> docker.io/kisudsoe/postgwas-env:latest



### Mount local volume to container

Mount volume, Ref: https://headsigned.com/posts/mounting-docker-volumes-with-docker-toolbox-for-windows/

In Oracle VM Virtual Box Manager -> Settings (default machine) -> Shared Folders -> Add new shared folder.

* Don't forget to check `Auto-mount` and `Make Permanent` options

Transient Folders: `C:\Users\kisud\OneDrive\Suh's Lab\Postgwas_v3` path named as `Postgwas_v3`

```CMD
docker-machine restart # Run first time only

# option 1: Run postgwas-env image
docker run -it -v "C:\Users\kisud\OneDrive\Suh's Lab\Postgwas_v3:/postgwas/data" ^
  kisudsoe/postgwas-env:latest ^
  /bin/bash

# option 2: Run postgwas image
docker run -it -v "C:\Users\kisud\OneDrive\Suh's Lab\Postgwas_v3:/postgwas/data" ^
  kisudsoe/postgwas:latest ^
  /bin/bash
```



## Update source codes

Run image with interactive mode:

```CMD
docker run -it -v "C:\Users\kisud\OneDrive\Suh's Lab\Postgwas_v3:/data" ^
  kisudsoe/postgwas:latest /bin/bash
```

In the container, copy the updated source to image

```CMD
cp -r /data/src /data/postgwas-exe.r /
```

From another CMD window, save image using Container ID, Ref: https://galid1.tistory.com/323

```CMD
docker ps
```

> CONTAINER ID        IMAGE                      COMMAND             CREATED             STATUS              PORTS               NAMES
> ae9fe688d36d        kisudsoe/postgwas:latest   "/bin/bash"         18 minutes ago      Up 18 minutes                           silly_kowalevski

```CMD
docker stop 3956737172ce
docker commit -a "jjy" 3956737172ce kisudsoe/postgwas:latest
```



### Test the Docker image

```CMD
docker pull kisudsoe/postgwas # This is optional
docker run -it kisudsoe/postgwas:latest
```

```bash
container$ cd /postgwas
container$ Rscript postgwas-exe.r --help
```

> Version: 2020-06-26
>
> Usage:
>     Rscript postgwas-exe.r --gwas <Function> --base <base file(s)> --out <out folder> [options:--p.criteria]
>     Rscript postgwas-exe.r --ldlink <Function> --base <base file(s)> --out <out folder> [options:--popul --r2 --dprime]
>     Rscript postgwas-exe.r --dbdown <Function> --out <out folder> [option:--hg]
>     Rscript postgwas-exe.r --dbfilt <Function> --base <base file/folder> --out <out folder> [option:--hg]
>     Rsciprt postgwas-exe.r --bedtools <Function> --base <base SNP BED file> --out <out folder>
>     ...
>
>
> Function calls:
>     --gwas    A function for GWAS Catalog data.
>     --ldlink  A function for LDlink data.
>     --dbdown  A function for downloading databases.
>     --dbfilt  A function for filtering data.
>     --dbvenn  A function for venn analysis.
>     --dbgene  A function for gene analysis.
>     --dbcomp  A function for PCA analysis for datasets.
>
> Global arguments:
>     --base    <Base input file path>
>               This is mendatory.
>     --out     <Out folder path>
>               This is mendatory.
>     --debug   <default: FALSE>
>               TRUE/FALSE: Rich description for debugging.
>               This is optional.
>
> Running functions with "--help" argument prints [Function] usage information.



### Upload new Docker image to Docker Hub

```CMD
# v4-2020-08-15
docker tag kisudsoe/postgwas:latest kisudsoe/postgwas:6
docker login
docker push kisudsoe/postgwas:latest
docker push kisudsoe/postgwas:6
```

Ref Docker: mirnylab/distiller_env



## * Incorporating SNPsea in Docker image

Install SNPsea requirements, Ref: https://snpsea.readthedocs.io/en/latest/installation.html

Run image:

```CMD
docker run -it -v "C:\Users\kisud\OneDrive\Suh's Lab\Postgwas_v3:/data" ^
  kisudsoe/postgwas:latest /bin/bash
```

Save updated image

```CMD
docker ps
```

> CONTAINER ID        IMAGE                      COMMAND             CREATED             STATUS              PORTS               NAMES
> 8fbcde8f38d8        kisudsoe/postgwas:latest   "/bin/bash"         3 hours ago         Up 3 hours                              festive_black

```CMD
docker stop 88530ebe5cfb
docker commit -a "jjy" 88530ebe5cfb kisudsoe/postgwas:latest
```

Upload image to Docker Hub server

```CMD
# v5-2020-08-14
docker tag kisudsoe/postgwas:latest kisudsoe/postgwas:5
docker login
docker push kisudsoe/postgwas:latest
docker push kisudsoe/postgwas:5
```



# 3. Run PostGWAS pipeline

## Start up the Postgwas image

Download Docker image `Postgwas:latest`.

```bash
$ docker pull kisudsoe/postgwas
```

Run the image

```bash
$ docker run -it -v \
	"/Postgwas_v3:/postgwas/data" \
	kisudsoe/postgwas:latest \
	/bin/bash
```

Source code located in `/postgwas` directory and local directory mounted as `/source`



## Run Postgwas pipeline

This process run in the image: `root@106653e2d631:/#`

Check running the function.

```bash
Rscript postgwas/postgwas-exe.r --help
```







# Archives

# #. Run Docker image through Singularity

Ref: https://sylabs.io/guides/3.6/user-guide/

### Install Singularity 3.6 in Windows

Ref: https://sylabs.io/guides/3.6/admin-guide/installation.html#installation-on-windows-or-mac

Install system dependencies

* Git for Windows https://gitforwindows.org/
* VirtualBox for Windows https://www.virtualbox.org/wiki/Downloads
* Vagrant for Windows https://www.vagrantup.com/downloads.html
* Vagrant Manager for Windows https://www.vagrantmanager.com/downloads/



### Singularity Vagrant Box

Run **Git Bash** and create and enter a directory to be used with your Vagrant VM.

```bash
mkdir singularity && \
    cd singularity
```

(Option) Destroy old VM and delete the Vagrantfile

```bash
vagrant destroy && \
    rm Vagrantfile
```

Bring up the Virtual Machine

```bash
export VM=sylabs/singularity-3.5-ubuntu-bionic64 && \
    vagrant init $VM && \
    vagrant up && \
    vagrant ssh
```

Check Singularity version and run:

```bash
vagrant@vagrant:~$ singularity version
vagrant@vagrant:~$ singularity help
```

> 3.5.1
>
> Linux container platform optimized for High Performance Computing (HPC) and
> Enterprise Performance Computing (EPC)
>
> Usage:
>   singularity [global options...]
>
> Description:
>   Singularity containers provide an application virtualization layer enabling
>   mobility of compute via both application and environment portability. With
>   Singularity one is capable of building a root file system that runs on any
>   other Linux system where Singularity is installed.
>
> Options:
>   -c, --config string   specify a configuration file (for root or
>                         unprivileged installation only) (default
>                         "/usr/local/etc/singularity/singularity.conf")
>   -d, --debug           print debugging information (highest verbosity)
>   -h, --help            help for singularity
>       --nocolor         print without color output (default False)
>   -q, --quiet           suppress normal output
>   -s, --silent          only print errors
>   -v, --verbose         print additional information
>
> Available Commands:
>   build       Build a Singularity image
>   cache       Manage the local cache
>   capability  Manage Linux capabilities for users and groups
>   config      Manage various singularity configuration (root user only)
>   delete      Deletes requested image from the library
>   exec        Run a command within a container
>   help        Help about any command
>   inspect     Show metadata for an image
>   instance    Manage containers running as services
>   key         Manage OpenPGP keys
>   oci         Manage OCI containers
>   plugin      Manage Singularity plugins
>   pull        Pull an image from a URI
>   push        Upload image to the provided URI
>   remote      Manage singularity remote endpoints
>   run         Run the user-defined default command within a container
>   run-help    Show the user-defined help for an image
>   search      Search a Container Library for images
>   shell       Run a shell within a container
>   sif         siftool is a program for Singularity Image Format (SIF) file manipulation
>   sign        Attach digital signature(s) to an image
>   test        Run the user-defined tests within a container
>   verify      Verify cryptographic signatures attached to an image
>   version     Show the version for Singularity
>
> Examples:
>   $ singularity help <command> [<subcommand>]
>   $ singularity help build
>   $ singularity help instance start



### Get Docker image through Singularity

```bash
vagrant@vagrant:~$ singularity pull docker://kisudsoe/postgwas:latest
#singularity pull docker://godlovedc/lolcow
```

> INFO:    Converting OCI blobs to SIF format
> INFO:    Starting build...
> Getting image source signatures
> Copying blob 5bed26d33875 done
> Copying blob f11b29a9c730 done
> Copying blob 930bda195c84 done
> Copying blob 78bf9a5ad49e done
> Copying blob 8bd8fcbbf4e9 done
> Copying blob 793bc4eb372b done
> Copying blob d96f1d58ad07 done
> Copying blob d7b3260324e1 done
> Copying blob 3c3c295cf442 done
> Copying config 45f5420dcb done
> Writing manifest to image destination
> Storing signatures
> 2020/08/01 20:34:24  info unpack layer: sha256:5bed26d33875e6da1d9ff9a1054c5fef3bbeb22ee979e14b72acf72528de007b
> 2020/08/01 20:34:25  info unpack layer: sha256:f11b29a9c7306674a9479158c1b4259938af11b97359d9ac02030cc1095e9ed1
> 2020/08/01 20:34:25  info unpack layer: sha256:930bda195c84cf132344bf38edcad255317382f910503fef234a9ce3bff0f4dd
> 2020/08/01 20:34:25  info unpack layer: sha256:78bf9a5ad49e4ae42a83f4995ade4efc096f78fd38299cf05bc041e8cdda2a36
> 2020/08/01 20:34:25  info unpack layer: sha256:8bd8fcbbf4e90950638093835b7fafde39e9e121ddb531a39cd2f37b31c8f1aa
> 2020/08/01 20:34:27  info unpack layer: sha256:793bc4eb372b927d6ee5f26b59bf55740e50d720bc2406e8b66f99515fc834de
> 2020/08/01 20:34:27  info unpack layer: sha256:d96f1d58ad07cfa6022ec91fda0cf3e2b3d4a2f3072657435319c78f71277f5e
> 2020/08/01 20:34:44  info unpack layer: sha256:d7b3260324e17c77442de149b3cdf709f0369161e2f1dde8fc9570724d328dd5
> 2020/08/01 20:34:51  info unpack layer: sha256:3c3c295cf442ad39a5d27f1f9d08f8a682c36eaffa02177dac035fcf543b5fe7
> INFO:    Creating SIF file...



### Run docker image through Singularity

Is it really need?



### (Optional) Pull latest Postgwas image from Docker Hub

```bash
docker pull kisudsoe/postgwas
```

> Using default tag: latest
> latest: Pulling from kisudsoe/postgwas
> 5bed26d33875: Pull complete                                                                                                                f11b29a9c730: Pull complete                                                                                                                930bda195c84: Pull complete                                                                                                                78bf9a5ad49e: Pull complete                                                                                                                8bd8fcbbf4e9: Pull complete                                                                                                                793bc4eb372b: Pull complete                                                                                                                d96f1d58ad07: Pull complete                                                                                                                d7b3260324e1: Pull complete                                                                                                                3c3c295cf442: Pull complete                                                                                                                
> Digest: sha256:e308f2d4b32a6b05d4923441832f707b6d43f20f556529a28da775f8bc5ef034
> Status: Downloaded newer image for kisudsoe/postgwas:latest
> docker.io/kisudsoe/postgwas:latest



# #. Set up Nextflow



## Install Nextflow in local

* Be aware of target folder path. If a special character or a space is included, it causes an Error.

```bash
$ java -version
$ sudo curl -s https://get.nextflow.io | bash
$ mv nextflow "/mnt/c//Users/kisud/OneDrive/Suh's Lab/Postgwas_v3"
$ cd "/mnt/c//Users/kisud/OneDrive/Suh's Lab/Postgwas_v3"
$ ./nextflow run hello
```

> N E X T F L O W  ~  version 20.07.1
> Pulling nextflow-io/hello ...
> downloaded from https://github.com/nextflow-io/hello.git
> Launching `nextflow-io/hello` [focused_darwin] - revision: 96eb04d6a4 [master]
> executor >  local (4)
> [4e/9eeb86] process > sayHello (2) [100%] 4 of 4 ✔
> Hola world!
>
> Hello world!
>
> Bonjour world!
>
> Ciao world!



## Nextflow tutorials

Ref: https://www.nextflow.io/docs/latest/getstarted.html

Nextflow scripting language is an extension of the Groovy programming language.

* Ref: https://www.nextflow.io/docs/latest/script.html



### Tutorial 1

This is a code in `nextflow/tutorial.nf`:

```java
#!/usr/bin/env nextflow

params.str = 'Hello world!'

process splitLetters {

    output:
    file 'chunk_*' into letters

    """
    printf '${params.str}' | split -b 6 - chunk_
    """
}


process convertToUpper {

    input:
    file x from letters.flatten()

    output:
    stdout result

    """
    cat $x | tr '[a-z]' '[A-Z]'
    """
}

result.view { it.trim() }
```

To run this code, there are seveal commands like below:

```bash
./nextflow run tutorial.nf
./nextflow run tutorial.nf -resume
./nextflow run tutorial.nf --str 'Bonjour le monde'
```



### Tutorial 2

Processes are executed independently and are isolated from each other, i.e. they do not share a common (writable) state. The only way they can communicate is via asynchronous FIFO queues, called *channels* in Nextflow.

Any process can define one or more channels as *input* and *output*.

```java
// Script parameters
params.query = "/some/data/sample.fa"
params.db = "/some/path/pdb"

db = file(params.db)
query_ch = Channel.fromPath(params.query)

process blastSearch {
    input:
    file query from query_ch

    output:
    file "top_hits.txt" into top_hits_ch

    """
    blastp -db $db -query $query -outfmt 6 > blast_result
    cat blast_result | head -n 10 | cut -f 2 > top_hits.txt
    """
}

process extractTopHits {
    input:
    file top_hits from top_hits_ch

    output:
    file "sequences.txt" into sequences_ch

    """
    blastdbcmd -db $db -entry_batch $top_hits > sequences.txt
    """
}
```



## Using Docker with Nextflow

Ref: https://www.nextflow.io/blog/2016/docker-and-nextflow.html



### Nextflow for Docker

Ref: https://www.nextflow.io/docs/latest/docker.html

Ref: https://www.docker.com/

It is possible to specify a different Docker image for each process definition in your pipeline script.

```java
process foo {
  container 'image_name_1'

  '''
  do this
  '''
}

process bar {
  container 'image_name_2'

  '''
  do that
  '''
}
```

Alternatively, the same containers definitions can be provided by using the `nextflow.config` file as shown below:

```java
process {
    withName:foo {
        container = 'image_name_1'
    }
    withName:bar {
        container = 'image_name_2'
    }
}
docker {
    enabled = true
}
```



### Configuration: Scope *docker*

The `docker` configuration scope controls how [Docker](https://www.docker.com/) containers are executed by Nextflow.

The following settings are available:

| Name          | Description                                                  |
| ------------- | ------------------------------------------------------------ |
| enabled       | Turn this flag to `true` to enable Docker execution (default: `false`). |
| envWhitelist  | Comma separated list of environment variable names to be included in the container environment. |
| legacy        | Uses command line options removed since version 1.10.x (default: `false`). |
| sudo          | Executes Docker run command as `sudo` (default: `false`).    |
| tty           | Allocates a pseudo-tty (default: `false`).                   |
| temp          | Mounts a path of your choice as the `/tmp` directory in the container. Use the special value `auto` to create a temporary directory each time a container is created. |
| remove        | Clean-up the container after the execution (default: `true`). For details see: [https://docs.docker.com/engine/reference/run/#clean-up—rm](https://docs.docker.com/engine/reference/run/#clean-up---rm) . |
| runOptions    | This attribute can be used to provide any extra command line options supported by the `docker run` command. For details see: https://docs.docker.com/engine/reference/run/ . |
| registry      | The registry from where Docker images are pulled. It should be only used to specify a private registry server. It should NOT include the protocol prefix i.e. `http://`. |
| fixOwnership  | Fixes ownership of files created by the docker container.    |
| engineOptions | This attribute can be used to provide any option supported by the Docker engine i.e. `docker [OPTIONS]`. |
| mountFlags    | Add the specified flags to the volume mounts e.g. mountFlags = ‘ro,Z’ |





# #. Run PostGWAS pipeline

Ref: https://github.com/kisudsoe/PostGWAS-tools

Ref: https://github.com/mirnylab/distiller-nf

Ref: http://github.com/nextflow-io/hello



## Clone Github repository

Ref: https://www.nextflow.io/docs/latest/sharing.html#publishing-your-pipeline

* Download latest version of PostGWAS-tools

```bash
$ mkdir pe-5e8
$ ./nextflow clone kisudsoe/PostGWAS-tools ./pe-5e8
```



## Run nextflow pipeline

```bash
$ ./nextflow run pe-5e8
```

