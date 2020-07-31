# Docker and Nextflow

Post-GWAS analysis code running environment developing using Docker and Nextflow

https://www.nextflow.io/docs/latest/getstarted.html



# Requirements

## Install Docker or Docker toolbox in windows

Docker Ref: https://docs.docker.com/docker-for-windows/

Docker toolbox, Ref: https://docs.docker.com/toolbox/toolbox_install_windows/



### Troubleshoot the Docker daemon

Ref: https://docs.docker.com/config/daemon/

```bash
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



## Install Docker in Linux

Ref: https://phoenixnap.com/kb/how-to-install-docker-on-ubuntu-18-04

```bash
sudo apt update -y
sudo apt remove docker docker-engine docker.io # remove old version
sudo apt install docker.io

sudo systemctl start docker # Start docker
sudo systemctl enable docker # Automate docker

docker version
```



## Set up Docker environment

### Set up Docker Environment

Ref: https://velog.io/@lazysoul/Docker-Basic-Usage

Docker 개발환경 팁, Ref: https://jhb.kr/368

Run Docker Quickstart Terminal:

```bash
$ docker images
$ docker run -v "C:/Users/kisud/oneDrive/Suh's Lab/Postgwas_v3:/postgwas" alpine ls /postgwas
$ docker run -it ubuntu
```

> REPOSITORY          TAG                 IMAGE ID            CREATED             SIZE
> alpine              latest              a24bb4013296        2 months ago        5.57MB
> ubuntu              latest              4e5021d210f6        4 months ago        64.2MB
> ubuntu-upstart      latest              b28219773b9b        4 years ago         253MB



### Install Java 8 or later, upto 11

Ref: https://www.nextflow.io/

* Make sure java 8 (or later, upto 11) is installed.

```bash
sudo apt update
sudo apt -y upgrade
sudo apt install openjdk-11-jre-headless
```



### Build image from Dockerfile

Docker document, Ref: https://docs.docker.com/engine/reference/commandline/image_save/

Build docker image, Ref: https://subicura.com/2017/02/10/docker-guide-for-beginners-create-image-and-deploy.html

Install R packages, Ref: https://stackoverflow.com/questions/45289764/install-r-packages-using-docker-file

Cannot install R packages error, Ref: https://stackoverflow.com/questions/41601496/cannot-install-r-packages-in-docker-image

Avoiding user interaction with tzdata, Ref: https://askubuntu.com/questions/909277/avoiding-user-interaction-with-tzdata-when-installing-certbot-in-a-docker-contai



Generate a `Dockerfile` with this contents:

```bash
FROM ubuntu:latest
MAINTAINER kisudsoe@gmail.com
RUN apt -y update
RUN apt -y upgrade
RUN apt install openjdk-11-jre-headless
...
```

Automated build by Docker Hub, Ref: https://shine-yeolmae.tistory.com/34

* Link github and docker hub
* 

(Optional) Local build of image

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



### Upload the Environment image to Docker Hub

```bash
$ docker login
$ docker tag kisudsoe/postgwas-env:latest kisudsoe/postgwas-env:1

$ docker push kisudsoe/postgwas-env:latest
$ docker push kisudsoe/postgwas-env:1
```



## Update Docker with source codes

### Pull latest Docker environment

```bash
#docker rmi -f 67febf44105e
$ docker pull kisudsoe/postgwas-env
```



### Run postgwas-env:latest

```bash
$ docker run -it kisudsoe/postgwas-env:latest
$ mkdir /postgwas
```



### Update source codes

Mount volume, Ref: https://headsigned.com/posts/mounting-docker-volumes-with-docker-toolbox-for-windows/

In Oracle VM Virtual Box Manager -> Settings (default machine) -> Shared Folders -> Add new shared folder.

* Don't forget to check `Auto-mount` and `Make Permanent` options

Transient Folders: `C:\Users\kisud\OneDrive\Suh's Lab\Postgwas_v3` path named as `Postgwas_v3`

```bash
$ docker-machine restart
$ docker run -it -v \
  "/Postgwas_v3:/source" \
  kisudsoe/postgwas-env:latest \
  /bin/bash
```

Move updated source files

```bash
container$ cp -r /source/src /source/postgwas-exe.r /postgwas/
```



### Save as new image

Save image, Ref: https://galid1.tistory.com/323

```bash
$ docker ps
```

> CONTAINER ID        IMAGE                          COMMAND             CREATED             STATUS              PORTS               NAMES
> bb10195afd1c        kisudsoe/postgwas-env:latest   "/bin/bash"         12 minutes ago      Up 12 minutes                           unruffled_dewdney

```bash
$ docker stop bb10195afd1c
$ docker commit -a "jjy" bb10195afd1c kisudsoe/postgwas:latest
```



### Test Docker image

```bash
$ docker run -it kisudsoe/postgwas:latest
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



### Upload the new Docker image to Docker Hub

```bash
$ docker tag kisudsoe/postgwas:latest kisudsoe/postgwas:1

$ docker login
$ docker push kisudsoe/postgwas:latest
$ docker push kisudsoe/postgwas:1
```

Docker: mirnylab/distiller_env



## Install Nextflow in local

* Be aware of target folder path. If a special character or a space is included, it causes an Error.

```bash
$ java -version
$ sudo curl -s https://get.nextflow.io | bash
$ mv nextflow target/folder/path/
$ cd target/folder/path/
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



# Nextflow tutorials



## Tutorials

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





# Run PostGWAS pipeline

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

