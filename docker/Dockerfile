FROM ubuntu:latest

MAINTAINER Suhlab-SKim

RUN apt -y -qq update && \
    apt -y -qq upgrade

# Set locale
RUN DEBIAN_FRONTEND=noninteractive apt -y install tzdata

# Install dependencies
RUN apt -y -qq install \
		openjdk-11-jre-headless \
		r-base \
		wget \
		bedtools \
		libcurl4-openssl-dev libxml2-dev \
		libssl-dev

# Install R packages
RUN R -e "install.packages(c('dplyr','data.table','eulerr','circlize','LDlinkR','RSQLite','reshape','ggplot2','plyr','RSQLite'), dependencies=T, repos='http://cran.us.r-project.org/')" && \
	R -e "if (!requireNamespace('BiocManager',quietly=T)) install.packages('BiocManager')" && \
	R -e "Sys.setenv(R_INSTALL_STAGED = FALSE)" && \
    R -e "install.packages('XML', repos = 'http://www.omegahat.net/R')" && \
	R -e "BiocManager::install('biomaRt')" && \
    R -e "BiocManager::install('limma')"

# Install PostGWAS-tools
ADD https://github.com/kisudsoe/PostGWAS-tools/archive/0.1.tar.gz
RUN unzip 0.1.tar.gz

# Version number contained in image
ADD VERSION .
