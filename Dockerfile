FROM ubuntu:latest
WORKDIR /

MAINTAINER Suhlab-Seungsoo_Kim

RUN apt -y -qq update && \
	apt -y -qq upgrade

# Set locale & R v4.0 repository
RUN DEBIAN_FRONTEND=noninteractive apt -y install tzdata

# Install dependencies
RUN apt -y -qq install \
	openjdk-11-jre-headless \
	wget \
	bedtools \
	libcurl4-openssl-dev \
	libxml2-dev \
	libssl-dev \
	vim \
	software-properties-common \
	libcairo2-dev \
	libxt-dev \
	zsh
RUN apt-key adv --keyserver keyserver.ubuntu.com --recv-keys E298A3A825C0D65DFD57CBB651716619E084DAB9 && \
	add-apt-repository 'deb https://cloud.r-project.org/bin/linux/ubuntu focal-cran40/' && \
	apt -y -qq install r-base

# Install R packages
RUN R -e "install.packages(c('dplyr','plyr','tidyr','data.table','eulerr','circlize','LDlinkR','reshape','ggplot2','RSQLite','argparser'), dependencies=T, repos='http://cran.us.r-project.org/')" && \
	R -e "if (!requireNamespace('BiocManager',quietly=T)) install.packages('BiocManager')" && \
	R -e "Sys.setenv(R_INSTALL_STAGED = FALSE)" && \
	R -e "install.packages('XML', repos = 'http://www.omegahat.net/R')" && \
	R -e "BiocManager::install(version = '3.12')" && \
	R -e "BiocManager::install('biomaRt')" && \
	R -e "BiocManager::install('limma')" && \
	R -e "BiocManager::install('regioneR')" && \
	R -e "BiocManager::install('ComplexHeatmap')" && \
	R -e "BiocManager::install('fgsea')" && \
	R -e "BiocManager::install('hypeR')"

# Add Version number
ADD VERSION .

# Install PostGWAS-tools
COPY . .
WORKDIR /
