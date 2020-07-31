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
		bedtools

# Install R packages
RUN R -e "install.packages(c('dplyr','data.table','eulerr'), dependencies=T, repos='http://cran.us.r-project.org/')" && \
	R -e "if (!requireNamespace('BiocManager',quietly=T)) install.packages('BiocManager')" && \
    R -e "BiocManager::install('biomaRt')" && \
    R -e "BiocManager::install('limma')"