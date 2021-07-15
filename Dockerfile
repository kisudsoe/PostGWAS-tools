FROM ubuntu:latest
WORKDIR /

MAINTAINER Suhlab-Seungsoo_Kim

RUN apt -y -qq update --fix-missing && \
	apt -y -qq upgrade

# Set locale & R v4.0 repository
RUN DEBIAN_FRONTEND=noninteractive apt -y -qq install tzdata

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
	sqlite3 \
	parallel

# Install R
RUN apt-key adv --keyserver keyserver.ubuntu.com --recv-keys E298A3A825C0D65DFD57CBB651716619E084DAB9 && \
	add-apt-repository 'deb https://cloud.r-project.org/bin/linux/ubuntu focal-cran40/' && \
	apt -y -qq install r-base \
	libcairo2-dev \
	libxt-dev

# Install R packages
RUN R -e "install.packages(c('dplyr','plyr','tidyr','tidyverse','data.table','eulerr','circlize','LDlinkR','reshape','ggplot2','RSQLite','argparser','future.apply','corrplot'), dependencies=T, repos='http://cran.us.r-project.org/')" && \
	R -e "if (!requireNamespace('BiocManager',quietly=T)) install.packages('BiocManager')" && \
	R -e "Sys.setenv(R_INSTALL_STAGED = FALSE)" && \
	R -e "install.packages('XML', repos = 'http://www.omegahat.net/R')" && \
	#RUN R -e "BiocManager::install(version = '3.12')" && \
	R -e "BiocManager::install(c('biomaRt','limma','regioneR','ComplexHeatmap','fgsea','hypeR','clusterProfiler'))"
#RUN R -e "library(devtools); install_bitbucket('ibi_group/disgenet2r')" 

# Add Version number
ADD VERSION .

# Install PostGWAS-tools
COPY . .
WORKDIR /
