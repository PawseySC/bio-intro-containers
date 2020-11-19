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
