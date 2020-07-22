# syntax=docker/dockerfile:experimental
# Copyright (c) Jupyter Development Team.
# Copyright (c) CMRI-ProCan Team.
# Distributed under the terms of the Modified BSD License.

FROM ubuntu:18.04

USER root

ENV DEBIAN_FRONTEND noninteractive

COPY packages.txt .

RUN apt-get update \
 && apt-get install -y gnupg \
 && echo deb https://cloud.r-project.org/bin/linux/ubuntu bionic-cran35/ >> /etc/apt/sources.list \
 && gpg --keyserver hkp://keyserver.ubuntu.com:80 --recv-keys E298A3A825C0D65DFD57CBB651716619E084DAB9 \
 && gpg -a --export E298A3A825C0D65DFD57CBB651716619E084DAB9 | apt-key add - \
 && apt-get install -y ca-certificates \
 && apt-get update && apt-get -yq dist-upgrade

RUN cat packages.txt | xargs apt-get install -y --no-install-recommends \
 && apt-get clean \
 && rm -rf /var/lib/apt/lists/*

RUN rm packages.txt

COPY Rpackages.txt .

RUN R -e "packages <- read.table('./Rpackages.txt', stringsAsFactors=FALSE); install.packages(packages[,1]);"

RUN rm Rpackages.txt

RUN R -e "install.packages('BiocManager'); BiocManager::install(c('sva', 'preprocessCore'))"

#Install MKL
RUN wget https://apt.repos.intel.com/intel-gpg-keys/GPG-PUB-KEY-INTEL-SW-PRODUCTS-2019.PUB \
 && apt-get update \
 && apt-get install -y gpg-agent \
 && apt-key add GPG-PUB-KEY-INTEL-SW-PRODUCTS-2019.PUB \
 && rm GPG-PUB-KEY-INTEL-SW-PRODUCTS-2019.PUB \
 && sh -c 'echo deb https://apt.repos.intel.com/mkl all main > /etc/apt/sources.list.d/intel-mkl.list' \
 && apt-get update \
 && apt-get install -y intel-mkl-64bit-2019.5-075 \
 && echo "/opt/intel/lib/intel64" > /etc/ld.so.conf.d/mkl.conf \
 && echo "/opt/intel/mkl/lib/intel64" >> /etc/ld.so.conf.d/mkl.conf \
 && ldconfig \
 && echo "MKL_THREADING_LAYER=GNU" >> /etc/environment

#Actualy use it
RUN update-alternatives --install /usr/lib/x86_64-linux-gnu/libblas.so libblas.so-x86_64-linux-gnu /opt/intel/mkl/lib/intel64/libmkl_rt.so 50 \
 && update-alternatives --install /usr/lib/x86_64-linux-gnu/libblas.so.3 libblas.so.3-x86_64-linux-gnu /opt/intel/mkl/lib/intel64/libmkl_rt.so 50 \
 && update-alternatives --install /usr/lib/x86_64-linux-gnu/liblapack.so liblapack.so-x86_64-linux-gnu /opt/intel/mkl/lib/intel64/libmkl_rt.so 50 \
 && update-alternatives --install /usr/lib/x86_64-linux-gnu/liblapack.so.3 liblapack.so.3-x86_64-linux-gnu /opt/intel/mkl/lib/intel64/libmkl_rt.so 50

# Install Tini
RUN wget --quiet https://github.com/krallin/tini/releases/download/v0.10.0/tini \
 && echo "1361527f39190a7338a0b434bd8c88ff7233ce7b9a4876f3315c22fce7eca1b0 *tini" | sha256sum -c - \
 && mv tini /usr/local/bin/tini \
 && chmod +x /usr/local/bin/tini

# Configure container startup
ENTRYPOINT ["tini", "--"]

# Overwrite this with 'CMD []' in a dependent Dockerfile
CMD ["/bin/bash"]
