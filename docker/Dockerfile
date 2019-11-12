# syntax=docker/dockerfile:experimental
# Copyright (c) Jupyter Development Team.
# Copyright (c) CMRI-ProCan Team.
# Distributed under the terms of the Modified BSD License.

FROM ubuntu:18.04

USER root

ENV DEBIAN_FRONTEND noninteractive

COPY docker/packages.txt .

RUN apt-get update \
 && apt-get -yq dist-upgrade \
 && cat packages.txt | xargs apt-get install -y --no-install-recommends \
 && apt-get clean \
 && rm -rf /var/lib/apt/lists/*

RUN rm packages.txt

COPY docker/Rpackages.txt .

RUN R -e "packages <- read.table('./Rpackages.txt', stringsAsFactors=FALSE); install.packages(packages[,1]);"

RUN rm Rpackages.txt

RUN R -e "source('https://bioconductor.org/biocLite.R'); biocLite('sva'); biocLite('preprocessCore')"

RUN --mount=type=secret,id=gitpassword cd / \
 && git clone https://rshah-cmri:$(cat /run/secrets/gitpassword)@bitbucket.org/cmriprocan/ruviii_c.git  \
 && cd ruviii_c \
 && git remote set-url origin https://bitbucket.org/cmriprocan/ruviii_c.git \
 && R CMD INSTALL . \
 && cd .. \
 && rm -rf ./ruviii_c

RUN cd / \
 && git clone https://github.com/rohan-shah/ruv.git \
 && R CMD INSTALL /ruv \
 && rm -r /ruv

# Install Tini
RUN wget --quiet https://github.com/krallin/tini/releases/download/v0.10.0/tini \
 && echo "1361527f39190a7338a0b434bd8c88ff7233ce7b9a4876f3315c22fce7eca1b0 *tini" | sha256sum -c - \
 && mv tini /usr/local/bin/tini \
 && chmod +x /usr/local/bin/tini

# Configure container startup
ENTRYPOINT ["tini", "--"]

# Overwrite this with 'CMD []' in a dependent Dockerfile
CMD ["/bin/bash"]
