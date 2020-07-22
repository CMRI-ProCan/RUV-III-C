# syntax=docker/dockerfile:experimental
# Copyright (c) Jupyter Development Team.
# Copyright (c) CMRI-ProCan Team.
# Distributed under the terms of the Modified BSD License.

FROM cmriprocan/ruviii_c-base:0.0.2

USER root

ARG bitbucket_user

ENV DEBIAN_FRONTEND noninteractive

RUN --mount=type=secret,id=bitbucket_pw cd / \
 && git clone https://${bitbucket_user}:$(cat /run/secrets/bitbucket_pw)@bitbucket.org/cmriprocan/ruviii_c.git  \
 && cd ruviii_c \
 && git remote set-url origin https://bitbucket.org/cmriprocan/ruviii_c.git \
 && R CMD INSTALL . \
 && cd .. \
 && rm -rf ./ruviii_c

# Install Tini
RUN wget --quiet https://github.com/krallin/tini/releases/download/v0.10.0/tini \
 && echo "1361527f39190a7338a0b434bd8c88ff7233ce7b9a4876f3315c22fce7eca1b0 *tini" | sha256sum -c - \
 && mv tini /usr/local/bin/tini \
 && chmod +x /usr/local/bin/tini

# Configure container startup
ENTRYPOINT ["tini", "--"]

# Overwrite this with 'CMD []' in a dependent Dockerfile
CMD ["/bin/bash"]
