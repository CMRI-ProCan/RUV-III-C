FROM r-base:3.6.1

RUN R -e "install.packages(c('RSpectra', 'progress', 'ruv', 'Rcpp', 'RcppEigen'))"

RUN R -e "install.packages('testthat')"

#Update docker to correct version

RUN apt update \
 && apt install -y gnupg curl git \
 && export VER=18.09.3 \
 && curl -L -o /tmp/docker-$VER.tgz https://download.docker.com/linux/static/stable/x86_64/docker-$VER.tgz \
 && tar -xz -C /tmp -f /tmp/docker-$VER.tgz \
 && mv /tmp/docker/* /usr/bin

