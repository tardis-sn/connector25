FROM ubuntu:focal
RUN apt-get update && \
    DEBIAN_FRONTEND=noninteractive apt-get install -y \
    curl \
    make \
    gfortran \
    libblas-dev \
    liblapack-dev \
    && \
    apt-get autoremove --purge -y && \
    apt-get autoclean -y && \
    rm -rf /var/cache/apt/* /var/lib/apt/lists/*
RUN curl https://stellarcollapse.org/codes/SNEC-1.01-20161007.tar.gz | \
    tar xzf -
# comment out Mac OS X flags and uncomment Linux flags
RUN sed -i '7s/^#//; 9s/^/#/' /SNEC-1.01/make.inc
WORKDIR /SNEC-1.01
RUN make