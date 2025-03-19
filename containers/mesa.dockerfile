FROM ubuntu:focal
# prereqs for MESASDK
RUN apt-get update && \
    DEBIAN_FRONTEND=noninteractive apt-get install -y \
      binutils \
      curl \
      make \
      perl \
      libx11-dev \
      zlib1g \
      zlib1g-dev \
      tcsh \
      git \
      git-lfs \
      python3 \
      ruby \
      && \
    apt-get autoremove --purge -y && \
    apt-get autoclean -y && \
    rm -rf /var/cache/apt/* /var/lib/apt/lists/*
# MESASDK from Zenodo
RUN curl -L https://zenodo.org/records/13768913/files/mesasdk-x86_64-linux-24.7.1.tar.gz | \
    tar xzf - -C /opt/
ENV MESASDK_ROOT=/opt/mesasdk
# Need to set shell to bash for source to work as part of MESA install
SHELL ["/bin/bash", "-c"]
# Need to download mesa-24.08.1.zip from MESA website first 
# and store in the same directory as this Dockerfile
RUN --mount=type=bind,source=mesa-24.08.1.zip,target=mesa-24.08.1.zip \
unzip mesa-24.08.1.zip
RUN echo $MESASDK_ROOT
ENV MESA_DIR=/mesa-24.08.1
# Remove the lines that check sudo install and source MESASDK to install MESA
RUN cd $MESA_DIR && sed -i '82,103d' "install" && source $MESASDK_ROOT/bin/mesasdk_init.sh && ./install
CMD ["/bin/bash", "-c", "source /opt/mesasdk/bin/mesasdk_init.sh && exec /bin/bash"]

