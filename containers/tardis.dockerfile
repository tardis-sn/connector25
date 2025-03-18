FROM condaforge/miniforge3
RUN apt-get update && \
    DEBIAN_FRONTEND=noninteractive apt-get install -y \
      curl
RUN curl -L https://github.com/tardis-sn/tardis/archive/refs/tags/release-2025.03.16.tar.gz | \
    tar xz
RUN curl -L https://github.com/tardis-sn/tardis/releases/download/release-2025.03.16/conda-linux-64.lock
WORKDIR /tardis-release-2025.03.16
RUN conda create --file conda-linux-64.lock --name tardis
# handles conda init and conda activate because of bashrc changes
RUN conda init && . /root/.bashrc && conda activate tardis
RUN pip install .