# MESA

[MESA](https://mesastar.org/) is a versatile open-source 1D stellar evolution code. For the Connector 25 project, we [install](https://docs.mesastar.org/en/latest/installation.html) MESA version [24.08.1](https://zenodo.org/records/13353788) and MESA SDK version [24.7.1/24.10.1/25.3.1 (Linux (Intel/AMD) / Mac OS (Intel) / Mac OS (ARM))](http://user.astro.wisc.edu/~townsend/static.php?ref=mesasdk) for bit-for-bit reproducibility.

## MESA container

### Build

We also provide a Docker definitions file to allow the production of a 
containerized version of MESA using the versions described above. 
To build the container, [install Docker](https://docs.docker.com/engine/install/).
Then download the zip for MESA version [24.08.1](https://zenodo.org/records/13353788) 
and run the command `sudo docker build -f /path/to/mesa.dockerfile -t mesa-latest .` 
in the same directory as `mesa-24.08.1.zip`.

### Run with Docker

To run the built container with Docker to check it has built correctly, 
run `sudo docker run -it mesa-latest`. This will drop you into a bash prompt. 
MESA is located in `/mesa-24.08.1` and the SDK in `/opt/mesadsk`. 
Run `exit` to end the container. 

To use the container to run the MESA step of the Connector with Docker, you will
need to bind your copy of the Connector repository to the container. Run the
command `sudo docker run -v /path/to/connector25:/connector25-it  mesa-latest`. 
Now the path `/connector25` within the container prompt will have access to your 
local copy of the Connector repository.

### Run with Singularity/Apptainer

To run the built container on an HPCC you will likely need to use Singularity
or Apptainer (these are basically the same thing as of March 2025). 
Pull the container from DockerHub using the command 
`singularity pull docker://mesa-connector`.
Run the container using `singularity run mesa-connector-latest`. Singularity
automatically mounts your home directory, so you can navigate to your clone of 
the Connector repository within the container prompt.

## STELLA container

### Run with Docker

Stella is part of the MESA container. 

# STIR

# TARDIS

Standard installation of TARDIS can be accomplished using our locked environment
file, with instructions [here](https://tardis-sn.github.io/tardis/installation.html). 
The release of TARDIS used as part of the Connector is [XXXX](Link to github release).

## TARDIS container

### Build

We also provide a Docker definitions file to allow the production of a 
containerized version of TARDIS using the version described above. 
To build the container, [install Docker](https://docs.docker.com/engine/install/) and run the command `sudo docker build -f /path/to/tardis.dockerfile -t tardis-latest .`.

### Run with Docker

To run the built container with Docker to check it has built correctly, 
run `sudo docker run -it tardis-latest`. This will drop you into a bash prompt. 
TARDIS is located in `/tardis-release-XXXX`. Run `exit` to end the container.

To use the container to run the TARDIS step of the Connector with Docker, you will
need to bind your copy of the Connector repository to the container. Run the
command `sudo docker run -v /path/to/connector25:/connector25-it  mesa-latest`.


### Run with Singularity/Apptainer

To run the built container on an HPCC you will likely need to use Singularity
or Apptainer (these are basically the same thing as of March 2025). 
Pull the container from DockerHub using the command 
`singularity pull docker://tardis-connector`.
Run the container using `singularity run tardis-connector-latest`. Singularity
automatically mounts your home directory, so you can navigate to your clone of 
the Connector repository within the container prompt.