# connector25
TARDIS Connector 25

tardis_run.py documentation --

tardis_run.py is a pared down script created in the 2025 tardis connector pipeline that runs tardis using a specificied configuration file, synthesizes a spectrum, and then dumps the output into an hdf. It is builting using tardis workflows, including the InnerVelocitySolver workflow which determines a photospheric inner boundary where the rosseland mean opacity is 2/3rds. This could be changed in the script if desired. 

Note that the tardis_config.yml specified must be set up appropriately for a tardis run. 

Inside the script, workflow.run() runs the entire tardis simulation. 

## Installation

This script is intended to be used in the tardis 2025 connector docker container, but will also work with an unspecific tardis installation of release 2025.03.19 https://github.com/tardis-sn/tardis/releases/tag/release-2025.03.19.


## Usage

python tardis_run.py([tardis_config.yml], [output.hdf])
