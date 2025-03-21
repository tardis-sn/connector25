# Running the Connector on SLURM HPCC

Follow the instructions in the [installation guide](installation.md), checking
that the containers run using the Singularity or Apptainer system on your HPCC.
Make sure you login using X-Window forwarding (e.g., `ssh -X username@hpcc.example.com`).

## MESA

`salloc -n 16 -t 8:00:00 --x11` will request 16 CPUs for 8 hours including an
X11 binding so you can inspect the simulation as it runs. Once your job launches,
you can run the container. Once in the container, `cd` to 
`connector25/mesa-connector/template_binary_dev`. Run `./mk` then `./rn`. You
should see MESA begin to compute the evolution. The X-window based plots should
also begin to display after some time.

## STIR



## MESA + STELLA

`salloc -n 16 -t 8:00:00 --x11` will request 16 CPUs for 8 hours including an
X11 binding so you can inspect the simulation as it runs. Once your job launches,
you can run the container. Once in the container, `cd` to 
`connector25/mesa-connector/stella`. Run XXX then XXX. You
should see MESA begin to compute the evolution. The X-window based plots should
also begin to display after some time.

## TARDIS

`salloc -n 32 -t 1:00:00` will request 32 CPUs for 1 hour. Once your job launches,
you can run the container. Once in the container, `cd` to 
`connector25/tardis-connector/`. Run `conda activate tardis` then
`ipython tardis_run.py`. You should see a progress bar and text log output.
