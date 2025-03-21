# Two arguments are config model path and output file name
from tardis.workflows.v_inner_solver import InnerVelocitySolverWorkflow
from tardis.io.configuration.config_reader import Configuration
from pathlib import Path
import sys

config = Configuration.from_yaml(Path(sys.argv[1]))
workflow = InnerVelocitySolverWorkflow(
    config, tau=2.0 / 3, mean_optical_depth="rosseland", csvy=True
)
workflow.run()
spectrum = workflow.spectrum_solver.spectrum_real_packets
spectrum.to_hdf(sys.argv[2])
