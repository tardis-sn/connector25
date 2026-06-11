import os
import numpy as np
import re
import sys
import h5py


def SNEC_output_parser(outfile, cache=True):
    """
    Parse the file and return a dictionary with time values as keys
    and numpy arrays with two columns as values.

    Parameters:
    -----------
    filename : str
        Path to the input file

    Returns:
    --------
    data: np.array [time, mass, quantity]
    """
    # Check for cached binary version
    cache_file = outfile + ".npz"
    if os.path.exists(cache_file):
        cached = np.load(cache_file)
        return {float(k): cached[k] for k in cached.files}

    data_dict = {}

    with open(outfile, "r") as f:
        content = f.read()

    # Split by "Time =" to get each time block
    blocks = re.split(r'"Time =', content)

    for block in blocks[1:]:  # Skip the first empty split
        lines = block.strip().split('\n')

        # Extract time value from first line
        time_line = lines[0]
        time_value = float(time_line.split()[0])

        # Parse data lines
        data_lines = []
        for line in lines[1:]:
            line = line.strip()
            if line and (not line.startswith('"')) and (not line.startswith('#')):  # Skip empty lines and closing quotes
                parts = line.split()
                if len(parts) >= 2:
                    try:
                        col1 = float(parts[0].replace('E+0', 'E+'))
                        col2 = float(parts[1].replace('E+0', 'E+'))
                        data_lines.append([col1, col2])
                    except ValueError:
                        continue

        # Convert to numpy array
        if data_lines:
            data_dict[time_value] = np.array(data_lines)

    if cache:
        # Save to binary cache
        np.savez(cache_file, **{str(k): v for k, v in data_dict.items()})

    return data_dict


# For hdf5 output below
def read_SNEC_hdf5(filename):
    """
    Parse a SNEC HDF5 output file into a flat dictionary.
    Parameters
    ----------
    filename : str
        Path to output.h5

    Returns
    -------
    dict
        Keys are field/scalar names (e.g. 'fields/rho', 'scalars/lum_observed').
        Values are numpy arrays:
          - 1-D arrays of length imax for fields (modeflag=1 snapshots)
          - 1-D arrays of length n_timesteps for scalars (modeflag=2 time series)
        The root 'time' scalar (last written simulation time) is included as 'time'.
    """
    data = {}

    with h5py.File(filename, "r") as f:
        _recurse(f, "", data)

    return data


def _recurse(node, prefix, data):
    """Walk the HDF5 tree and collect all datasets into data."""
    for key in node.keys():
        item = node[key]
        path = f"{prefix}/{key}".lstrip("/")
        if isinstance(item, h5py.Dataset):
            arr = item[()]
            # squeeze length-1 dims (e.g. the root 'time' scalar stored as shape (1,))
            if arr.shape == (1,):
                arr = arr[0]
            data[path] = arr
        elif isinstance(item, h5py.Group):
            _recurse(item, path, data)



if __name__ == "__main__":
    filename = sys.argv[1] if len(sys.argv) > 1 else "output.h5"

    data = read_snec_hdf5(filename)

    print(f"Read {len(data)} entries from {filename}\n")
    print(f"{'Key':<40}  {'Shape':<20}  {'Dtype'}")
    print("-" * 70)
    for key, val in sorted(data.items()):
        val = np.atleast_1d(val)
        print(f"{key:<40}  {str(val.shape):<20}  {val.dtype}")
