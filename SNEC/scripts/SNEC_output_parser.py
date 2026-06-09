import os
import numpy as np
import re
import sys


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
