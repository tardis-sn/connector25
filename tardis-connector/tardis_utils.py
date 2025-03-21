import yaml

from tardis.util.base import is_valid_nuclide_or_elem


def write_tardis_csvy(tardis_sample_csvy_path, modify_csvy_headers, df_csv, new_csvy_path):
    """
    Purpose:
    ---------
    Write the TARDIS model csvy file for a specific time step.

    ----------
    Parameters:
        modify_csvy_headers: dict
            The dictionary that contains the to-be modified headers
        new_csvy_path: str
            The path to the new csvy file
    """
    # Read the sample csvy file
    with open(tardis_sample_csvy_path, "r") as file:
        csvy_lines = file.readlines()

    # Find the lines between "---" and datatype -- these are the headers
    start_index = csvy_lines.index("---\n")
    end_index = csvy_lines.index("datatype:\n")

    # Parse the lines as YAML
    yml_lines = csvy_lines[start_index:end_index]
    yml_data = yaml.safe_load("".join(yml_lines))

    # Modify the header dictionary as needed
    for key, value in modify_csvy_headers.items():
        yml_data[key] = value

    # add the datatype fields
    fields = get_fields_names(df_csv.columns.to_list())
    yml_data["datatype"] = {"fields": fields}

    # Convert the yml data back to lines
    yml_lines = yaml.dump(yml_data, sort_keys=False).splitlines()
    yml_lines = [line + "\n" for line in yml_lines]

    # Convert the csv data to lines
    fields_columns = [field["name"] for field in fields]
    csv_lines = (
        df_csv[fields_columns].to_csv(index=False, float_format="%.5e", sep=",").splitlines()
    )
    csv_lines = [line + "\n" for line in csv_lines]

    # Save the updated csvy data
    updated_csvy_lines = csvy_lines[: start_index + 1] + yml_lines + ["---\n"] + csv_lines
    with open(new_csvy_path, "w") as file:
        file.writelines(updated_csvy_lines)


def get_fields_names(column_names):
    """get the field names for csvy file based the column names
       also get rid of the isotopes that is not valid in tardis database

    Parameters
    ----------
    column_names : _type_
        _description_

    Returns
    -------
    _type_
        _description_
    """
    fields = [
        {
            "name": "velocity",
            "unit": "cm/s",
            "desc": "velocities of shell outer bounderies.",
        },
        {
            "name": "density",
            "unit": "g/cm^3",
            "desc": "density within shell with corresponding outer velocity.",
        },
    ]

    if ("velocity" not in column_names) or ("density" not in column_names):
        ValueError("velocity and density are required fields.")
    else:
        column_names.remove("velocity")
        column_names.remove("density")

    if "t_rad" in column_names:
        fields.append(
            {
                "name": "t_rad",
                "unit": "K",
                "desc": "radiative temperature within shell with corresponding outer velocity.",
            }
        )
        column_names.remove("t_rad")

    if "dilution_factor" in column_names:
        fields.append(
            {
                "name": "dilution_factor",
                "desc": "dilution factor within shell with corresponding outer velocity.",
            }
        )
        column_names.remove("dilution_factor")

    # assume the rest are fractional abundance
    for element in column_names:
        # check the column is a valid element
        if is_valid_nuclide_or_elem(element):
            fields.append(
                {
                    "name": element,
                    "desc": f"fractional {element} abundance",
                }
            )
        else:
            print(f"{element} is not valid nuiclide in tardis database.")

    return fields


def write_tardis_config(
    tardis_sample_config_path, modify_parameters, new_config_path, csvy_model_path=None
):
    """
    Purpose:
    ---------
    Write the TARDIS config file for a specific time step.

    ----------
    Parameters:
        modified_parameters: dict
            The dictionary that contains the to-be modified parameters
        new_config_path: str
            The path to the new config file
    """
    # load in the sample config yml
    with open(tardis_sample_config_path, "r") as file:
        config = yaml.safe_load(file)

    # Modify the config dictionary as needed
    for key1, params in modify_parameters.items():
        for key2, value in params.items():
            config[key1][key2] = value

    if csvy_model_path is not None:
        config["csvy_model"] = csvy_model_path

    # Save the modified config back to a new YAML file
    with open(new_config_path, "w") as file:
        yaml.safe_dump(config, file, sort_keys=False)
