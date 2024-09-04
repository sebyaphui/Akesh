"""
Handles all data entering the simulation (command line arguments, taking data from the config file etc.). Doesn't encompass/ interact with conf.py though
"""

from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter
import warnings
import pandas as pd
from types import SimpleNamespace
import os
import json
import numpy as np
import importlib
from dataclasses import dataclass


@dataclass
class ConfigType:
    simtype: str
    subtype: str


def get_file_dir():
    return os.path.join(os.path.dirname(__file__))


def get_config_type(simtype: str, config: pd.DataFrame) -> ConfigType:
    subtype = config["subtype"]

    if subtype == "":
        # Default
        subtype = "lattice"

    print(subtype)

    config_type = ConfigType(simtype, subtype)
    return config_type


def verify_and_extract_config(config, config_type: ConfigType):
    """
    Extract all the `required` and `optional` arguments and put _only_ these into the namespace. \n
    Raise errors if `required` arguments are not found. \n
    Present warnings (can these be suppressed?) if `ignored` arguments are provided. \n
    With these aspects, you guarantee you know exactly what parameters are needed, and the program only handles what it needs.
    """

    required, optional, ignored = get_argument_priorities_from_config_type(config_type)

    for arg in required:
        if config[arg] == "":
            raise ValueError(f"Required argument {arg} is not specified.")

    for arg in ignored:
        if config[arg] != "":
            print(
                f"Ignored argument {arg} has been specified, but will not be used in the simulation."
            )
        del config[arg]

    print(config)

    return config


def get_argument_priorities_from_config_type(config_type: ConfigType):
    """
    Get the priority of each argument, based on the configuration type. We have the following priorities:

    `0`: `required`
    `1`: `optional`
    `2`: `ignored`

    """
    if config_type.subtype == "":
        config_type.subtype = (
            "lattice"  # By default, we assume this is the desired subtype
        )

    required = []
    optional = []
    ignored = []

    args_df = pd.read_csv(get_file_dir() + "/config/arguments.csv", dtype="object")
    args = args_df[args_df["simtype"] == config_type.simtype][
        args_df["subtype"] == config_type.subtype
    ]  # We assume uniqueness in this file

    if args.empty:
        raise ValueError(
            f"There is no configuration with simtype {config_type.simtype} and subtype {config_type.subtype}"
        )

    args_series = next(args.iterrows())[
        1
    ]  # Extract the first row, pulling the second element out of the (index, row) object
    del args_series["simtype"]
    del args_series["subtype"]

    for key, value in args_series.items():
        if value == "0":
            required.append(key)
        elif value == "1":
            optional.append(key)
        elif value == "2":
            ignored.append(key)
        else:
            raise ValueError(
                f"Argument priorities not configured correctly for simtype {config_type.simtype} and subtype {config_type.subtype}. {key} had unexpected value {value}, rather than 0, 1, or 2."
            )

    return required, optional, ignored


def get_sim_configuration(simtype, config_name=None):
    df = pd.read_csv(
        get_file_dir() + "/config/configurations.csv",
        usecols=lambda x: x != "description",
        dtype="string",
        keep_default_na=False,
    )
    df = df.set_index("name")
    print(df.index.tolist())
    if config_name is None:
        config_name = input("Choose the config: ")

    if config_name == "":
        config_name = "Default"

    full_config = get_config_df_from_row(df, config_name)
    config_type = get_config_type(simtype, full_config)
    full_config["subtype"] = (
        config_type.subtype
    )  # In case we have to overwrite "" with our default substitute
    config = verify_and_extract_config(full_config, config_type)
    return get_config_namespace_from_df(config)


def get_config_namespace_from_df(config: pd.DataFrame):
    return SimpleNamespace(**config)


def try_to_numeric(value):
    if value == "":
        return value
    else:
        try:
            return pd.to_numeric(value)
        except ValueError:
            return value


def get_config_df_from_row(df: pd.DataFrame, config_name):
    """
    Extract the row from `configurations.csv` associated with your specified config.
    """
    try:
        config = df.loc[config_name]
    except KeyError:
        raise KeyError("Simulation name invalid.")
    config = config.apply(try_to_numeric)
    config["name"] = config.name
    return config


def read_phenotypes(filename):
    filepath = get_file_dir() + "/" + filename
    with open(filepath, "r") as file:
        phenotypes = [line.rstrip("\n") for line in file]
    return phenotypes


def get_matrix_function_from_config(matrix_config_path):

    with open(get_file_dir() + "/" + matrix_config_path) as file:
        config = json.load(file)

    if "from" not in config.keys():
        raise KeyError("Invalid configuration file: missing 'from' attribute")

    if config["from"] == "path":
        if "delimiter" not in config.keys() or config["delimiter"] == "":
            config["delimiter"] = " "

        matrix = np.loadtxt(config["path"], delimiter=config["delimiter"]) # Assuming you have a CTL as each row, and a tumour as each column, you need to transpose.

        def get_matrix(sim):
            return matrix

        # return get_matrix
    elif config["from"] == "function":
        if "where" not in config.keys() or config["where"] == "":
            try:
                get_matrix = globals()[config["function"]]
            except KeyError:
                raise KeyError(
                    "The specified function does not exist in the global scope here."
                )
        else:
            try:
                module = importlib.import_module(
                    config["where"]
                )  # Supports relative imports
            except ModuleNotFoundError:
                raise ModuleNotFoundError(
                    "The specified module containing the matrix function does not exist."
                )
            get_matrix = getattr(module, config["function"])

    else:
        raise NotImplementedError("Invalid 'from' value.")

    return get_matrix


def set_up_and_get_arguments():
    warnings.simplefilter("ignore")

    parser = ArgumentParser(formatter_class=ArgumentDefaultsHelpFormatter)
    parser.add_argument(
        "-sf", "--save_figure", help="Whether to save a graph of the simulation."
    )
    parser.add_argument(
        "-ow",
        "--overwrite",
        help="Whether to overwrite a pre-existing simulation if it exists.",
    )
    parser.add_argument(
        "-c",
        "--config",
        help="Name of the configuration in configurations.csv to load.",
    )
    args = vars(parser.parse_args())

    save_fig = args["save_figure"]
    overwrite = args["overwrite"]
    config_name = args["config"]

    return save_fig, overwrite, config_name


def reset_interrupt():
    path = get_file_dir() + "/config/conf.py"
    with open(path, "r") as file:
        text = file.read()

    text = text.replace("interrupt = True", "interrupt = False")

    with open(path, "w") as file:
        file.write(text)
