# -*- coding: utf-8 -*-
"""
Aliasing module. Module to controls the correspondece between
the Lagrangian model output variables and the x,y,z of the MYCOASTFTLE module.
Each Lagrangian model produces the outputs with their internal variable names.
The alias dictionary "links" those variables with the
internal x,y,z of the MYCOAST_FTLE module.

 - model_input (string): Name of the Lagrangian model used.

    Inputs supported in netcdf format:
        - Pylag
        - LAGAR

Example:
    - pylag_output.nc
        output_variables: xpos, ypos, zpos
    -alias_dict:
        alis_dictionary={'xpos': 'x', 'ypos': 'y', 'zpos': 'z',
                         'time': 'time'}

In order to work properly, the lagrangian model should produced netcdf outputs
with dimensions [particle] or [id] and [time].

The lagrangian simulation should be started in grid format filling the whole
domain of the velocity field.

o---------o---------o---------o---------o
|         |         |         |         |
|         |         |         |         |
|         |         |         |         |
o---------o---------o---------o---------o
|         |         |         |         |
|         |         |         |         |
|         |         |         |         |
o---------o---------o---------o---------o
|         |         |         |         |
|         |         |         |         |
|         |         |         |         |
o---------o---------o---------o---------o
|         |         |         |         |
|         |         |         |         |
|         |         |         |         |
o---------o---------o---------o---------o


"""

import xarray as xr

class Aliasing:
    def get_alias(model_input=None, alias=None):
        """
        Get the alias dictionary (from the model input or from the json file)
        to rename the dataset variables names of your lagrangian model with the
        MYCOAST internal variables (x,y,z) to compute the FTLE.

        Args:
            - alias (dict): Dictionary to link the variable.
        """

        if model_input == 'pylag':
            alias = {'x': 'x', 'y': 'y', 'z': 'z',
                     'time': 'time', 'particles': 'particle_id'}
        elif model_input == 'lagar':
            alias = {'lon': 'x', 'lat': 'y', 'depth': 'z',
                     'time': 'time', 'id': 'particle_id'}
        elif isinstance(alias, dict):
            return alias
        else:
            raise ValueError("""provide an alias dictionary or model to link
              your variables with my [x, y, z, time, particle_id] variables""")

        return alias


    def rename_dataset(alias: dict, input_file: str) -> xr.Dataset:
        """
        Rename the variables/dimensions in the dataset using alias.

        Args:
            alias (dict): Dictionary with variable/dimension names
            input_file (str): Path to netCDF dataset.

        Returns:
            renamed_dataset (TYPE): Dataset with variables/dimension renamed.

        """

        raw_dataset = xr.open_dataset(input_file)
        renamed_dataset = raw_dataset.rename(alias)

        return renamed_dataset
