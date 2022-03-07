# -*- coding: utf-8 -*-
"""
ArrayToGrid module. This module contains a class to transform array data
indexing into grid indexing if the initial grid of particles verify the
constrits mentioned. The module turns the dataset with dimension:

    [time,id/particles] ---grid_shape----> [time,z0,y0,x0]

Where [z0,y0,x0] are the dim (labels) to mark the initial position os the grid
and allow to identify the initial neighbours of the particles to compute
the deformation.

To perform this change the grid_shape variable with dimensions [nz,ny,nx]
must be provided.


"""

import xarray as xr


class ArrayToGrid:

    def __init__(self):

        self.vars_labels = ['z', 'y', 'x']
        self.coords_labels = ['time', 'z0', 'y0', 'x0']
        self.ds = []

    def array_to_grid(self, ds: xr.Dataset, grid_shape):
        """

        This functions turns the input dataset comming from lagrangian model
        input into a compatible structured dataset to compute FTLE.

        Args:
            - ds(xr.Dataset) : netcdf xarray dataset with dimensions [id,time]
            - grid_shape(list): grid of shape of points.


        Returns:
            - ds(xarray): Return a xarray-dataset with dimensions
            [time,z0,y0,x0]
            FTLE computation.

        """

        n_points = ds.x.isel(time=0).size
        print('-> INPUT >> Points array -> Structured grid')
        print('-> INPUT >>', n_points, '->', grid_shape)

        nvars = 3
        grid_data = [ds.z, ds.y, ds.x]
        ds['x0'] = ds.x.isel(time=0).values.reshape(grid_shape)[0, 0, :]
        ds['y0'] = ds.y.isel(time=0).values.reshape(grid_shape)[0, :, 0]
        ds['z0'] = ds.z.isel(time=0).values.reshape(grid_shape)[:, 0, 0]
        coords_data = [ds.time.data, ds.z0.data, ds.y0.data, ds.x0.data]

        # Setting the grid coordinates
        coords = dict(zip(self.coords_labels,
                          zip(self.coords_labels, coords_data)))

        # Transform the variables into a structure grid
        variables_in_grid_form = [var_to_reshape.values.reshape(
            [ds.time.size]+grid_shape) for var_to_reshape in grid_data]

        # Transform the variables into a dataset format.
        variables_in_ds_form = dict(zip(self.vars_labels, zip(
             nvars*[self.coords_labels], list(variables_in_grid_form))))

        ds_output = xr.Dataset(variables_in_ds_form, coords=coords)
        ds_output = land_data_to_nan(ds_output)
        ds_output = squeeze_z_dim(ds_output)

        return ds_output


def squeeze_z_dim(ds: xr.Dataset):
    """ Turns the Lagrangian input dataset into a compatible structured
        dataset to compute FTLE.

        Args:
            - ds(xr.Dataset) : netcdf xarray dataset with dimensions [id,time]
            - grid_shape(list): grid of shape of points.
    """

    if ds['z0'].size == 1:
        ds = ds.squeeze()
        print('-> INPUT >> Squeezing z-dim')
    return ds


def land_data_to_nan(ds: xr.Dataset):
    """ Mask those values with no movement in all the timesteps from initial
    time instant.

    Args:
        - ds(xr.Dataset): netcdf xarray dataset with dimensions [id,time]
        - grid_shape(list): grid of shape of points.

    Returns:
        - ds(xr.Dataset): netcdf xarray dataset with dimensions [id,time]

    """
    mask = (ds.isel(time=0) == ds).all(dim='time')
    return ds.where(mask == False)
