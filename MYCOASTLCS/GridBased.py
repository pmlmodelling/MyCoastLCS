# -*- coding: utf-8 -*-
""" Module to create cell-grids. A cell-grid is a domain division into cells to
perform operations with particles on cell volumes/areas done by cell centers 
and cell edges."""

import numpy as np
import xarray as xr


class GridBased:

    def __init__(self, nbins: int or list, bins_option: str, static):
        """
        Grid constructor.

        Args:
            nbins (int or list): Integer/s with the number of bins per dim.
            bins_option (str): string with: "origin, domain,custom"

        Options:
            - "origin": It uses the initial grid position to define the domain
            limits to apply the "nbins" binning.
            - "domain": It uses the final particle positions to define the
            domain limits to apply the "nbins" binning.
            - "custom": It uses custom user domain limits to apply the "nbins"
            binning.

        """
        self.measure = ''
        self.abbrev = ''
        self.bins_option = bins_option
        self.nbins = nbins
        self.coords_labels = ['time', 'z_c', 'y_c', 'x_c']
        self.dims = ['time', 'z_c', 'y_c', 'x_c']
        self.static = static

    def get_bins(self, ds: xr.Dataset):
        """
        Gets the bins from an input xr.Dataset.

        It uses the "nbins" to create the bins with the "bins_option".

        Options:
            - "origin": It uses the initial grid position to define the domain
            limits to apply the "nbins" binning.
            - "domain": It uses the final particle positions to define the
            domain limits to apply the "nbins" binning.
            - "custom": It uses custom user domain limits to apply the "nbins"
            binning.

        Args:
            ds (xr.Dataset): DESCRIPTION.

        Returns:
            None.

        """

        if self.bins_option == 'origin':
            if self.nbins:
                x_bins = np.linspace(ds.x0.min(), ds.x0.max(), self.nbins)
                y_bins = np.linspace(ds.y0.min(), ds.y0.max(), self.nbins)
                z_bins = np.linspace(ds.z0.min(), ds.z0.max(), self.nbins)
            else:
                centers = ds.z0.values, ds.y0.values, ds.x0.values
                bins = []
                for center in centers:
                    if center.size == 1:
                        bins.append(np.ones(1))
                        print('->' + self.abbrev,
                              '  >> warning: just 1 dim in z depth. I Can''t estimate bins.')
                    else:
                        dr = [(center[1]-center[0])/2.]
                        bins.append(np.append(center-dr, center[-1]+dr))

            z_bins, y_bins, x_bins = bins

        elif self.bins_option == 'domain':
            if isinstance(self.nbins, int):
                x_bins = np.linspace(ds.x.min(), ds.x.max(), self.nbins)
                y_bins = np.linspace(ds.y.min(), ds.y.max(), self.nbins)
                z_bins = np.linspace(ds.z.min(), ds.z.max(), self.nbins)
            elif isinstance(self.nbins, list):
                x_bins = np.linspace(ds.x.min(), ds.x.max(), self.nbins[2])
                y_bins = np.linspace(ds.y.min(), ds.y.max(), self.nbins[1])
                z_bins = np.linspace(ds.z.min(), ds.z.max(), self.nbins[0])
        elif self.bins_option == 'custom':
            z_bins = np.linspace(*self.nbins[0])
            y_bins = np.linspace(*self.nbins[1])
            x_bins = np.linspace(*self.nbins[2])

        if (z_bins.size > 1) and (self.static is False):
            self.bins = (z_bins, y_bins, x_bins)
        elif (z_bins.size <= 1) and (self.static is False):
            self.bins = (y_bins, x_bins)
            self.coords_labels = ['time', 'y_c', 'x_c']
            self.dims = ['time', 'y_c', 'x_c']
        elif (z_bins.size > 1) and (self.static is True):
            self.bins = (z_bins, y_bins, x_bins)
            self.coords_labels = ['z_c', 'y_c', 'x_c']
            self.dims = ['z_c', 'y_c', 'x_c']
        elif (z_bins.size <= 1) and (self.static is True):
            self.bins = (y_bins, x_bins)
            self.coords_labels = ['y_c', 'x_c']
            self.dims = ['y_c', 'x_c']
        return

    @staticmethod
    def check_2d_data(ds: xr.Dataset) -> bool:
        """
        Check if a xr.Dataset is 2D (spatial).

        Args:
            ds (xr.Dataset): Input dataset

        Returns:
            bool: Flag. True(is 2D)

        """

        flag_2d = np.all(np.abs(ds.z.isel(time=-1) - ds.z.isel(time=0)) < 1e-6)
        return flag_2d

    def get_centers(self):
        """
        Get the center point based on bins defining each cell.

        Returns:
            None.

        """
        self.centers = [(self.bins[i][:-1] + self.bins[i][1:]) / 2.
                        for i in range(0, len(self.bins))]

    def print_bins_info(self):
        """Print bins information"""

        print('-> ' + self.abbrev,
              '>> x [min,max,nbins] -> [',
              self.bins[-1].min(), ',',
              self.bins[-1].max(), ',',
              self.bins[-1].size, ']')
        print('-> ' + self.abbrev,
              '>> y [min,max,nbins] -> [',
              self.bins[-2].min(), ',',
              self.bins[-2].max(), ',',
              self.bins[-2].size, ']')
        if len(self.bins) > 2:
            print('-> ' + self.abbrev,
                  '>> z [min,max,nbins] -> [',
                  self.bins[0].min(), ',',
                  self.bins[0].max(), ',',
                  self.bins[0].size, ']')

    def init_grid(self, ds_input: xr.Dataset):
        """
        Initialize the grid from a provided dataset.

        Args:
            ds_input (xr.Dataset): Input dataset.

        Returns:
            None.

        """
        self.get_bins(ds_input)
        self.get_centers()
        self.print_bins_info()

    def init_dataset(self, ds_input):
        if self.bins_option == 'origin':
            ds_output = xr.Dataset({})
            if 'time' in self.dims:
                ds_output = ds_output.assign_coords(ds_input.coords)
            else:
                # remove time dimension for static measures
                ds_output = ds_output.assign_coords(ds_input.coords)
                ds_output = ds_output.drop_dims('time')
        else:
            print('-> ' + self.abbrev, '>> Creating coords to store counting grid.')
            if 'time' in self.dims:
                coords = dict(zip(self.coords_labels,
                                  zip(self.coords_labels,
                                      [ds_input.time] + self.centers)))
            else:
                coords = dict(zip(self.coords_labels,
                                  zip(self.coords_labels, self.centers)))

            ds_output = xr.Dataset({}, coords=coords)
        return ds_output

    def to_dataset(self, ds_output: xr.Dataset, data: np.array):
        """
        Write numerical array into a dataset.

        Args:
            ds_output (xr.Dataset): Dataset to append numerical data.
            data (np.array): Numerical data.

        Returns:
            ds_output (xr.Dataset): Dataset with variable appended

        """
        ds_output[self.name] = (self.dims, data)
        return ds_output