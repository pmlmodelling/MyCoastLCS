# -*- coding: utf-8 -*-
""" Module to compute concentrations. It gets the raw counts in each cell
obtained using diferent bining domain options.
"""

import numpy as np
import xarray as xr
from .GridBased import GridBased


class Concentrations(GridBased):

    def __init__(self, nbins: int or list, bins_option: str):
        """
        Concentration initializer.

        It inherits from GridBased.

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

        Returns:
            None.


        """
        GridBased.__init__(self, nbins, bins_option, False)
        self.name = 'concentrations'
        self.abbrev = 'CONC'

    def get_counts(self, ds: xr.Dataset) -> np.array:
        """
        Gets the raw number of counts per cell at each timestep.

        Args:
            ds (xr.Dataset): Input dataset with particle positions

        Returns:
            concentrations (np.array): Array with number of particles per cell.

        """

        print('-> CONC  >> Computing..')

        n_tzyx = [ds.time.size] + list(map(np.size, self.centers))
        concentrations = np.zeros((n_tzyx))

        for i in range(0, n_tzyx[0]):
            print('-> CONC  >> ', (i/n_tzyx[0])*100., '%', end="\r")
            if concentrations.ndim == 4:
                r = np.c_[ds.z[i].values.flatten(),
                          ds.y[i].values.flatten(),
                          ds.x[i].values.flatten()]
            elif concentrations.ndim == 3:
                r = np.c_[ds.y[i].values.flatten(),
                          ds.x[i].values.flatten()]
            concentrations[i], _ = np.histogramdd(r, bins=self.bins)
        return concentrations

    def get_concentrations(self, ds_input: xr.Dataset) -> xr.Dataset:
        """
        Get the concentration and append the result to the input dataset.

        Args:
            ds_input (xr.Dataset): Input dataset with particle positions.

        Returns:
            ds_output (xr.Dataset): Output dataset with concentrations.

        """
        self.init_grid(ds_input)
        ds_output = self.init_dataset(ds_input)
        concentrations = self.get_counts(ds_input)
        self.to_dataset(ds_output, concentrations)
        return ds_output
