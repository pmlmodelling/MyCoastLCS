""" Module to compute residence time. It computes the time that a particle
spent on a cell. """

import numpy as np
import xarray as xr
from .GridBased import GridBased


class ResidenceTime(GridBased):

    def __init__(self, nbins, bins_option):
        """
        Residence time initializer.

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
        GridBased.__init__(self, nbins, bins_option, True)
        self.name = 'residence_time'
        self.abbrev = 'RESD'

    def get_counts(self, ds: xr.Dataset) -> np.array:
        """Computes the average residence time. For each particle, it
        aproximates the time that a particles spents on a cell.

        Args:

            ds (xr.Dataset): Description

        Returns:
            residence_time(np.array): array with average time spent at cell.

        """
        print('-> RESD  >> Computing... ')

        n_tzyx = list(map(np.size, self.centers))
        time_in_cell = np.zeros((n_tzyx))
        mask_dif_id = np.zeros((n_tzyx))

        # Assuming that all particles have same dt
        dt = np.array((ds.time[1]-ds.time[0]).values, dtype='timedelta64[s]')
        dt = dt/np.timedelta64(1, 's')
        nz, ny, nx = ds.z0.size, ds.y0.size, ds.x0.size
        n_p = nz*nx*ny

        counter = 0

        # We move the check of the dimension to the outest part of the loop
        # to avoid excesive number of check
        if time_in_cell.ndim == 3:
            for i in range(0, nz):
                for j in range(0, ny):
                    for k in range(0, nx):
                        if nz > 1:
                            r = np.c_[ds.z[:, i, j, k].values,
                                      ds.y[:, i, j, k].values,
                                      ds.x[:, i, j, k].values]
                        else:
                            r = np.c_[ds.z[:, j, k].values,
                                      ds.y[:, j, k].values,
                                      ds.x[:, j, k].values]
                        # If particle is NaN go to next
                        if np.all(r == np.nan):
                            continue
                        counts, _ = np.histogramdd(r, bins=self.bins)
                        time_in_cell += counts*dt  # Counts the time
                        mask_dif_id += (counts > 0)
                        counter += 1
                    print('-> RESD  >> ', '%3.2f' % (counter/n_p*100.),
                          '%', end="\r", flush=True)

        elif time_in_cell.ndim == 2:
            for i in range(0, nz):
                for j in range(0, ny):
                    for k in range(0, nx):
                        r = np.c_[ds.y[:, j, k].values,
                                  ds.x[:, j, k].values]
                        if np.all(r == np.nan):  # If nan particle, go to next
                            continue
                        counts, _ = np.histogramdd(r, bins=self.bins)
                        time_in_cell += counts*dt
                        mask_dif_id += (counts > 0)
                        counter += 1
                    print('-> RESD  >> ', '%3.2f' % (counter/n_p*100.),
                          '%', end="\r", flush=True)

        residence_time = time_in_cell/mask_dif_id  # average residence time
        return residence_time

    def get_residence_time(self, ds_input: xr.Dataset) -> xr.Dataset:
        """
        Get the residence time and append the result to the input dataset.

        Args:
            ds_input (xr.Dataset): Input dataset with particle positions.

        Returns:
            ds_output (xr.Dataset): Output dataset with concentrations.

        """
        self.init_grid(ds_input)
        ds_output = self.init_dataset(ds_input)
        resd = self.get_counts(ds_input)
        self.to_dataset(ds_output, resd)
        return ds_output
