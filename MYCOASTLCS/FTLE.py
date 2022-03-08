# -*- coding: utf-8 -*-
"""
Module to compute the Finite Time Lyapunov Exponents using the Cauchy Green
Finite Time deformation tensor. For more information about the method [a]_:


FTLE-definition
===============
Finite-Time Lyapunov Exponents have been used exten sively to quantify mixing
and especially to extract persisting transport patterns in the flow,
referred to as LCS. The FTLE at a given location measures the maximum
stretching rate of an infinitesimal fluid parcel over the interval
:math:`[t_0 , t_0 + T]` starting at the point :math:`r_0` at time :math:`t_0`.


Numerical Computation
=====================
The standard method to compute the FTLE field starts with the initialization
of a grid of material points equally spaced in a structured grid at 
:math:`\mathbf{r_0}_{ij}`
at time :math:`t_0`.

This simplifies gradient computation and does not require a priory knowledge
of the flow topology, which makes the discussion most broadly applicable.
The initial locations of these points, as opposed to the final
locations, represent the locations at which FTLE will be computed – the FTLE
grid. Then, we using the particles trajectories given by the model solution
for each material point for a finite time interval, :math:`[t_0 ,t_0 + T]`,
we compute the numerical gradient of the flow map.

We use second order accurate central differences in the interior points and
first order (forward or backwards) differences at the boundaries.

"""
import xarray as xr
import numpy as np


class FTLE:
    """
    FTLE
    ====

    It provides the methods to compute the FTLE from a structured gridded
    dataset for 2D or 3D data.

    """

    def __init__(self, spherical_flag: bool, integration_time_index: int):
        """Init the object with the setup provided in order to extract ftle.

        Args:
            - spherical (bool): True/false. cartesian(false) or spherical(true)
            - integration_time(int): Integer number with the time index of the
          the last position of the particles.
        """

        self.spherical_flag = bool(spherical_flag)
        self.integration_time = integration_time_index

    def get_integration_time(self, ds: xr.Dataset) -> [float, int]:
        """ Gets the integration in seconds from the integration time index.

        Args:
            - ds (xr.Dataset): dataset with lagrangian simulation.

        Returns:
            - T(float): integration time in seconds
            - timeindex(int): integration time index in time dataset axis

        """
        if isinstance(self.integration_time, int):
            timeindex = self.integration_time
            T = (ds.time[timeindex]-ds.time[0]).values.astype('timedelta64[s]').astype('f8')
        if T == 0.0:
            T = (ds.time[timeindex]-ds.time[0]).values.astype('f8')*1e-9
        elif isinstance(self.integration_time, float):
            timeindex = -1
            T = self.integration_time

        return T, timeindex

    def get_ftle_2d_cartesian(self, ds: xr.Dataset) -> np.array:
        """ Gets the 2D FTLE field in cartesian coordinates.

            Shadden, Shawn & Lekien, Francois & Marsden, Jerrold. (2005).
            Definition and properties of Lagrangian coherent structures from
            finit-time Lyapunov exponents in two-dimensional aperiodic flows.
            Physica D. 212. 271-304. 10.1016/j.physd.2005.10.007.

        Args:
            - ds (xr.Dataset): dataset with lagrangian simulation.

        Returns:
            - ftle (np.array) : 2d dimensional array.

        """
        T, timeindex = self.get_integration_time(ds)
        x_T = ds.x.isel(time=timeindex).values.squeeze()
        y_T = ds.y.isel(time=timeindex).values.squeeze()
        dxdy, dxdx = np.gradient(x_T, ds.y0.values, ds.x0.values)
        dydy, dydx = np.gradient(y_T, ds.y0.values, ds.x0.values)
        ny, nx = dxdy.shape
        ftle = np.zeros([ny, nx])
        J = np.zeros([2, 2])
        for i in range(0, ny):
            for j in range(0, nx):
                J = np.array([[dxdx[i, j], dxdy[i, j]],
                              [dydx[i, j], dydy[i, j]]])
                C = np.dot(np.transpose(J), J)
                eig_lya, _ = np.linalg.eigh(C)
                ftle[i, j] = (1./np.abs(T))*np.log(np.sqrt(eig_lya.max()))
        return ftle

    def get_ftle_2d_spherical(self, ds: xr.Dataset) -> np.array:
        """ Gets the 2D FTLE field in in lat - lon coordinates.

        Args:
            - ds (xr.Dataset): grid structured dataset with lagrangian
            simulation.

        Returns:
            - ftle (np.array) : 2d dimensional array.

        """
        T, timeindex = self.get_integration_time(ds)
        R = 6370000.
        x_T = ds.x.isel(time=timeindex).values.squeeze()
        y_T = ds.y.isel(time=timeindex).values.squeeze()
        dxdy, dxdx = np.gradient(x_T, ds.y0, ds.x0)
        dydy, dydx = np.gradient(y_T, ds.y0, ds.x0)
        ny, nx = np.shape(dxdx)
        ftle = np.zeros([ny, nx])
        theta = ds.y.isel(time=timeindex).squeeze().values
        for i in range(0, ny):
            for j in range(0, nx):
                J = np.array([[dxdx[i, j], dxdy[i, j]],
                              [dydx[i, j], dydy[i, j]]])
                M = np.array(
                    [[R*R*np.cos(theta[i, j]*np.pi/180.), 0], [0., R*R]])
                C = np.dot(np.dot(np.transpose(J), M), J)
                eig_lya, _ = np.linalg.eigh(C)
                ftle[i, j] = (1./np.abs(T))*np.log(np.sqrt(eig_lya.max()))
        return ftle

    def get_ftle_3d_cartesian(self, ds: xr.Dataset) -> np.array:
        """ Gets the 3D FTLE field in cartesian coordinates.

        Args:
            - ds (xr.Dataset): grid structured dataset with lagrangian
            simulation.

        Returns:
            - ftle (np.array) : 3d dimensional array.

        """
        T, timeindex = self.get_integration_time(ds)
        x_T = ds.x.isel(time=timeindex).values.squeeze()
        y_T = ds.y.isel(time=timeindex).values.squeeze()
        z_T = ds.z.isel(time=timeindex).values.squeeze()
        dxdz, dxdy, dxdx = np.gradient(x_T, ds.z0, ds.y0, ds.x0)
        dydz, dydy, dydx = np.gradient(y_T, ds.z0, ds.y0, ds.x0)
        dzdx, dzdy, dzdz = np.gradient(z_T, ds.z0, ds.y0, ds.x0)
        nz, ny, nx = dxdz.shape
        ftle = np.zeros([nz, ny, nx])
        J = np.zeros([3, 3])
        for i in range(0, nz):
            for j in range(0, ny):
                for k in range(0, nx):
                    J = np.array([[dxdx[i, j, k], dxdy[i, j, k], dxdz[i, j, k]],
                                  [dydx[i, j, k], dydy[i, j, k], dydz[i, j, k]],
                                  [dzdx[i, j, k], dzdy[i, j, k], dzdz[i, j, k]]])
                    C = np.dot(np.transpose(J), J)
                    eig_lya, _ = np.linalg.eig(C)
                    ftle[i, j, k] = (1./np.abs(T))*np.log(np.sqrt(eig_lya.max()))

        return ftle

    def to_dataset(self, ds: xr.Dataset, ftle: np.array) -> xr.Dataset:
        """ writes the ftle field in the dataset.

        Args:
            - ds (xr.Dataset): grid structured dataset with lagrangian
            simulation.
            - ftle (array) : 3d dimensional array.

        Returns:
            -  ds (xr.Dataset): grid structured dataset with ftle field

        """
        T, timeindex = self.get_integration_time(ds)
        if T > 0:
            ds['FTLE_forward'] = (ds.x.isel(time=timeindex).dims, ftle)
        if T < 0:
            ds['FTLE_backward'] = (ds.x.isel(time=timeindex).dims, ftle)
        return ds

    def get_ftle(self, ds: xr.Dataset, to_dataset=True):
        """

        It computes the FTLE (Finite Time Lyapunov Exponents) using the
        Cauchy Green finite time deformation tensor described, Shadden (2005).


        Args:
            - to_dataset (bool, optional): By default, it added the computed
            FTLE field, to the output dataset

        Returns:
            -  ds (xr.Dataset): grid structured dataset with ftle field

        """
        T, _ = self.get_integration_time(ds)
        t0 = ds.time.isel(time=0).values

        print('-> FTLE  >> Computing...')
        print('-> FTLE  >> t0 (starting time): ', t0)
        print('-> FTLE  >> T (advection time): ', T, 's')

        flag_3d = hasattr(ds, 'z0') & (ds.z0.size > 1)

        if (self.spherical_flag is False) and (flag_3d is False):
            ftle = self.get_ftle_2d_cartesian(ds)
        elif (self.spherical_flag is True) and (flag_3d is False):
            ftle = self.get_ftle_2d_spherical(ds)
        elif (self.spherical_flag is False) and (flag_3d is True):
            ftle = self.get_ftle_3d_cartesian(ds)
        elif((self.spherical_flag is True) and (flag_3d is True)):
            print('-> No spherical 3D FTLE available at the moment')
            print('-> Choose cartersian-2d or spherical_2d or 3d')
            return

        if to_dataset is True:
            self.to_dataset(ds, ftle)

        return ftle

    def explore_ftle_timescale(self, ds: xr.Dataset) -> xr.Dataset:
        """
        It computes the FTLE for all timesteps instead of a given one.
        The output produced will help you to explore the timescale of the
        deformation in order to infer the attributes for LCS and FTLE
        computation.


        Args:
            - ds (xr.Dataset): grid structured dataset with ftle field

        Returns:
            - ds (xr.Dataset) ds_output with FTLE computed for all timesteps.

        """
        print('-> FTLE >> Exploring time-scale T (advection time)')
        ftle = np.zeros_like(ds.x.values)
        nsteps = ds.time.size

        for i in range(0, nsteps):
            self.integration_time = i
            ftle[i] = self.get_ftle(ds, to_dataset=False)
            print('-> FTLE  >> ' + 'time index: ' + str(i) + '/' + str(nsteps) + '\r')
        ds['FTLE'] = (ds.x.dims, ftle)
        return ds
