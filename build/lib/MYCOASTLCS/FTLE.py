#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 10 11:00:45 2020

@author: gfnl143
"""

import xarray as xr
import numpy as np


class FTLE:
    """ 

    Attributes:
        - alias (dict): Dictionary to link the variable outputs of your lagrangian model with the standard x,y,z variables to compute the FTLE.
        - ds (xarray.dataset): Netcdf input file in xarray.dataset format.
        - ds_output (xarray.dataset): Netcdf input file in xarray.dataset format.
        - grid_coords_labels (list): Internal coords names used by MYCOASTFTLE module.
        - grid_shape (list): Dimensions of the grid used to perform the advection (optional).
        - model_input (string): Name of the Lagrangian model used. 
        - spherical (boolean): True/false if the model is in cartesian or spherical coordinates.
        - time (ds.Datarray):Datarray of 'time' aliases positions from netcdf input.
        - vars_labels (list): Private variables names used by MYCOASTFTLE module.
        - x (xarray.Datarray): Datarray of 'x' aliases positions from netcdf input.
        - x0 (xarray.Datarray): Datarray of 'x0' grid coords.  
        - y (xarray.Datarray): Datarray of 'y' aliases positions from netcdf input.
        - y0 (xarray.Datarray): Datarray of 'y0' grid coords.
        - z (xarray.Datarray): Datarray of 'z' aliases positions from netcdf input.
        - z0 (xarray.Datarray):Datarray of 'z0' grid coords.



    """

    def __init__(self, spherical, integration_time):
        """
        
        Init the object with the setup provided in order to extract ftle.
        
        Args:
            - alias (dict): Dictionary to link the variable names of your lagrangian model with the private variables with the standard x,y,z variables used to compute the FTLE.
            - spherical (boolen): True/false if the model is in cartesian or spherical coordinates.
            - model_input (string): Name of the Lagrangian model used. 
            - grid_shape (list): Dimensions of the grid of initial conditions used to perform the advection (optional)

            Inputs supported:
                - Pylag
                - LAGAR

        Example:
            - pylag_ouput.nc
                output_variables: xpos, ypos, zpos
            -alias_dict:
                alis_dictionary={'x':'xpos','y':'ypos','z':'zpos'}   

        Deleted Parameters:
            **output_dict: Description

        """

        # Here we stablish the correspondece between the Lagrangian model input and also
        # the private variables x,y,z of the MYCOASTFTLE module.
        # Each model produces the outputs with their internal variable names.

        self.spherical = spherical
        self.integration_time = integration_time


        # grid_shape variable is used for those Lagrangian models which provided the
        # the evolution of particle trajectories in a array of the form. x[time,particles_id],
        # y[time,particles_id]. Other models produce ouputs in the format x[time,x_0,y_0]
        # the var_labels and and grid_coords_labels, are the internal (or private) variables
        # used to compute the FTLE and LCS, in 2D or 3D dimensions.
        # time index on lagrangian used to compute the FTLE
    
    def get_integration_time(self,ds):
        if type(self.integration_time) == int:
            timeindex = self.integration_time
            T = (ds.time[timeindex]-ds.time[0]).values.astype('timedelta64[s]').astype('f8')
        if T == 0.0:
            T = (ds.time[timeindex]-ds.time[0]).values.astype('f8')*1e-9
        elif type(self.integration_time) == float:
            timeindex = -1
            T = self.integration_time      
        return T,timeindex
    
    
    def get_ftle_2d_cartesian(self,ds):
        T,timeindex = self.get_integration_time(ds)
        dxdy, dxdx = np.gradient(ds.x.isel(time=timeindex).squeeze(), ds.y0.values, ds.x0.values)
        dydy, dydx = np.gradient(ds.y.isel(time=timeindex).squeeze(), ds.y0.values, ds.x0.values)
        ny, nx = dxdy.shape
        ftle = np.zeros([ny, nx])
        J = np.zeros([2, 2])
        for i in range(0, ny):
            for j in range(0, nx):
                J = np.array([[dxdx[i, j], dxdy[i, j]],
                              [dydx[i, j], dydy[i, j]]])
                C = np.dot(np.transpose(J), J)
                eigLya, _ = np.linalg.eigh(C)
                ftle[i, j] = (1./T)*np.log(np.sqrt(eigLya.max()))
        return ftle
    
    
    def get_ftle_2d_spherical(self,ds):
        T,timeindex = self.get_integration_time(ds)
        R = 6370000.
        dxdy, dxdx = np.gradient(ds.x.isel(time=timeindex).values.squeeze(), ds.y0, ds.x0)
        dydy, dydx = np.gradient(ds.y.isel(time=timeindex).values.squeeze(), ds.y0, ds.x0)
        ny, nx = np.shape(dxdx)
        ftle  = np.zeros([ny, nx])
        theta = ds.y.isel(time=timeindex).squeeze().values
        for i in range(0, ny):
            for j in range(0, nx):
                J = np.array([[dxdx[i, j], dxdy[i, j]],
                              [dydx[i, j], dydy[i, j]]])
                M = np.array(
                    [[R*R*np.cos(theta[i, j]*np.pi/180.), 0], [0., R*R]])
                C = np.dot(np.dot(np.transpose(J), M), J)
                eigLya, _ = np.linalg.eigh(C)
                ftle[i, j] = (1./T)*np.log(np.sqrt(eigLya.max()))
        return ftle


    def get_ftle_3d_cartesian(self,ds):
       T,timeindex = self.get_integration_time(ds)
       dxdz, dxdy, dxdx = np.gradient(ds.x.isel(time=timeindex).values.squeeze(), ds.z0, ds.y0, ds.x0)
       dydz, dydy, dydx = np.gradient(ds.y.isel(time=timeindex).values.squeeze(), ds.z0, ds.y0, ds.x0)
       dzdx, dzdy, dzdz = np.gradient(ds.z.isel(time=timeindex).values.squeeze(), ds.z0, ds.y0, ds.x0)
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
                   eigLya, _ = np.linalg.eig(C)
                   ftle[i, j, k] = (1./T)*np.log(np.sqrt(eigLya.max()))
       
       return ftle
    
    
    def to_dataset(self,ds,ftle):
        T,timeindex = self.get_integration_time(ds)
        if T>0:
            ds['FTLE_forward'] = (ds.x.isel(time=timeindex).dims, ftle)
        if T<0:
            ds['FTLE_backward'] = (ds.x.isel(time=timeindex).dims, ftle)
        return ds
    
    
    def get_ftle(self,ds,to_dataset=True):
        """

        It computes the FTLE (Finite Time Lyapunov Exponents) using the 
        Cauchy Green finite time deformation tensor described, Shadden (2005).
        

        Args: 
            - T (float, optional): If the dataset has no time attribute in datetimeformat,
            you can provide the advection time.
            - timeindex (int, optional): Index time to select the to compute the FTLE.
            By default is the last time step of the Lagrangian advection.
            - to_dataset (bool, optional): By default, it added the computed FTLE field,
            to the output dataset

        Returns:
            self: Added the FTLE_forward or FTLE_backward field to the ds_output.

        """
        T,timeindex = self.get_integration_time(ds)
        t0 = ds.time.isel(time=0).values
        
        print('-> FTLE  >> Computing...')
        print('-> FTLE  >> t0 (starting time): ', t0)
        print('-> FTLE  >> T (advection time): ' ,T, 's')
        
        flag_3d = hasattr(ds, 'z0') & (ds.z0.size > 1)

        if (self.spherical == False) and (flag_3d == False):
            ftle = self.get_ftle_2d_cartesian(ds)
        elif (self.spherical == True) and (flag_3d == False):
            ftle = self.get_ftle_2d_spherical(ds)
        elif (self.spherical == False) and (flag_3d == True):
            ftle = self.get_ftle_3d_cartesian(ds)
        elif((self.spherical == True) and (flag_3d == True)):
            print('-> No spherical 3D FTLE available at the moment')
            print('-> Choose another option: cartersian-2d or spherical_2d or 3d')
            return 

        if to_dataset == True:
            self.to_dataset(ds,ftle)

        return ds
            

    def explore_ftle_timescale(self,ds):
        """

        It computes the FTLE for all timesteps instead of a given one.
        The output produced will help you  to explore  the timescale of the deformation 
        in order to infer the attributes for LCS and FTLE computation.


        Args:
            - to_dataset (bool, optional): By default, it added the computed FTLE field,
            to the output dataset.
 

        Returns:
            self(MYCOASTFTLE): ds_output with FTLE computed for all timesteps.

        """

        ftle = np.zeros_like(ds.x.values)
        nsteps = ds.time.size

        for i in range(0, nsteps):
            self.integration_time = i
            ftle[i] = self.get_ftle(ds,to_dataset=False)
            print('-> FTLE >> field for step'+str(i)+'Percentage:'+str(100.*i/nsteps)+'\n')
        ds['FTLE'] = (ds.x.dims, ftle)
        return 