#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 10 12:58:52 2020

@author: gfnl143
"""

import xarray as xr

class ArrayToGrid:
    
    def __init__(self):
        
        self.x = []
        self.y = []
        self.z = []
        self.time = []
        
        self.x0 = []
        self.y0 = []
        self.z0 = []
        
        
        self.vars_labels = ['x', 'y','z']
        self.coords_labels = ['time', 'z0','y0', 'x0']
        self.ds = []
        
    
    def array_to_grid(self,ds,grid_shape):
        """
    
        This functions turns the input dataset comming from lagrangian model
        input into a compatible structured dataset to compute FTLE.
        New models will be available in the future.
    
    
        Args:
            - self(MYCOAST_FTLE) : MYCOAST_FTLE instance
    
        Returns:
            - self(MYCOAST_FTLE): Return a MYCOAST_FTLE instance ready for FTLE computation.
    
        """
    
        n_points = ds.x.isel(time=0).size 
        print('-> INPUT >> Points array -> Structured grid')
        print('-> INPUT >>', n_points,'->',grid_shape)
        
        nvars = 3
        grid_data = [ds.x, ds.y, ds.z]
        ds['x0'] = ds.x.isel(time=0).values.reshape(grid_shape)[0, 0, :]
        ds['y0'] = ds.y.isel(time=0).values.reshape(grid_shape)[0, :, 0]
        ds['z0'] = ds.z.isel(time=0).values.reshape(grid_shape)[:, 0, 0]
        coords_data = [ds.time, ds.z0, ds.y0, ds.x0]
            
        # Setting the grid coordinates
        coords = dict(zip(self.coords_labels, zip(self.coords_labels, coords_data)))
          
         # Transform the variables into a structure grid
        variables_in_grid_form = [var_to_reshape.values.reshape(
            [ds.time.size]+grid_shape) for var_to_reshape in grid_data]
         
        #Transform the variables into a dataset format.
        variables_in_ds_form = dict(zip(self.vars_labels, zip(
             nvars*[self.coords_labels], list(variables_in_grid_form))))
        
        ds_output = xr.Dataset(variables_in_ds_form, coords=coords)
        
        if grid_shape[0] == 1:
            ds_output = ds_output.squeeze()
            print('-> INPUT >> Squeezing z-dim')

        return ds_output
          