#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 10 11:07:09 2020

@author: gfnl143
"""

import xarray as xr

class Aliasing:

    def __init__(self,model_input):
        
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
        
        self.model_input = model_input
        self.alias = []
        self.ds_output = []
        self.ds = []

        # Here we stablish the correspondece between the Lagrangian model input and also
        # the private variables x,y,z of the MYCOASTFTLE module.
        # Each model produces the outputs with their internal variable names.
        # The alias dictionary "links" those variables with the
        # internal x,y,z of the MYCOAST_FTLE module.
        # There are Two defaults dictionaries for LAGAR and PYLAG.
    
        
    def get_alias(self,alias):
        if self.model_input == 'pylag':
            self.alias = {'xpos': 'x', 'ypos': 'y','zpos': 'z', 'time': 'time'}
        elif self.model_input == 'lagar':
            self.alias = {'lon': 'x', 'lat': 'y', 'depth': 'z', 'time': 'time'}           
        else:
            self.alias = alias
    
    def get_input(self, netcdf_file):
        
        ds = xr.open_dataset(netcdf_file).fillna(0)
        ds = ds.rename(self.alias)
        
        return ds
    
    
        

        
    

        
    

    
    
    
    