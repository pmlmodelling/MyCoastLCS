# -*- coding: utf-8 -*-

from .Aliasing import Aliasing
from .ArrayToGrid import ArrayToGrid
from .FTLE import FTLE
from .LCS import LCS
import json
import glob
import xarray as xr
import os

class Common:
    
    def __init__(self):
        self.json_file = []
        self.alias = []
        self.model = []
        self.grid_shape = []
        self.integration_time_index = []
        self.spherical_flag = []
        self.ftle_thrsh = []
        self.nr_neighb = []
        self.ridge_points_flag = []
        self.ftle_LCS_only = True
        self.disk_or_mem = 'disk'
        
    def read_json(self,case_json):
               
        json_file = open(case_json)
        self.setup_file = json.load(json_file)
        json_file.close()
        
        print('-> SETUP >>', json.dumps(self.setup_file,indent=4))
        
        self.alias = self.setup_file['common']['alias']
        self.model = self.setup_file['common']['model']
        self.grid_shape = self.setup_file['common']['grid_shape']
        
        if 'FTLE' in self.setup_file:
            self.integration_time_index = self.setup_file['FTLE']['integration_time_index']
            self.spherical_flag = self.setup_file['FTLE']['spherical_flag']
            
        if 'LCS' in self.setup_file:
            self.eval_thrsh = self.setup_file['LCS']['eval_thrsh']
            self.ftle_thrsh = self.setup_file['LCS']['ftle_thrsh']
            self.area_thrsh = self.setup_file['LCS']['area_thrsh']
            self.nr_neighb = self.setup_file['LCS']['nr_neighb']
            self.ridge_points_flag = self.setup_file['LCS']['ridge_points_flag']
    
    
    def save_ftle_lcs_data(self,ds):
        var_names = ['LCS_forward','FTLE_forward','LCS_backward','FTLE_backward']
        var_inside_ds = [var for var in var_names if var in ds]
        t0 = ds.time.isel(time=0)
        ds = ds[var_inside_ds]
        ds = ds.assign_coords(time=t0)
        ds = ds.expand_dims('time')
        return ds
     
    def process_ftle_one_file(self,input_file, output_file):
        inputdata = Aliasing(self.model)
        inputdata.get_alias(self.alias)
        ds = inputdata.get_input(input_file)
        
        array_ds = ArrayToGrid()
        grid_ds = array_ds.array_to_grid(ds,self.grid_shape)
        
        if 'FTLE' in self.setup_file:
            FTLE_extractor = FTLE(self.spherical_flag, self.integration_time_index)
            FTLE_extractor.get_ftle(grid_ds)
        
        if 'LCS' in self.setup_file:
            LCS_extractor = LCS(self.eval_thrsh, self.ftle_thrsh, self.area_thrsh, self.nr_neighb, self.ridge_points_flag)
            LCS_extractor.get_lcs(grid_ds)
        
        if self.ftle_LCS_only == True:
            grid_ds = self.save_ftle_lcs_data(grid_ds)
    
        grid_ds.to_netcdf(output_file)
        return grid_ds
    
    
    def run_ftle_lcs(self,input_file_path_pattern, output_file):
        nc_file_list = sorted(glob.glob(input_file_path_pattern))
        if len(nc_file_list) == 0.:
             print('-> There is not file to process. Exiting')
             return
        
        elif len(nc_file_list) == 1:
             print('-> INPUT >> Processing file: ',nc_file_list[0])
             self.process_ftle_one_file(nc_file_list[0], output_file)
        
        else:
            ds_step_list = []
            step_file_list = []
            step = 0
            for nc_ftle_field in nc_file_list:
                print('\n')
                print('-> INPUT >> Processing file:', step+1, 'of ',len(nc_file_list),'>>', nc_ftle_field)
                ds_step = self.process_ftle_one_file(nc_ftle_field, str(step).zfill(3)+'.nc')
                step_file_list.append(str(step).zfill(3)+'.nc')
                ds_step_list.append(ds_step)
                step = step + 1
            print('\n')
            print('-> OUT  >> Merging all steps into:',output_file)
            xr.concat(ds_step_list,dim ='time').to_netcdf(output_file)
            if os.path.exists(output_file):
                for step_file in step_file_list:
                    os.remove(step_file)
        return

