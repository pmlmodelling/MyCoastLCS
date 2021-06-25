# -*- coding: utf-8 -*-

import json
import glob
import xarray as xr
import os
from . import Aliasing
from .ArrayToGrid import ArrayToGrid
from .FTLE import FTLE
from .LCS import LCS
from .Concentrations import Concentrations
from .ResidenceTime import ResidenceTime


class Common:

    def __init__(self):
        self.json_file = []
        self.alias = {}
        self.model = ''
        self.grid_shape = []
        self.ftle_LCS_only = True
        self.disk_or_mem = 'disk'

    def read_json(self, case_json):

        json_file = open(case_json)
        self.setup_file = json.load(json_file)
        json_file.close()

        print('-> SETUP >>', json.dumps(self.setup_file, indent=4))

        if 'model' in self.setup_file['common']:
            self.model = self.setup_file['common']['model']

        elif 'alias' in self.setup_file['common']:
            self.alias = self.setup_file['common']['alias']

        self.grid_shape = self.setup_file['common']['grid_shape']

    def save_ftle_lcs_data(self, ds):
        var_names = ['LCS_forward', 'FTLE_forward',
                     'LCS_backward', 'FTLE_backward']
        var_inside_ds = [var for var in var_names if var in ds]
        t0 = ds.time.isel(time=0)
        ds = ds[var_inside_ds]
        ds = ds.assign_coords(time=t0)
        ds = ds.expand_dims('time')
        return ds

    def process_one_file(self, input_file, output_file):
        alias = Aliasing.get_alias(self.model, self.alias)
        ds = Aliasing.rename_dataset(alias, input_file)

        array_ds = ArrayToGrid()
        grid_ds = array_ds.array_to_grid(ds, self.grid_shape)

        output_filenames = self.get_output_filenames(input_file)

        if 'CONC' in self.setup_file:
            Count_extractor = Concentrations(**self.setup_file['CONC'])
            conc_ds = Count_extractor.get_concentrations(grid_ds)
            conc_ds.to_netcdf(output_filenames['CONC'])
            conc_ds.close()

        if 'RESD' in self.setup_file:
            RESD_extractor = ResidenceTime(**self.setup_file['RESD'])
            resd_ds = RESD_extractor.get_residence_time(grid_ds)
            resd_ds.to_netcdf(output_filenames['RESD'])
            resd_ds.close()

        if 'FTLE' in self.setup_file:
            FTLE_extractor = FTLE(**self.setup_file['FTLE'])
            if self.setup_file['FTLE']['integration_time_index'] == 'all':
                FTLE_extractor.explore_ftle_timescale(grid_ds)
                grid_ds = grid_ds.drop(['x', 'y', 'z'])  # Remove duplicated.
                grid_ds.to_netcdf(output_filenames['FTLE'])  # Save all measure.
                return
            else:
                FTLE_extractor.get_ftle(grid_ds)

        if 'LCS' in self.setup_file:
            LCS_extractor = LCS(**self.setup_file['LCS'])
            LCS_extractor.get_lcs(grid_ds)

            grid_ds = grid_ds.drop(['x', 'y', 'z'])  # Remove duplicated vars
            grid_ds.to_netcdf(output_filenames['FTLE'])  # Save all measure

        # Only FTLE/LCS vars are stored in the ds to be concatenated in
        # time
        if self.ftle_LCS_only is True:
            ds_ftle_lcs = self.save_ftle_lcs_data(grid_ds)

        return ds_ftle_lcs, output_filenames

    def get_output_filenames(self, input_filename):
        output_filenames = {}
        base_filename = os.path.basename(input_filename).split('.')[0]
        if 'FTLE' in self.setup_file:
            output_filenames['FTLE'] = base_filename + '_ftle.nc'
        if 'CONC' in self.setup_file:
            output_filenames['CONC'] = base_filename + '_conc.nc'
        if 'RESD' in self.setup_file:
            output_filenames['RESD'] = base_filename + '_resd.nc'

        print('-> OUT  >>', json.dumps(output_filenames, indent=4))
        return output_filenames

    def run_ftle_lcs(self, input_file_path_pattern, output_file):
        nc_file_list = sorted(glob.glob(input_file_path_pattern), key=os.path.getmtime)
        if len(nc_file_list) == 0.:
            print('-> There is not file to process. Exiting')
            return

        if len(nc_file_list) == 1:
            print('-> INPUT >> Processing file: ', nc_file_list[0])
            self.process_one_file(nc_file_list[0], output_file)

        else:
            ds_step_list = []
            step_file_list = []
            step = 0
            for nc_ftle_field in nc_file_list:
                print('\n')
                print('-> INPUT >> Processing file:', step+1, 'of ',
                      len(nc_file_list), '>>', nc_ftle_field)
                ds_ftle, ds_ftle_filename = self.process_one_file(nc_ftle_field, str(step).zfill(3) + '.nc')
                if 'FTLE' in self.setup_file:
                    fname = ds_ftle_filename['FTLE']
                    step_file_list.append(fname)
                    ds_step_list.append(ds_ftle)
                step = step + 1
            print('\n')
            if 'FTLE' in self.setup_file:
                output_file = os.path.basename(output_file).split('.')[0] + 'ftle.nc'
                print('-> OUT  >> Merging ftle steps into:', output_file)
                xr.concat(ds_step_list, dim='time').to_netcdf(output_file)
                if os.path.exists(output_file):
                    for step_file in step_file_list:
                        os.remove(step_file)
            else:
                print('-> There is no merging of the concentration and residence times calculations \n')

        return
