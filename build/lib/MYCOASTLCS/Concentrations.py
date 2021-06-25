import numpy as np
import xarray as xr
from MYCOASTLCS import base


class Concentrations:
    
    def __init__(self, base, bins_option, nbins):
        
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
        self.nbins = nbins
        self.bins_option = bins_option
        self.bins = []
        
        # Here we stablish the correspondece between the Lagrangian model input and also
        # the private variables x,y,z of the MYCOASTFTLE module.
        # Each model produces the outputs with their internal variable names.
        # The alias dictionary "links" those variables with the
        # internal x,y,z of the MYCOAST_FTLE module.
        # There are Two defaults dictionaries for LAGAR and PYLAG.
        if self.nbins:
            self.nbis = nbins
        else:
            self.nbins = 100
                       
            
    def get_bins(self):
        """
        It computes the raw residence time. 
        This function counts the timesteps that a particle spetn in a box.
        It doesn't matter if it is just a particle or a bunch of them. 


        Args:
            timeindex (TYPE, optional): Description
            bins_option (str, optional): Description
            nbins (int, optional): Description
            to_dataset (bool, optional): Description

        Returns:
            TYPE: Description

        """
        if self.bins_option == 'domain':
            x_bins = np.linspace(self.x.min(), self.x.max(), self.nbins)
            y_bins = np.linspace(self.y.min(), self.y.max(), self.nbins)
        elif self.bins_option == 'origin':
            if self.nbins:
                x_bins = np.linspace(self.x0.values.min(),self.x0.values.max(),self.nbins)
                y_bins = np.linspace(self.y0.values.min(),self.y0.values.max(),self.nbins)
            else:
                 x_bins = self.x0.values
                 y_bins = self.y0.values
        elif self.bins_option == 'personalized':
            x_bins = np.linspace(*self.nbins[0])
            y_bins = np.linspace(*self.nbins[1])
                 
        self.bins = (x_bins, y_bins)

        print('**********************************')
        print('Computing concentraions in domain:')
        print('**********************************')
        print('-> x :',x_bins.min(),' ',x_bins.max())
        print('-> y :',y_bins.min(),' ',y_bins.max())
        print('-> nbins:',self.nbins)
        return
    
    def get_centers(self):
        self.centers = [(self.bins[0][:-1] + self.bins[0][1:]) / 2., (self.bins[1][:-1] + self.bins[1][1:]) / 2.]
        
        
    def get_concentrations(self, to_dataset=True):
        """

        It computes the raw residence time. 
        This function counts the timesteps that a particle spetn in a box.
        It doesn't matter if it is just a particle or a bunch of them. 


        Args:
            timeindex (TYPE, optional): Description
            bins_option (str, optional): Description
            nbins (int, optional): Description
            to_dataset (bool, optional): Description

        Returns:
            TYPE: Description


        """
        
        nt,ny,nx = self.ds.time.size,self.centers[1].size,self.centers[2].size
        concentrations = np.zeros((nt,ny,nx))

        for i in range(0, self.time.size):
            concent, _, _ = np.histogram2d(self.base.x.isel(time=i).values.ravel(), self.base.y.isel(time=i).values.ravel(), bins=self.bins)
            concentrations[i, :, :] = concent.transpose()

        coords = {'time': self.time, 'y_c': ( 'y_c', self.y_centers), 'x_c': ('x_c', self.x_centers)}
        dims = ['time', 'y_c', 'x_c']
        
        if to_dataset == True:
            self.da_output['concentrations'] = xr.DataArray(concentrations, coords, dims=dims)
        else:
            return concentrations

        
    def get_netcdf(self,file_output):
        self.ds_output = xr.Dataset(self.da_output)
        self.ds_output.to_netcdf(file_output)



