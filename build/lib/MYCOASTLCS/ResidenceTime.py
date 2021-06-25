import numpy as np
import xarray as xr


class Concentrations:
    
    def __init__(self, bins_option, nbins):
        
        
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
    
    
    def get_bins(self,base):
        
        if self.bins_option == 'domain':
            x_bins = np.linspace( base.x.min(), base.x.max(), self.nbins)
            y_bins = np.linspace( base.y.min(), base.y.max(), self.nbins)
        elif self.bins_option == 'origin':
            if self.nbins:
                x_bins = np.linspace( base.x0.values.min(),base.x0.values.max(), self.nbins)
                y_bins = np.linspace( base.y0.values.min(),base.y0.values.max(), self.nbins)
            else:
                 x_bins = base.x0.values
                 y_bins = base.y0.values
        elif self.bins_option == 'personalized':
            x_bins = np.linspace(*self.nbins[0])
            y_bins = np.linspace(*self.nbins[1])
                 
        self.bins = (x_bins, y_bins)

        x_centers = (x_bins[:-1] + x_bins[1:]) / 2.
        y_centers = (y_bins[:-1] + y_bins[1:]) / 2.
        
        self.centers = [x_centers,y_centers]

        print('**********************************')
        print('Computing concentraions in domain:')
        print('**********************************')
        print('-> x :',x_bins.min(),' ',x_bins.max())
        print('-> y :',y_bins.min(),' ',y_bins.max())
        print('-> nbins:',self.nbins)
        
        return
                       

    def get_concentrations(self, base, to_dataset=True):
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
        
        concentrations = np.zeros((self.ds.time.size, self.centers[1], self.centers[2]))

        for i in range(0, base.time.size):
            concent, _, _ = np.histogram2d(base.x.isel(time=i).values.ravel(), base.y.isel(time=i).values.ravel(), bins=self.bins)
            concentrations[i, :, :] = concent.transpose()

        coords = {'time': base.time, 'y_c': ('y_c', self.centers[1]), 'x_c': ('x_c', self.centers[2])}
        dims = ['time', 'y_c', 'x_c']
        
        if to_dataset == True:
            self.da_output['concentrations'] = xr.DataArray(concentrations, coords, dims=dims)
        else:
            return concentrations

        
    def get_netcdf(self,file_output):
        self.ds_output = xr.Dataset(self.da_output)
        self.ds_output.to_netcdf(file_output)



