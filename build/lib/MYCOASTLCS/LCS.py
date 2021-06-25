import xarray as xr
import numpy as np
from skimage.feature import hessian_matrix, hessian_matrix_eigvals
from skimage.measure import label

        
class LCS:

    def __init__(self, eval_thrsh='infer', ftle_thrsh='infer', area_thrsh=100, nr_neighb=8, ridge_points=False, to_dataset=True):

        self.eval_thrsh = eval_thrsh
        self.ftle_thrsh = ftle_thrsh
        self.area_thrsh = area_thrsh
        self.nr_neighb = nr_neighb
        self.ridge_points = ridge_points

    def get_lcs_mask_2d(self,ds):
    
        """

         Extract points that sit on the dominant ridges of FTLE 2D data
         A ridge point is defined as a point, where the gradient vector is orthogonal
         to the eigenvector of the smallest (negative) eigenvalue of the Hessian
         matriv (curvature).
         The found points are filtered to extract only points on strong FTLE
         separatrices. Therefore points are only taken, if:

         a) FTLE value is high at the point's position
         b) the negative curvature is high
       
        
        Args:
            - eval_thrsh (str or float, optional): scalar. Selects zones with small Hessian eigenvalue smaller than eval_thrsh. use 'infer' to obtain the threshold from the 95 percentile of the data of the FTLE field.
            - FTLE_thrsh (str or float, optional): scalar. Selects connected areas (with 4 or 8 neighbors) larger than area_thrsh. use 'infer' to obtain the threshold from the 95 percentile of the data of the FTLE field.
            - area_thrsh (float, optional):  scalar. Selects connected areas (with 4 or 8 neighbors) larger than area_thrsh
            - nr_neighb (int, optional): scalar. Connection neighbours (4 or 8)
            - ridge_points (bool, optional): x0,y0 exact ridge poisition if 1. (matrix coordinates)  
            - to_dataset (bool, optional): Logical mask for ridges in the FTLE field LCS_forward and LCS_backward to the outputted dataset.

            
        Example:
             Define variables
             eval_thrsh = -0.005;        # Selects zones with small Hessian eigenvalue smaller than eval_thrsh
             FTLE_thrsh = 0.07;          # Selects zones with FTLE larger than FTLE_thrsh
             area_thrsh = 10;            # Selects connected areas (with 4 or 8 neighbors) larger than area_thrsh
             nr_neighb = 8;              # conection neighbours (4 or 8)

        Returns:
            combined_bin: logical mask for ridges in FTLE field
            x0: Positions of points on FTLE ridges
            y0: Positions of points on FTLE ridges
        """

        if 'FTLE_forward' in ds.keys():
            ftle = ds['FTLE_forward'].fillna(0).squeeze().values
            m, n = ftle.shape
        elif 'FTLE_backward' in ds.keys():
            ftle = ds['FTLE_backward'].fillna(0).squeeze().values
            m, n = ftle.shape

        if self.ftle_thrsh == 'infer':
            self.ftle_thrsh = np.percentile(ftle, 95)

        # Gradient and Hessian matrix (2nd derivatives) from finite differences
        [dy, dx] = np.gradient(ftle)

        # Eigenvalues of Hessian matrix (analytically)
        # EVal = min( 1/2*(a+d) + sqrt(b.*c + 4*(a-d).^2), 1/2*(a+d) - sqrt(b.*c + 4*(a-d).^2) );
        # psu = (a-EVal-c)./(d-EVal-b);
        # Smaller (negative) Eigenvector of Hessian matrix (analytically)
        # EVecx = 1./sqrt(1^2 + psu.^2);
        # EVecy = psu./sqrt(1^2 + psu.^2);

        # Make 2D hessian
        hxx, hxy, hyy = hessian_matrix(ftle, sigma=3,order='xy')

        i1, i2 = hessian_matrix_eigvals([hxx, hxy, hyy])

        Lambda2 = np.zeros_like(hxx)
        Lambda1 = np.zeros_like(hxx)
        EVecx = np.zeros_like(hxx)
        EVecy = np.zeros_like(hxx)

        for i in range(0, m):
            for j in range(0, n):
                eigvalue, eigenvec = np.linalg.eigh(
                    np.matrix([[hxx[i, j], hxy[i, j]], [hxy[i, j], hyy[i, j]]]))
                Lambda2[i, j], Lambda1[i, j], = eigvalue
                EVecx[i, j], EVecy[i, j] = eigenvec[1, 0], eigenvec[1, 1]

        EVal = np.minimum(Lambda2, Lambda1)
        EVal = np.nan_to_num(EVal)

        # Compute the direction of the minor eigenvector
        # Define ridges as zero level lines of inner product 
        # of gradient and eigenvector of smaller (negative) eigenvalue
        ridge = EVecx*dx + EVecy*dy

        
        if self.eval_thrsh == 'infer':
            self.eval_thrsh = np.percentile(EVal, 95)

        # Define filter masks with high negative curvature
        Eval_neg_bin = EVal < self.eval_thrsh

        # High FTLE value
        FTLE_high_bin = ftle > self.ftle_thrsh

        # Combined mask
        combined_bin = FTLE_high_bin*Eval_neg_bin

        # Remove small connected areas in combined mask

        # label areas
        L = label(combined_bin, connectivity=self.nr_neighb)
        # filter mask
        Lmax = L.max()
        for i in range(1, Lmax):
            if ((L[L == i]).size < self.area_thrsh):  # || (Lmax_ind > L_f/2)
                L[L == i] = 0

        # set new mask
        combined_bin[L > 0] = 1
        combined_bin[L == 0] = 0
        
        print('-> LCS   >> Area threshold:', self.area_thrsh)
        print('-> LCS   >> FTLE threshold:', self.ftle_thrsh)
        print('-> LCS   >> EVal threshold:', self.eval_thrsh)
        print('-> LCS   >> Potential LCS:', Lmax)

        
        return combined_bin, ridge

    
    def to_dataset(self, ds, combined_bin): 
        if 'FTLE_forward' in ds.keys():
            ds['LCS_forward'] = ( ds['FTLE_forward'].dims, combined_bin)
        elif 'FTLE_backward' in ds.keys():
            ds['LCS_backward'] = (ds['FTLE_backward'].dims, combined_bin)


    def get_lcs_ridge_points(self,combined_bin,ridge):
        x0 = 0
        y0 = 0

        m, n = ridge.shape
        X = np.arange(1, n)
        Y = np.arange(1, m)
        cnt = 0
        x0 = []
        y0 = []
        for i in np.arange(0, m-1):
            for j in np.arange(0, n-2):
                if combined_bin[i, j] == 1:
                    x = X[j:j+2]
                    z = np.array([ridge[i, j], ridge[i, j+1]])
                    m1 = (z[1]-z[0])/(x[1]-x[0])  # Calculamos la pendiente
                    # Calculamos el punto de corte con el eje y,
                    n1 = z[0] - m1*x[0]
                    # round zero point (-n1/m1) to the nearest j
                    if(np.abs(-n1/m1-x[0]) < np.abs(-n1/m1-x[1])):
                        js = j
                    else:
                        js = j+1
                    if(x[0] < -n1/m1 and -n1/m1 < x[1] and combined_bin[i, js] == 1):
                        cnt = cnt+1
                        x0.append(-n1/m1)
                        y0.append(Y[i])

        for i in np.arange(0, m-2):
            for j in np.arange(0, n-1):
                if combined_bin[i, j] == 1:
                    y = Y[i:i+2]
                    z = np.array([ridge[i, j], ridge[i+1, j]])
                    m1 = (z[1]-z[0])/(y[1]-y[0])
                    n1 = z[0] - m1*y[0]
                    # round zero point (-n1/m1) to the nearest i
                    if(np.abs(-n1/m1-y[0]) < np.abs(-n1/m1-y[1])):
                        isi = i
                    else:
                        isi = i+1

                    if(y[0] < -n1/m1 and -n1/m1 < y[1] and combined_bin[isi, j] == 1):
                        cnt = cnt+1
                        x0.append(X[j])
                        y0.append(-n1/m1)

                return x0, y0, combined_bin
        
    def get_lcs(self,ds):
        combined_bin, ridge = self.get_lcs_mask_2d(ds)
        self.to_dataset(ds,combined_bin)
        if self.ridge_points == True:
           self.get_lcs_ridge_points(combined_bin,ridge)
                 
          
        
    

   