# -*- coding: utf-8 -*-
""" LCS Module: This module contains the function to extract the Lagrangian
Coherent structures as ridges in the 2D FTLE field.

The detection of LCS based on the FTLE ridge extraction introduced in [a]_

A ridge of FTLE is essentially a curve where it is locally maximized in
the transverse direction, leading to a “second-derivative ridge” definition
that required that a ridge was a curve (more generally, hypersurface in high
dimensions), requiring the following conditions:

1. The first derivative of the FTLE field must be zero in the normal direction,
:math:`\mathbf{n}·\\nabla(FTLE) = 0.

2. The second derivative of the field must be negative and minimum in the
normal direction,
:math:`ε_{min}(\\nabla^2(FTLE)) < 0`,
where :math:`ε_{min}` is the minimum eigenvalue of :math:`( \\nabla ^2(FTLE))`.

3. The largest FTLE :math:`μ_{max} = μ_1` must be positive. There should be attraction
within the surface and repulsion normal to it,
:math:`μ_1 > 0 > μ_2`
"""

import xarray as xr
import numpy as np
from skimage.feature import hessian_matrix, hessian_matrix_eigvals
from skimage.measure import label


class LCS:

    def __init__(self, eval_thrsh='infer', ftle_thrsh='infer', area_thrsh=100,
                 nr_neighb=8, ridge_points_flag=False, to_dataset=True):

        self.eval_thrsh = eval_thrsh
        self.ftle_thrsh = ftle_thrsh
        self.area_thrsh = area_thrsh
        self.nr_neighb = nr_neighb
        self.ridge_points_flag = ridge_points_flag

    def get_lcs_mask_2d(self, ds):
        """

         Extract points that sit on the dominant ridges of FTLE 2D data
         A ridge point is defined as a point, where the gradient vector is
         orthogonal to the eigenvector of the smallest (negative) eigenvalue of
         the Hessian matrix (curvature).
         The found points are filtered to extract only points on strong FTLE
         separatrices. Therefore points are only taken, if:

         a) FTLE value is high at the point's position
         b) the negative curvature is high

        Args:
            - eval_thrsh (str or float, optional): scalar. Selects zones with
            small Hessian eigenvalue smaller than eval_thrsh. use 'infer' to
            obtain the threshold from the 95 percentile of the data of the
            FTLE field.
            - FTLE_thrsh (str or float, optional): scalar. Selects connected
            areas (with 4 or 8 neighbors) larger than area_thrsh. use 'infer'
            to obtain the threshold from the 95 percentile of the data of the
            FTLE field.
            - area_thrsh (float, optional):  scalar. Selects connected areas
            (with 4 or 8 neighbors) larger than area_thrsh
            - nr_neighb (int, optional): scalar. Connection neighbours (4 or 8)
            - ridge_points_flag (bool, optional): x0,y0 exact ridge poisition if 1.
            (matrix coordinates)
            - to_dataset (bool, optional): Logical mask for ridges in the FTLE
            field LCS_forward and LCS_backward to the outputted dataset.


        Example:
             Define variables
             eval_thrsh = -0.005;  # Selects zones with small Hessian eigenvalue smaller than eval_thrsh
             FTLE_thrsh = 0.07;    # Selects zones with FTLE larger than FTLE_thrsh
             area_thrsh = 10;      # Selects connected areas (with 4 or 8 neighbors) larger than area_thrsh
             nr_neighb = 8;        # conection neighbours (4 or 8)

        Returns:
            ridge_mask: logical mask for ridges in FTLE field
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

        # Make 2D hessian
        hxx, hxy, hyy = hessian_matrix(ftle, sigma=3, order='xy')

        i1, i2 = hessian_matrix_eigvals([hxx, hxy, hyy])

        Lambda2 = np.zeros_like(hxx)
        Lambda1 = np.zeros_like(hxx)
        EVecx = np.zeros_like(hxx)
        EVecy = np.zeros_like(hxx)

        for i in range(0, m):
            for j in range(0, n):
                eigvalue, eigenvec = np.linalg.eigh(
                    np.matrix([[hxx[i, j], hxy[i, j]],
                               [hxy[i, j], hyy[i, j]]]))
                Lambda2[i, j], Lambda1[i, j], = eigvalue
                EVecx[i, j], EVecy[i, j] = eigenvec[1, 0], eigenvec[1, 1]

        EVal = np.minimum(Lambda2, Lambda1)
        # EVal = np.nan_to_num(EVal)

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
        ridge_mask = FTLE_high_bin*Eval_neg_bin

        # Remove small connected areas in combined mask

        # label areas
        L = label(ridge_mask, connectivity=self.nr_neighb)
        # filter mask
        Lmax = L.max()
        for i in range(1, Lmax):
            if ((L[L == i]).size < self.area_thrsh):  # || (Lmax_ind > L_f/2)
                L[L == i] = 0

        # set new mask
        ridge_mask[L > 0] = 1
        ridge_mask[L == 0] = 0

        print('-> LCS   >> Area threshold:', self.area_thrsh)
        print('-> LCS   >> FTLE threshold:', self.ftle_thrsh)
        print('-> LCS   >> EVal threshold:', self.eval_thrsh)
        print('-> LCS   >> Potential LCS:', Lmax)

        return ridge_mask, ridge

    def to_dataset(self, ds: xr.Dataset, ridge_mask: np.array):
        """
        Add to the dataset the mask with potential LCS detected from FTLE.

        Args:
            ds (xr.Dataset): Dataset containing FTLE fields.
            ridge_mask (np.array): boolean array. True = LCS candidate

        Returns:
            None.

        """
        if 'FTLE_forward' in ds.keys():
            ds['LCS_forward'] = (ds['FTLE_forward'].dims, ridge_mask)
            ds['LCS_forward'] = ds.LCS_forward.where(~np.isnan(ds.FTLE_forward))
        elif 'FTLE_backward' in ds.keys():
            ds['LCS_backward'] = (ds['FTLE_backward'].dims, ridge_mask)
            ds['LCS_backward'] = ds.LCS_backward.where(~np.isnan(ds.FTLE_backward))

    def get_lcs_ridge_points(self, ridge_mask: np.array) -> [list, list]:
        """
        Extract the geometric x,y ridge postions from the LCS mask.

        Args:
            ridge_mask (np.array): boolean array. True = LCS candidate

        Returns:
            x_ridge (list): list of x ridge points
            y_ridge (list): list of y ridge points

        """

        m, n = ridge_mask.shape
        X = np.arange(1, n)
        Y = np.arange(1, m)
        cnt = 0
        x_ridge = []
        y_ridge = []
        for i in np.arange(0, m-1):
            for j in np.arange(0, n-2):
                if ridge_mask[i, j] == 1:
                    x = X[j:j+2]
                    z = np.array([ridge_mask[i, j], ridge_mask[i, j+1]])
                    m1 = (z[1]-z[0])/(x[1]-x[0])  # Calculamos la pendiente
                    # Calculamos el punto de corte con el eje y,
                    n1 = z[0] - m1*x[0]
                    # round zero point (-n1/m1) to the nearest j
                    if(np.abs(-n1/m1-x[0]) < np.abs(-n1/m1-x[1])):
                        js = j
                    else:
                        js = j+1
                    if(x[0] < -n1/m1 and -n1/m1 < x[1] and
                       ridge_mask[i, js] == 1):
                        cnt = cnt+1
                        x_ridge.append(-n1/m1)
                        y_ridge.append(Y[i])

        for i in np.arange(0, m-2):
            for j in np.arange(0, n-1):
                if ridge_mask[i, j] == 1:
                    y = Y[i:i+2]
                    z = np.array([ridge_mask[i, j], ridge_mask[i+1, j]])
                    m1 = (z[1]-z[0])/(y[1]-y[0])
                    n1 = z[0] - m1*y[0]
                    # round zero point (-n1/m1) to the nearest i
                    if(np.abs(-n1/m1-y[0]) < np.abs(-n1/m1-y[1])):
                        isi = i
                    else:
                        isi = i+1

                    if(y[0] < -n1/m1 and
                       -n1/m1 < y[1] and  ridge_mask[isi, j] == 1):
                        cnt = cnt+1
                        x_ridge.append(X[j])
                        y_ridge.append(-n1/m1)

                return x_ridge, y_ridge, ridge_mask

    def get_lcs(self, ds:xr.Dataset):
        """Get the LCS from a xr.Dataset containing FTLE field"

        Args:
            ds (xr.Dataset): DESCRIPTION.

        Returns:
            None.

        """

        ridge_mask, ridge = self.get_lcs_mask_2d(ds)
        self.to_dataset(ds, ridge_mask)
        if self.ridge_points_flag is True:
            self.get_lcs_ridge_points(ridge_mask, ridge)
