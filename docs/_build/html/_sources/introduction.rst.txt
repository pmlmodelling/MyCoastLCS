.. highlight:: rst


MYCOASTLCS-Summary
===================

Introduction
------------
MyCoastLCS has been developed to enable the visualisation, identification and quantification of pollution hotspots of non-reactive, buoyant, slow diffusing and short lived substances. These assumptions can approximate the behaviour of plastic debris as well as sewage waste over the time scale of a few days. The tools specifically calculates the spatial distribution and time evolution of Finite Time Lyapunov Exponent (FTLE) fields, extracts Lagrangian Coherent Structures (LCS) and estimates spatially discretised time evolving concentrations and residence times. 

The tools have been developed by Plymouth Marine Laboratory (PML) and Universidade de Santiago de Compostela (USC). The system comprises of a python/cython package that includes a lagrangian particle simulator (Pylag), the MyCoastLLCS tool and a series of scripts to generate initial seeding locations for PyLag and to plot the results and compile a set of fixed  graphs into a html and/or pdf file. 
We provide PyLag to enable the generation of lagrangian trajectories required to calculate the Finite Time Lyapunov Exponents and pathways and concentrations of point source pollutants. However, MyCoastLCS tool can be used with other lagrangian models provided they produce a netcdf file with variables that include particle positions at discrete time intervals in 3 dimensions and associated particle ids. MyCoastLCS has so far been tested with LAGAR and PyLag. 

The Finite Time Lyapunov exponentes and the Lagrangian Coherent Structures are computed using the  `ridge detection method: <https://shaddenlab.berkeley.edu/uploads/LCS-tutorial/LCSdef.html>`_ for PYLAG/LAGAR Lagrangian output models. 

At the moment, this module support the following features for Pylag and LAGAR models:
 
- Computation of FTLE in 2d and 3d cartesian coordinates
- Computation of FTLE in 2d spherical cooordinates.
- Extraction of LCS using FTLE ridge method for 2d FTLE.

The module can be used with other lagrangian models. It requires a netCDF dataset alias dictionary to provided the link between the lagrangian model variable names and the internal variables x,y,z used to compute the FTLE and other Lagrangian measures. The advection performed (in order to compute the FTLE) should be done with a grid of initial conditions regularly spaced in 2d and/or 3d.

Authoring
---------
.. image:: mycoast_logo.png
	:scale: 65 %
.. image:: usc_logo.jpg
	:scale: 25 %
.. image:: pml_logo.png
	:scale: 65 %

.. bibliographic fields (which also require a transform):

| **Authors:** Dr. Ángel Daniel Garaboa Paz, Prof. Vicente Pérez Muñuzuri, Dr. James Clark and Dr. Ricardo Torres. 
| **Organization:** University of Santiago de Compostela, Plymouth Marine Laboratory. 
| **Date:** March-2019. 
| **Status:** Alpha release. 
| **Revision:**: 0000. 
| **Version:** 0.1. 
| **Copyright:** GPLv3. 


Plattform
---------
| **Language:** Python 3
| **OS:** Linux or Windows
| **Requirements:** Anaconda / pip
| **Python packages required:** 
| - xarray
| - numpy
| - scipy
| - tqdm (fancy progress bars)
