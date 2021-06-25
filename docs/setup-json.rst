.. highlight:: rst

main Module
===========

It runs the core :file:`Common.py` of the MYCOASTLCS package. It reads the setup file, runs and connects the different inputs and outputs of the different MYCOASTLCS modules.


Reading the config file
-----------------------

It reads the JSON file described in `LAGAR Setup`_ passed as a argument from the command line,

::

    $ python -m MYCOASTLCS -j setup.json -i 'Pylag_0*.nc' -o output.nc


Setup json template
===================


Before running, the first step is to configure the :file:`template_setup.json`. This file contains different
dictionaries to be read by the MYCOASTLCS main method. Each dictionary inside the JSON file configures a module or class.



Setup structure
------------------

.. inheritance-diagram:: MYCOASTLCS.common
   :parts: 1


The dictionary controls manages two main things: 
    1) Link between the grid of particles by your Lagrangian Model and MYCOASTLCS. 
    2) The desire measures to compute.

::

    {
    "common":{
        "model":"pylag",
        "alias":{"xpos":"x",
                 "ypos":"y",
                 "zpos":"z",
                 "time":"time",
                 "particles":"particle_id"
                },
        "grid_shape":[1,500,500]
    },
    "FTLE":{
        "spherical_flag":0,
        "integration_time_index":-1
    },
    "LCS":{
        "eval_thrsh":"infer",
        "ftle_thrsh" : "infer",
        "area_thrsh" : 100,
        "nr_neighb" : 2,
        "ridge_points_flag" : 0
    },
    "CONC":{
        "bins_option":"origin",
        "nbins":[]
    },
    "RESD":{
        "bins_option":"origin",
        "nbins":[]
    }
    }


Common - keys
--------------

The **common** dictionary ontrols the correspondece between the Lagrangian model output variables and the x,y,z of the MYCOASTFTLE module.
Each Lagrangian model produces the outputs with their internal variable names. 

-**model**: String  with ame of the Lagrangian model used. At this moment just Pylag and LAGAR are supported.

- **alias**: In case that your model is not supported, you should specify the link between your model output variables and the MYCOAST internal variables. An example, if you module produces, *xpos, ypos, zpos, time, particles*, your dictionary should looks like:


        ::

            "alias":{"xpos":"x",
                     "ypos":"y",
                     "zpos":"z",
                     "time":"time",
                     "particles":"particle_id"
                    },


- **grid_shape**: If your model produces outputs with where the set of particles has a unique index dimension (it is flatted), for example: :math:`\mathbf{r}_i`, you should specify the *grid dimension* initially used: :math:`\mathbf{r}_{ijk}, i=1:n_z, j=1:n_y, k:1:n_x`.

            ::

                "grid_shape":[1,500,500]

In this example, our grid of initial condition, has 250000 particles: 1 depth layer (2D :math:`n_z=1`) and :math:`n_x*n_y = 500 x 500` particles on horizontal.

FTLE - keys
------------
We have two options to compute the FTLE. 

- **spherical_flag**: Flag to set when the coordinates of your particle positions are given in degrees(1) or meters(0).
- **integration_time_index**: It sets which timestep of your each netCDF input dataset consider to compute the FTLE. It follows Python convention, so -1 means to use last timestep with the final particle positions. 


LCS - keys
----------
The LCS keys sets the LCS extractions parameters from FTLE.

- **eval_thrsh**: Eigenvalue threshold. It sets the thresholds for the minimum eigenvalue of the hessian. More negative means more curved or more folded FTLE ridge. `MYCOASTLCS.LCS`. If you set to `infer`, it takes the 95 percentile of minimum eigenvalues of the hessian from the FTLE field.

- **ftle_thrsh**: Ftle threshold. It sets the thresholds to filter weak FTLE values. Those values below that threshold cannot be considered as LCS. This filters low strechting manifolds. If you set to `infer`, it takes the 95 percentile of the values from the FTLE field to perform the filtering.

- **area_thrsh**: Area threshold in pixels. It sets the number of points in *r_ijk*/"grid_shape" that should be contiguous in order to be considered as "large" structure enough. Strongly dependent of **nr_neighb**.

- **nr_neighb** : Nearest neighboring (2 or 4). It sets the way to consider points contiguos FTLE points. It2 you just need a point and a continguos point.


