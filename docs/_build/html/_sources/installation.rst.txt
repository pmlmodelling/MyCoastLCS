
Requirements
------------

Requirements
    1. Data requirements:
        a. The users need to have access to outputs from a forecast system. Currently this is restricted to products available through the CMEMS platform, FVCOM based systems, NEMO based systems and soon ROMS based systems. Partner’s systems can be used provided they adhere to WP6 model output standards. 
    2. System requirements: 
        a. A multi-core Linux system (stand alone or HPC type) is recommended for running Pylag within a parallel environment due to the high computational requirements for running 105 particles. The MyCoastLCS tool doesn’t require multi-core capabilities.  The tools have been tested in PML’s local HPC (720 cpu, intel compilers), and several multi-core linux systems with Fedora versions ranging from 27 to 30. 
    3. Software requirements:
        a. MPI passage module, C++ compiler, bash terminal (or equivalent), python 3 and conda (optional). The tools’ dependencies on specific python packages and versions are all handled automatically during the installation process (i.e. Pylag-tools, PyFVCOM, scypy).



How to get the code
-------------------

To gain access to PML gitlab server where PyLag and MyCoastLCS are hosted contact Ricardo Torres (rito@pml.ac.uk) or James Clark (jcl@pml.ac.uk). 

1. PyLag: 
    a. Precompiled binaries for 64bit Linux are accessible through a conda installation archived in a private conda channel. Detailed instructions can be found here (https://gitlab.ecosystem-modelling.pml.ac.uk/PyLag/PyLag). 
    b. The original code is maintained in PML’s gitlab repository and users can gain access upon registration (https://gitlab.ecosystem-modelling.pml.ac.uk/PyLag). 
2. MyCoastLCS:
    a. Precompiled binaries are distributed with PyLag through a conda private channel (see above). 
    b. The code can be accessed through PML’s gitlab server (https://gitlab.ecosystem-modelling.pml.ac.uk/ricardo.torres/mycoast-ftle). 



Installation
------------

To gain access to PML gitlab server `PML_gitlab`_ where PyLag and MyCoastLCS are hosted contact Ricardo Torres (rito@pml.ac.uk) or James Clark (jcl@pml.ac.uk). 

Pylag
	To install there are two main options: 

	1) Install Pylag from conda-forge repositories with precompiled binaries (recomended for end-users), 

	2) Install from repository and compile the modules by yourself (recomended for developers).


	Without dependence of the option of your choice, we strongly recomended to install anaconda in your system.

	Once you installed anaconda, and before install Pylag we recomend to create an specific environment "pylag" environment with python 3.7 to install pylag:

	Option 1): install from conda-forge (end-users)

	.. code-block:: bash

		conda create -n pylag python=3.7
		conda activate pylag
		Option 1): install from conda-forge (end-users)
		conda config --append channels conda-forge
		conda config --append channels JimClark
		conda install -n pylag -c JimClark pylag

	Option 2): install from repositories (developers)

	.. code-block:: bash

		git clone https://gitlab.ecosystem-modelling.pml.ac.uk/PyLag/PyLag.git
		cd Pylag
		git pull
		conda build . 
		conda install -n pylag --use-local pylag
		pip install 


MyCoastLCS
	The code can be accessed through PML’s gitlab server, `MYCOASTLCS_https`
		.. code-block:: bash

			git clone https://gitlab.ecosystem-modelling.pml.ac.uk/ricardo.torres/mycoast-ftle.git
			cd mycoast-ftle 
			git checkout --track origin/T_direction_styling
			python3.7 -m pip install .

	.. PML_gitlab: //gitlab.ecosystem-modelling.pml.ac.uk/
	.. MYCOASTLCS_https: //gitlab.ecosystem-modelling.pml.ac.uk/ricardo.torres/mycoast-ftle). 

	How to install the code

	    2. MyCoastLCS can be installed through conda with PyLag or as a site package with python. Instructions are included with the downloaded package. 


how to run MYCOASTLCS 
---------------------

At the time of writing the tools are run via command line instructions from within a bash terminal (or equivalent). Once all the tools are installed, two configuration files are needed (examples provided by PML): one for the PyLag simulations and one for MyCoastLCS setup and plotting. Other configurable templates include an HPC queue submission script (if using an HPC system to run PyLag, PML can help if needed) and the location of rivers and/or Waste Water source points (Waste Water Treatment Plants and Surface Drainage Outfalls) as a shapefile. 

.. image:: MyCoastLCS_scheme.png

The steps a user would follow include:

1. Generate a dense regularly spaced set of initial conditions. A python script is provided for this purpose. 

2. Generate a set of initial positions for point source locations of pollutants. The python script provided as an example considers 3 types of sources points: rivers, waste water treatment plants and surface drainage outfalls. Ideally, the source concentrations would vary depending on river flow, water treatment plant capacity and rain fall. These dependencies are not implemented in our example. 

3. Generate master configuration file to specify type of simulation such as forward or backward integration (PyLag), duration of simulation for FTLE integration (i.e. 6-12 hours), total length of simulation (depends on user’s operational system forecast window), location of pollution point sources and  plots to be generated among others.. 

4. Run Pylag with output from your forecast model of choice to generate the files to be used by MyCoastLCS and plotting scripts. While PyLag reads netcdf files that can be accessed from a local folder or a THREDDS or ERDDAP server via OpendAP it is recommended that the files are downloaded first for performance benefits. 

5. Run MyCoastLCS to obtain FTLE and LCS fields. The tool can optionally be used to estimate concentrations and residence times that can be used to generate a qualitative water quality indicator. Depending on the configuration of the workflow, a set of standard plots will be generate for your domain.
