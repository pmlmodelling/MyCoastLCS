import setuptools

setuptools.setup(
    name='MYCOASTLCS',
    version='0.1',
    author='Angel Daniel Garaboa Paz, Vicente Pérez Muñuzuri, James Clark, Ricardo Torres',
    author_email='angeldaniel.garaboa@usc.es',
    packages=['MYCOASTLCS','examples','docs'],
    license='GPLv3',
    url='',
    description='A python package to compute Lagrangian Coherent structres from lagrangian models',
    long_description=open('README.rst').read(),
    install_requires=[
              "numpy",
              "xarray",
              "scikit-image",
              "netcdf4"
      ],
    python_requires='>=3.6',
)
