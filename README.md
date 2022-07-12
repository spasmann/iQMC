# iQMC
Author: Sam Pasmann

iQMC is an iterative Quasi-Monte Carlo (iQMC) code for neutron transport.
The theory behind the code is outlined in [1]. iQMC uses Quasi-Monte Carlo 
methods to solve successive iterations of the standard Source Iteration for 
neutron transport.
 
## Basic Usage 

The source functions can be installed as a package. Clone the repo, pen terminal 
and within the QMC1D folder type:
```
    python -m build
```

Alternatively, simply navigate to the */problems/* folder to run one of the 
available problems. Ex:
```
    python garcia.py
```

This will run the Garcia problem as described in [2]. The parameters for initializing
the Garcia problem are written in */src/init_files/garcia_init.py*. 

This is the basic structure all problems should follow. A script in the 
*/problems/* folder and an initialization file in */src/init_files/*. If you 
want to run the two pre-defined problems you only need to edit the problem file.

The init_files create an object of the necessary parameters and can therefore be edited
in the problem file after initialization. Ex:
```
    init_data = GarciaInit()
    init_data.N = 2**10
```
A full list of initialization parameters is in the section below.

The parameters data object is then used to initialize the Source Iteration.
The Source Iteration algorithm is then executed with the *Run* function. Ex:
```
    # initialize source iteration
    SI = SourceIteration(init_data)
    # run source iteration
    SI.Run()
```
If *save_data = True* the basic simulation parameters and tally results are 
saved in a HDF5 file located in */saved_data/* and can be plotted using scripts
in */plotting/*.

For creating custom problems, a new init_file will need to be created and new
material data (cross sections) will need to be defined in */src/functions/materials.py*.
The init_file will need to contain all relevant variables listed in the section
below.


## Init File Variables

```
    N: int, Number of particles per iteration per source
    
    Nx: int, Number of spatial cells
    
    generator: string, Random number generator
    
    totalDim: int, Number of dimensions to be sampled in the problem
    
    RB: float, Right most boundary of problem
    
    LB: float, Left most boundary of problem
    
    G: int, Number of groups
    
    right: True or False, Has right boundary source
    
    left: True or False, Has left boundary source
    
    phi_right: float, strength of right boundary source
    
    phi_left: float, strength of left boundary source
    
    source: numpy array (Nx,G), strength of volumetric source
    
    material_code: string, specify cross section data
    
    avg_scalar_flux: True or False, toggle average spatial cell scalar flux tally
    
    edge_scalar_flux: True or False, toggle edge spatial cell scalar flux tally
    
    avg_agular_flux: True or False, toggle average spatial cell angular flux tally
    
    avg_current: True or False, toggle average spatial cell current tally
    
    edge_current: True or False, toggle edge spatial cell current tally
    
    save_data: True or False, save output data to HDF5 file
    
    mesh: Mesh object, create tally mesh given Nx and numpy array of RB
    
    material: Material object, create cross section data given material code
              geometry, and mesh
              
    # Optional Variables #
    
    true_flux: numpy array (Nx,G), a known analytic or benchmark solution for
               the scalar flux
```

## Source Iteration Variables

The initialization of the SourceIteration() object also contains variables which
the user may want to change after initialization and before the .Run() command.
```
    max_iter: int, maximum number of iterations 
    tol: float, desired convergence tolerance
```

## Predefined Problems

Currently there are two fully defined problems available garcia.py and multigroup.py.

**garcia.py** features an exponentialy decaying scattering cross section with space,
a left boundary source, and no volumetric source.

**multigroup.py** simulates an infinite medium multi-group problem given cross
sections from high-density polyethylene. Currently there is data available for 
a 12, 70, or 618 group problem. This problem also features a true_flux variable 
which represents the analytic solution.

## Random Number Generators

Currently the available *generator* variables are:
```
    "random" : numpy's pseudo-random number generator
    "sobol" : sobol sequence from scipy.qmc
    "halton" : halton sequence from scipy.qmc
```

## Saving Output Data

After the Source Iteration has completed, if *save_data = True*, the output data
will be stored in a HDF5 file in */saved_data/* with the following naming
convention:
```
    material_code-generator-N-Nx
```

## References 

[1] Pasmann, S., Variansyah, I., and McClarren, R. G. Convergent transport 
    source iteration calculations with quasi-monte carlo. vol. 124,
    American Nuclear Society, pp. 192–195.
    
[2] Garcia, R., Siewert, C., Radiative Transfer in Finite Inhomogeneous 
    Plane-Parallel Atmospheres, J. Quantitative Spectroscopy & Radiative 
    Transfer, 27, 2, pp. 141–148.
