# 1DQMC
Author: Sam Pasmann

1DQMC is a hyrbid Quasi-Monte Carlo (QMC) and deterministic code for neutron transport.
The theory behind the code is outlined in [1]. 1DQM uses Quasi-Monte Carlo 
methods to solve successive iterations of the standard Source Iteration for 
neutron transport. In the future, iterative methods like
Krylov Subspace methods will also be added as the base QMC code allows for 
plug-and-play compatibility with different solvers.
 
## Basic Usage 

To use 1DQMC, simply clone the repo to your local machine and navigate to the
*/problems/* folder to run one of the available problems. Ex:
```
    python garcia.py
```

This will run the Garcia problem as described in [2]. The parameters for initializing
the Garcia problem are written in */src/init_files/garcia_init.py*. 

This is the basic structure all problems should follow. A script in the 
*/problems/* folder and an initialization file in */src/init_files/*. If you 
wants to run the three pre-defined problems you only need to edit the problem file.
The init_files create an object of the necessary parameters and can therefore be edited
in the problem file after initialization. Ex:
```
    init_data = GarciaInit()
    init_data.N = 2**10
```
A full list of initialization parameters is below.

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

- N: Number of particles per iteration per source
- Nx: Number of spatial cells
- generator: Random number generator


        self.N = N
        self.Nx = Nx
        self.generator = generator
        self.totalDim = 3
        self.RB = 5
        self.LB = 0
        self.G = 1
        self.right = False
        self.left = True
        self.phi_left = 1.0
        self.source = np.zeros((self.Nx,self.G))
        self.material_code = "garcia_data"
        self.geometry = "slab"
        self.avg_scalar_flux = True
        self.edge_scalar_flux = False
        self.avg_angular_flux = False
        self.avg_current = False
        self.edge_current = False
        self.shannon_entropy = False
        self.save_data = True
        self.mesh = Mesh(self.Nx, np.array((self.RB,)))
        self.material = Material(self.material_code, self.geometry, self.mesh)
        
### Tallies

## Source Iteration Variables

## Available Problems

## Random Number Generators

## Saving Output Data

## Citations

### Dependencies
- Numpy
- Scipy
- h5py



[1] Pasmann, S., Variansyah, I., and McClarren, R. G. Convergent transport 
    source iteration calculations with quasi-monte carlo. vol. 124,
    American Nuclear Society, pp. 192–195.
    
[2] 