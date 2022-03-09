# 1DQMC
Author: Sam Pasmann

1DQMC is a hyrbid Quasi-Monte Carlo (QMC) and deterministic code for neutron transport.
The theory behind the code is outlined in [1]. 1DQM uses Quasi-Monte Carlo 
methods to solve successive iterations of the standard Source Iteration for 
neutron transport. In the future, iterative methods like
Krylov Subspace methods will also be added as the base QMC code allows for 
plug-and-play compatibility with different solvers.
 
## Basic Usage 

To use 1DQMC, simply clone the repo to your local 

### Dependencies
- Numpy
- Scipy
- h5py

## Available Problems

## Random Number Generators

## Tallies

## Init File Variables

## Source Iteration Variables

## Saving Output Data

## Citations



1-Dimensional Quasi-Monte Carlo code for neutron transport.

Source code is located in */src/functions*

*/src/init_files/* contains problem-specific initialization parameters


[1] Pasmann, S., Variansyah, I., and McClarren, R. G. Convergent transport 
    source iteration calculations with quasi-monte carlo. vol. 124,
    American Nuclear Society, pp. 192–195.