![alt text](https://raw.githubusercontent.com/spasmann/iQMC/main/post_process/figures/iQMC_logo_2.png)

[![License](https://img.shields.io/github/license/spasmann/iQMC)
[![Code of Conduct](https://img.shields.io/badge/code%20of%20conduct-contributor%20covenant-green.svg)](https://www.contributor-covenant.org/version/2/1/code_of_conduct/code_of_conduct.md)

Author: Sam Pasmann

This codebase serves as a testbed for the iterative Quasi-Monte Carlo (iQMC) method for neutron transport.
The theory behind iQMC is outlined in [1]. iQMC uses Quasi-Monte Carlo 
methods to solve successive iterations of the Source Iteration, Power Iteration and other advanced
linear solvers for neutron transport.
 
## Basic Usage 

Predefined problems may be run by editing and running any of the files in */scripts/*.
Adding new problems requires a new input and material file.


## Random Number Generators (MC) / Low Discrepency Sequences (QMC)

Currently the available *generator* variables are:
```
    "random" : numpy's pseudo-random number generator
    "sobol" : sobol sequence from scipy.qmc
    "halton" : halton sequence from scipy.qmc
```


## References 

[1] Pasmann, S., Variansyah, I., and McClarren, R. G. Convergent transport 
    source iteration calculations with quasi-monte carlo. vol. 124,
    American Nuclear Society, pp. 192–195.
    
[2] Garcia, R., Siewert, C., Radiative Transfer in Finite Inhomogeneous 
    Plane-Parallel Atmospheres, J. Quantitative Spectroscopy & Radiative 
    Transfer, 27, 2, pp. 141–148.
    
## Random Joke
![Jokes Card](https://readme-jokes.vercel.app/api)
