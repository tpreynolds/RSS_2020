# Continuous Funnel Generation Algorithm

This code accompanies the paper submitted to the RSS 2020 workshop titled "PAPER_TITLE" by Taylor P. Reynolds, Danylo Malyuta, Mehran Mesbahi and Behcet Acikmese from the University of Washington.

## Basics

Using a simple 2D planar quadrotor model, this code computes a funnel around a pre-set nominal trajectory in order to demonstrate the use of temporally-interpolated quadratic Lyapunov functions. 

## Installation Guide

To use this code, you will need to clone the repository and include the submodule YALMIP. Using
`git clone --recurse-submodules git@github.com:tpreynolds/RSS_2020.git`
should do the trick. *You will need to install an additional solver to use the code*. 

### Solver Specifics
The code was developed using MOSEK, which can be [installed from here](https://www.mosek.com/downloads/), and free academic licenses are available. 

Depending on where the solver is installed on your machine, you may need to edit the `utils/set_path.m` file. Open this file and change the value of `solverroot` on line 27 to be wherever you've just installed the solver you intend to use. The picture below is my setup for MOSEK. 
![alt text][solverroot]

[solverroot]: https://github.com/tpreynolds/RSS_2020/blob/master/figs/solverroot.png

*If you are not using MOSEK*, add a line to the `main.m` file (anywhere before line 28) that says `cfga.opts.solver = 'solvername';`. The `solvername` here must match whatever YALMIP calls the solver that you intend to use.  

## License 

This repository is free of charge to use and is openly distributed, but note that

1. The following reference must be used if this code is used in a published work, or in preparation/validation of the results in a published work.:
> Taylor P. Reynolds, Danlyo Malyuta, Mehran Mesbahi and Behcet Acikmese, "PAPER_TITLE," 2nd RSS Workshop on Robust Autonomy, July 2020.

2. Use of this repository is subject to the copyright and licenses of YALMIP and any solvers used.

3. Forks or versions of this code may not be re-distributed as a part of a commercial product unless agreed upon with the authors.

4. This repository is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

5. Forks or versions of this repository must include, and follow, this license in any distribution. 