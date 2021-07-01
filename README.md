# h2p-ricpad

An implementation of the Riccati-Padé method to solve the Schrödinger Equation
for the hydrogen molecule ion.

# Installation

This program is designed to run on GNU/Linux systems. 
For building it requires cmake, a working C++ compiler, and, additionally, 
the following libraries:

* `Boost`
* `Eigen`
* `Ginac`
* `GMP`
* `MPFR`

## Instructions for installation in Ubuntu and derivates

First, install the dependencies with

    sudo apt install build-essential cmake libeigen3-dev libboost-all-dev libgmp-dev libmpfr-dev libginac-dev

Then download this program with:

    git clone https://github.com/javierelpianista/h2p-ricpad
    
(this step can be avoided if downloaded from a repository).
    
A folder called `h2p-ricpad` will be created in the current directory.
CD to that folder, then follow the typical procedure to build with cmake

    mkdir build
    cd build
    cmake ..
    cmake --build .

If succesful, an executable named `h2p-ricpad` will be created in this folder.

# Usage

    Usage: h2p-ricpad MODE OPTIONS
    
    MODE can be either 'minimum' for computing U, A, and  R simultaneously for 
    the equilibrium internuclear distance,  'fixed', for computing U and A for 
    a given R,  'coupling', for computing the coupling constant A, given U and R,  
    'fixed-U', for computing U given A and R.
    Required options for two- and three-parameters methods:
      --Dmin arg            Minimum D value
      --U0 arg              Initial value of the electronic+nuclear energy
      --A0 arg              Initial value of the coupling constant
      --R0 arg              Initial value of the internuclear distance
    
    Non-mandatory options:
      --help                    Print this message
      --Dmax arg (=-1)          Maximum D value. If Dmax == -1, Dmax is assumed to 
                                be infinite
      --d arg (=0)              Value of d
      --double-D-lambda         Double the value of D for the lambda equation. 
                                Useful at small R.
      --double-D-mu             Double the value of D for the mu equation. Useful 
                                at large R.
      --m arg (=0)              Value of the quantum number associated with the 
                                angular part of the spheroidal equations
      --s arg (=0)              Value stating if the second spheroidal equation's 
                                eigenfuncions are even (s=0) or odd (s=1)
      --ndigits arg (=1000)     Number of digits for the numerical calculations
      --tol arg (=1E-100)       Tolerance for the Newton-Raphson method
      --h arg (=1E-200)         Step size for the Newton-Raphson method
      --hd arg (=1E-400)        Step size for the Newton-Raphson for the 
                                computation of numerical differentiation outside 
                                the Newton-Raphson iterations (Set this to at least
                                sqrt(h))
      --use-E                   Set this option if the provided value of U0 is the 
                                electronic energy, instead of the electronic + 
                                nuclear one
      --no-log                  Set this option if you don't want h2p to output to 
                                any files
      --log-nr                  Set this option to print out each Newton-Raphson 
                                iteration
      -o [ --output-file ] arg  Output file
      --nr-max-iter arg (=100)  Maximum number of Newton-Raphson iterations
      --target-digits arg (=-1) Target number of digits. If this number is reached,
                                the program stops.
    
# Examples

The folder `examples` contains simple scripts that show how to use `h2p-ricpad`.
Before using them, the `h2p-ricpad` executable should be either copied to or linked to
it. If the instructions for installation were followed, then this should suffice.

    cd examples
    ln -s ../build/h2p-ricpad .
    sh req.sh

The following files are provided:

* `req.sh`: computation of the equilibrium internuclear distance, electronic+nuclear energy, and separation constant for the ground state
* `UA.sh`: computation of the electronic+nuclear energy and separation constant for the ground state at a fixed internuclear distance.

# Computation of potential energy curves

The Python3 script `sweep.py` that allows to use the `h2p-ricpad` program to compute the potential energy curve for several states.
In order to use it, the user should open it and modify it according to their setup and preferences.
It requires the libraries `numpy` and `mpmath`.

It can be used by writing
    python sweep.py <filename>

Where filename is a `.csv` file with a list of sets of `l, m, I` numbers.
File `examples/sweep_example.csv` is provided as an example input file.
