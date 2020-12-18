# Usage

h2p MODE OPTIONS

Where MODE can be either `fixed` or `minimum`. There are some optional 
and mandatory options, depending on the mode chosen.
Both modes have options in common:

## Required options
* Dmin
* U0: initial value of the electronic + nuclear energy
* A0: initial value of the coupling constant
* R0: in fixed mode, R0 represents the current value of the internuclear 
distance, whereas in minimum mode it is used as the initial value for the 
computations.

## Options with default values
* d (= 0)
* Dmax (= -1)
* m (= 0)
* s (= 0)
* ndigits (= 1000)
* tol (= 1E-100)
* h (= 1E-200)
* h2 (= 1E-400)
