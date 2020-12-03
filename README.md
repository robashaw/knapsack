## Knapsack 1.0.0

[![Build Status](https://dev.azure.com/robertshaw383/knapsack/_apis/build/status/robashaw.knapsack?branchName=master)](https://dev.azure.com/robertshaw383/knapsack/_build/latest?definitionId=1&branchName=master)
[![codecov](https://codecov.io/gh/robashaw/knapsack/branch/master/graph/badge.svg)](https://codecov.io/gh/robashaw/knapsack)

Knapsack is a FORTRAN 08 package that implements multiple algorithms for the efficient enumeration of bosonic configurations subject to an energy criterion. These configurations can then be used to determine non-radiative rates (internal conversion) within the Franck-Condon and Hertzberg-Teller regimes. More detailed documentation will be available shortly, but in the meantime, example inputs for calculating rates can be found in the examples folder.

## Dependencies

- A modern FORTRAN compiler, at least FORTRAN 08 compliant - you may need to provide flags to specify the FORTRAN version. This has been tested with:
  * gfortran (GCC v10.2 and above)
  * ifort (v19.2 and above)
- CMake/CTest build tools (v3.12 and higher)

The assumed input format and units of gradients (for Franck-Condon rates) and the hessian (for Hertzberg-Teller rates) is that outputed by Gaussian 16; more details can be found in the documentation.

## Installation

To clone the repo and the BLT submodule (used for building and testing), do one of the following:
```
git clone --recurse-submodules [repo URL]
```
or clone the repo normally then in the empty BLT directory run
```
git submodule init
git submodule update
```

Then from the main repo directory, run the following:
```
cmake -Bbuild .. [OPTIONS]
cd build
make [OPTIONS]
```
To then run the test suite and install, use
```
make test
make install
```

## Documentation

Coming soon.

## Acknowledging usage

If you use this library in your program and find it helpful, that's great! Any feedback would be much appreciated. If you publish results using this library, please consider citing the following paper detailing the implementation:

Coming soon.
