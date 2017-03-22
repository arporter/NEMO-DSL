# README #

This directory contains the 'PSyKAl step 1' version of the tra_adv
benchmark as extracted from the NEMO source code by CMCC. The step 1
version makes as few changes to the source code as possible whilst
being able to generate OpenMP and OpenACC parallel code.

## Building ##

The code uses the [dl_timer](https://bitbucket.org/apeg/dl_timer)
library for timing. Provided you followed the steps in the top-level
[README](../../../README.md#obtaining-the-code) then the source code
for this library will have been obtained when you cloned the
repository. The library is built automatically as part of the
compilation process.

The compiler and compiler flags to use must be set through the following
environment variables:

* F90      - the command with which to invoke the Fortran compiler
* F90FLAGS - flags to pass to the compiler, e.g. -g

e.g. to build with Gnu Fortran I use:

    export F90=gfortran
    export F90FLAGS=-O3

## Running ##

The problem size and number of iterations/time-steps to perform are also
set using environment variables:

* JPI - x-dimension of the domain
* JPJ - y-dimension of the domain
* JPK - number of vertical levels
* IT - number of iterations to perform

e.g.:

    export JPI=130
    export JPJ=130
    export JPK=31
    export IT=100
