# mppnccombine-fast
[![DOI](https://zenodo.org/badge/131938571.svg)](https://zenodo.org/badge/latestdoi/131938571)

An accelerated version of the `mppnccombine` post-processing tool for MOM

Uses HDF5's raw IO functions to speed up collating large datasets - a 0.1
degree model goes from taking 4 hours to collate a compressed variable with
mppnccombine, to 6 minutes with mppnccombine-fast running with 16 processes

Documentation at https://mppnccombine-fast.readthedocs.io

## Build

`mppnccombine-fast` requires HDF5 version 1.10.2 or above, as well as NetCDF 4,
a C compiler and a MPI library.

Cmake is used for building. A Makefile is also provided for ease of use:

    make # Release build

    make BUILD_TYPE=Debug # Debug build

    make check # Run tests

    make PREFIX=/apps/mppncc-fast install # Install mppnccombine-fast to $PREFIX

    make doc # Build documentation in doc/_build/html

On Raijin the Makefile loads all required modules

The environment variables `$OPENMPI_ROOT`, `$HDF5_ROOT` and `$NETCDF_ROOT` may
be used to help locate libraries at other sites

## Use

Use like

    mpirun -n 2 mppnccombine-fast --output out.nc input.nc.0000 input.nc.0001 input.nc.0002

Files will be collated along all axes with a `domain_distribution` attribute

At least 2 MPI ranks need to be used (rank 0 writes the output file, other
ranks read). More can be used - input files will be balanced between the MPI
ranks.

A full list of options can be found by running

    mppnccombine-fast --help

