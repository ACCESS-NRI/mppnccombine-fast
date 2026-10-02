Installation
============

For installation instructions see the README at https://github.com/coecms/mppnccombine-fast

Usage
=====

``mppnccombine-fast`` is a MPI program, and requires at least two MPI ranks to run::

    mpirun -n 2 mppnccombine-fast --output output.nc input.nc.000 input.nc.001

Variables in the input files whos dimensions have a ``domain_distribution``
attribute will be collated. All other dimensions, variables and attributes will
be copied from the first input file.

The ``domain_distribution`` values are expected to be in the format provided by
the MOM model - an array of 4 integer values using 1-based array indices:

 * First index of this dimension in the full dataset
 * Last index of this dimension in the full dataset
 * First index of this dimension in this file's data
 * Last index of this dimension in this file's data

A ``domain_distribution`` of ``[1, 10, 5, 10]`` states that the full dimension
has a length of 10, and this file contains the 5 values starting at offset
``5``.

Globbing inputs
---------------

Input files may be listed either individually or as an escaped shell glob (both
to reduce the ``history`` attribute in the output file as well as to avoid
issues when there are thousands of input files)::

    mpirun -n 2 mppnccombine-fast --output output.nc input.nc.\*

Missing tiles
-------------

Some MOM runs use a mask table to skip tiles that are entirely land, so
no file is written for those tiles. Regions of collated variables that
no input file covers are filled with the variable's ``_FillValue``
where one is declared. Where none is declared (as in MOM6 restart
files), they are filled with zero, matching ``mppnccombine``.

Changing compression settings
-----------------------------

Chunk size and compression settings will by default come from the first input
file, though they can be overridden using flags. Note that the optimised
copying routines can only be used when the compression settings of an input
file matches those of the output, and when the input file's data chunks align
with the chunks in the output file (e.g. if a variable in the  output file has
chunk sizes ``[10, 15, 30]`` then the input file's offset in the full dataset
must be ``[m*10, n*15, o*30]`` where ``m``, ``n`` and ``o`` are integers).

If only some of the chunks in the input file align with the output these chunks
will use the fast path (so partial chunks on the edges of the dataset are fine).
