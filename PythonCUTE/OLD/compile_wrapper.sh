#!/bin/bash

# The include path to the python libraries
PYTHON_INC="-I/opt/apps/python/2.7.8/gcc-4.4.7/include/python2.7"

# Create wrapper
swig -python swigfile.i

# Compile wrapper
gcc -fPIC -c swigfile_wrap.c $PYTHON_INC

# Make CUTE
make clean; make
