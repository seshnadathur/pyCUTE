# pyCUTE
A Python-wrapped and modified version of the CUTE correlation function code 

The original CUTE (Correlation Utilities and Two point Estimation) code was written by David Alonso and is available at 
https://github.com/damonge/CUTE. This version is derived from that original, but with added functionalities and with a 
Python wrapper.

The additional functionalities over the original CUTE are:
- option to compute cross-correlation function for two different point populations
- for periodic boxes, the option to compute xi(sigma, pi) and xi(r, mu) in addition to the monopole; these 
functionalities are also present for cross-correlations in both the box version and the sky survey version of the code

Note: The computation of cross-correlations requires the input of two different data catalogues, D1 and D2, and (in the 
sky survey case), two corresponding random catalogues R1 and R2.

### Requirements
You will need SWIG installed (which in turn requires GSL and PCRE), see PythonCUTE/PYTHON_README. You need to change the
PYTHONINC and PYTHONLIB paths in the Makefile to match your system (for PythonCUTE, also change GSL_INC and GSL_LIB). 

### Acknowledgments
In addition to everyone who contributed to the original development of CUTE and is mentioned in that GitHub download, 
contributions to this version of pyCUTE were made by:
- Seshadri Nadathur
- Hans Winther
