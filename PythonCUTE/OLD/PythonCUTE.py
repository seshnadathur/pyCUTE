import CUTEPythonWrapper as cute

#===============================================
# Initialize [if the paramfile changes between 
# runs we need to call this function again]
#===============================================
paramfile = "test/param.ini"
cute.read_run_params(paramfile)

#===============================================
# Read in the random catalog
#===============================================
random_filename = "test/random.dat"
random_catalog = cute.read_random_catalog(random_filename)

#===============================================
# Run CUTE
#===============================================
cute.runCUTE(random_catalog)

exit(1)

#===============================================
# We can run it again...
#===============================================
cute.runCUTE(random_catalog)

#===============================================
# Free up memory for the random catalog
#===============================================
cute.free_random_catalog()

#===============================================
# Run CUTE without having an external random catalog
#===============================================
cute.runCUTE(None)

