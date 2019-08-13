import CUTEPython as cute
import numpy as np

###############################################################
# Definition of helper functions:
###############################################################
"""
 Run CUTE from within Python using CUTEPython
 
 Allows to load a random_catalog and pass it to CUTE to avoid 
 having to read this from file every time in case of multiple calls.
 
 Input:
    * Filename of CUTE parameterfile (standard CUTE format). If not provided
      then we assume the parameters have been set by first calling set_CUTE_parameters(...)
    * Random catalog in CUTE format
      If not provided then the random catalog
      is read from file inside CUTE
      The random catalog can be created by calling 
      cute.read_random_catalog(paramfile, random_filename)

 Output:
    * [x], [corr]. The [x] depends on the [corr_type] option
      in CUTE. [Dz] for corr_type = radial, [theta] for corr_type = angular
      and [r] for corr_type = monopole. See the CUTE readme for more info. 
      If corr_type is different then we results are just written to file and not
      returned in the arrays below (too lazy to implement this).
      [corr] is the correlation function (radial, angular or monopole).

 We can also fetch paircount from [result] if wanted, see CUTEPython.i for 
 availiable functions.

 Python has responsibillity for the memory of [result] and [random_catalog] (if 
 random_catalog != None). Deallocation should be handled automatically, but not tested
 so to be sure of no memory leaks we can always call cute.free_result_struct(result) 
 and cute.free_Catalog(random_catalog) 

 MPI support is implemented, but not well tested. To run with MPI
 run as OMP_NUM_THREADS=1 mpirun -np N python2.7 script.py and remember 
 to call cute.finalize_mpi() in the end of the script.

 TODO: 
   * Should implement returning results for corr_type > 2.

   * Should allow for also passing galaxy catalog from Python to CUTE.
     
   * Should allow to also pass a second set of randoms/galaxies in case of use_two_catalogs = True

"""
def runCUTE(paramfile = None, random_catalog = None):
  if(paramfile != None): 
    cute.read_run_params(paramfile)
  result = cute.make_empty_result_struct()
  cute.runCUTE(random_catalog,result)
  if(cute.get_corr_type() > 2): 
    return None, None
  else:
    nbins = result.get_nbins()
    x     = np.array([result.get_x(i)    for i in range(nbins)])
    corr  = np.array([result.get_corr(i) for i in range(nbins)])
    cute.free_result_struct(result)
    return x, corr

"""
  Set parameters in CUTE directly without having to do the detour of writing an
  external parameterfile
"""
def set_CUTE_parameters(
    data_filename="file_none",
    data_filename_2="file_none",
    random_filename="file_none",
    random_filename_2="file_none",
    RR_filename="file_none",
    output_filename="file_none",
    corr_type="radial/angular/monopole",
    omega_M=-10.0,
    omega_L=-10.0,
    w=-10.0,
    radial_aperture=0.0,
    dim1_max=-1,
    dim1_min_logbin=-1,
    dim2_max=-1,
    dim3_max=-1,
    dim3_min=-1,
    input_format=-1,
    dim1_nbin=-1,
    dim2_nbin=-1,
    dim3_nbin=-1,
    log_bin=-1,
    use_pm=-1,
    n_pix_sph=-1):

  cute.initialize_binner();
  
  # Strings
  cute.set_data_filename(data_filename)
  cute.set_data_filename_2(data_filename_2)
  cute.set_random_filename(random_filename);
  cute.set_random_filename_2(random_filename_2);
  cute.set_RR_filename(RR_filename)
  cute.set_output_filename(output_filename)
  cute.set_corr_type(corr_type)

  # Doubles
  cute.set_omega_M(omega_M)
  cute.set_omega_L(omega_L)
  cute.set_w(w)
  cute.set_radial_aperture(radial_aperture)
  cute.set_dim1_max(dim1_max)
  cute.set_dim1_min_logbin(dim1_min_logbin)
  cute.set_dim2_max(dim2_max)
  cute.set_dim3_max(dim3_max)
  cute.set_dim3_min(dim3_min)

  # Integers
  cute.set_input_format(input_format)
  cute.set_dim1_nbin(dim1_nbin)
  cute.set_dim2_nbin(dim2_nbin)
  cute.set_dim3_nbin(dim3_nbin)
  cute.set_log_bin(log_bin)
  cute.set_use_pm(use_pm)
  cute.set_n_pix_sph(n_pix_sph)

  # Check if parameters are good 
  err = cute.verify_parameters()

###############################################################
# Examples of use below:
###############################################################

#==============================================================
# Set parameters directly and run the code
#==============================================================
set_CUTE_parameters(
      data_filename   = "test/shell.dat",
      random_filename = "test/random.dat",
      input_format    = 0,
      output_filename = "test/corr_full_pm.dat",
      RR_filename     = "test/rr.dat",
      corr_type       = "angular",
      omega_M         = 0.315,
      omega_L         = 0.685,
      w               = -1,
      log_bin         = 0,
      dim1_max        = 10.,
      dim1_min_logbin = 0.01,
      dim1_nbin       = 30,
      dim2_max        = 0.1,
      dim2_nbin       = 30,
      dim3_min        = 0.4,
      dim3_max        = 0.7,
      dim3_nbin       = 1,
      radial_aperture = 1.,
      use_pm          = 0,
      n_pix_sph       = 2048)
x, corr = runCUTE()

#====================================================
# Read in and store a random catalog and run the code
# with this random catalog and the rest of the params
# from an external parameterfile
#====================================================
paramfile       = "test/param.ini"
random_filename = "test/random.dat"
random_catalog  = cute.read_random_catalog(paramfile, random_filename)
x, corr = runCUTE(paramfile, random_catalog)

#====================================================
# Run CUTE from an external parameterfile
#====================================================
x, corr = runCUTE(paramfile)

#====================================================
# If MPI then finalize this within CUTE and free up memory
#====================================================
cute.finalize_mpi()
cute.free_Catalog(random_catalog)

