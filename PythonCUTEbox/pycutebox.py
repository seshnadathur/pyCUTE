import CUTEboxPython as cutebox
import numpy as np

"""
 Run CUTEbox from within Python using CUTEboxPython

 Allows to load a random_catalog and pass it to CUTE to avoid
 having to read this from file every time in case of multiple calls.

 Input:
    * Filename of CUTE parameterfile (standard CUTE format). If not provided
      then we assume the parameters have been set by first calling set_CUTEbox_parameters(...)
    * Catalog(s) in CUTE format
      If not provided then the catalogs is read from file inside CUTE
      The catalog can be created by calling readCatalog(filename)

 Output:
    * For corr_type = 1 [x, corr, paircounts] where [x] is [r]
      For corr_type = 2 [x, corr, paircounts] where [x] is [rt,rl]
      For corr_type = 3 [x, corr, paircounts] where [x] is [r,mu]
      See the CUTE readme or src/io.c for more info
      [corr] is the correlation function (radial, angular, monopole etc.).

 paircounts we get depends on the options we use:
      [ D1D1 ]               (if not use_randoms and not do_CFF)
      [ D1D2 ]               (if not use_randoms and     do_CFF)
      [ D1D1, D1R, RR ]      (if     use_randoms and not do_CFF)
      [ D1D2, D1R, D2R, RR ] (if     use_randoms and     do_CFF)

 We can also fetch paircount from [result] if wanted, see CUTEboxPython.i for
 availiable functions.

 Python has responsibillity for the memory of [catalog] (if
 they are passed to CUTE). Deallocation of catalogs are not handled automatically,
 so call cutebox.freeCatalog(random_catalog) once you are done with it.

"""
def runCUTEbox(paramfile = None, galaxy_catalog = None, galaxy_catalog2 = None, random_catalog = None, verbose = True):

  if(paramfile is not None):
    cutebox.read_run_params(paramfile)

  # Check for errors in parameters
  err = cutebox.verify_parameters()
  if(err > 0): return

  result = cutebox.make_empty_result_struct()
  cutebox.runCUTEbox(galaxy_catalog,galaxy_catalog2,random_catalog,result,verbose)

  do_CCF = cutebox.get_do_CCF()
  use_randoms = cutebox.get_use_randoms()

  # Fetch results
  if(cutebox.get_corr_type() == 1):
    nx   = result.get_nx()
    x    = np.array([result.get_x(i)    for i in range(nx)])
    corr = np.array([result.get_corr(i) for i in range(nx)])

    if(do_CCF == 1):
      D1D2   = np.array([result.get_D1D2(i) for i in range(nx)])
      D1R1   = np.array([result.get_D1R1(i) for i in range(nx)])
      D2R1   = np.array([result.get_D2R1(i) for i in range(nx)])
      R1R1   = np.array([result.get_R1R1(i) for i in range(nx)])
      if(use_randoms == 1):
        cutebox.free_result_struct(result)
        return x, corr, [D1D2, D1R1, D2R1, R1R1]
      else:
        cutebox.free_result_struct(result)
        return x, corr, [D1D2]
    else:
      D1D1   = np.array([result.get_D1D1(i) for i in range(nx)])
      D1R1   = np.array([result.get_D1R1(i) for i in range(nx)])
      R1R1   = np.array([result.get_R1R1(i) for i in range(nx)])
      if(use_randoms == 1):
        cutebox.free_result_struct(result)
        return x, corr, [D1D1, D1R1, R1R1]
      else:
        cutebox.free_result_struct(result)
        return x, corr, [D1D1]

  elif(cutebox.get_corr_type() == 2 or cutebox.get_corr_type() == 3):
    nx = result.get_nx()
    ny = result.get_ny()
    x = np.array([result.get_x(i) for i in range(nx)])
    y = np.array([result.get_y(i) for i in range(ny)])
    corr = np.zeros((nx,ny),dtype='float64')

    if(do_CCF == 1):
      D1D2 = np.zeros((nx,ny),dtype='float64')
      D1R1 = np.zeros((nx,ny),dtype='float64')
      D2R1 = np.zeros((nx,ny),dtype='float64')
      R1R1 = np.zeros((nx,ny),dtype='float64')
      for i in range(nx):
        for j in range(ny):
          corr[i][j] = result.get_corr(j + ny * i)
          D1D2[i][j] = result.get_D1D2(j + ny * i)
          D1R1[i][j] = result.get_D1R1(j + ny * i)
          D2R1[i][j] = result.get_D2R1(j + ny * i)
          R1R1[i][j] = result.get_R1R1(j + ny * i)
      if(use_randoms == 1):
        cutebox.free_result_struct(result)
        return [x, y], corr, [D1D2, D1R1, D2R1, R1R1]
      else:
        cutebox.free_result_struct(result)
        return [x, y], corr, [D1D2]
    else:
      D1D1 = np.zeros((nx,ny),dtype='float64')
      D1R1 = np.zeros((nx,ny),dtype='float64')
      R1R1 = np.zeros((nx,ny),dtype='float64')
      for i in range(nx):
        for j in range(ny):
          corr[i][j] = result.get_corr(j + ny * i)
          D1D1[i][j] = result.get_D1D1(j + ny * i)
          D1R1[i][j] = result.get_D1R1(j + ny * i)
          R1R1[i][j] = result.get_R1R1(j + ny * i)
      if(use_randoms == 1):
        cutebox.free_result_struct(result)
        return [x, y], corr, [D1D1, D1R1, R1R1]
      else:
        cutebox.free_result_struct(result)
        return [x, y], corr, [D1D1]

"""
 Set parameters in CUTE either by providing a parameterfile or
 by setting them directly
 If paramfile is not None then we read the parameterfile
 NB: if we read the parameterfile and there are critical errors then
 C calls exit() (standard in CUTE). If we set them one by one then
 errors will be shown, but no call to exit()
"""
def set_CUTEbox_parameters(
    paramfile=None,
    data_filename="",
    data_filename2="",
    random_filename="",
    use_randoms=0,
    reuse_randoms=0,
    num_lines="all",
    input_format=2,
    output_filename="",
    corr_type=1,
    use_pm=1,
    box_size=1.0,
    use_tree=0,
    max_tree_order=0,
    max_tree_nparts=0,
    do_CCF=0,
    n_grid_side=0):

  if(paramfile is not None):
    cutebox.read_run_params(paramfile)
    return

  # Strings
  cutebox.set_data_filename(data_filename)
  cutebox.set_data_filename2(data_filename2)
  cutebox.set_random_filename(random_filename);
  cutebox.set_output_filename(output_filename)
  cutebox.set_num_lines(num_lines)

  # Doubles
  cutebox.set_box_size(box_size)

  # Integers
  cutebox.set_corr_type(corr_type)
  cutebox.set_input_format(input_format)
  cutebox.set_use_randoms(use_randoms)
  cutebox.set_reuse_randoms(reuse_randoms)
  cutebox.set_use_pm(use_pm)
  cutebox.set_use_tree(use_tree)
  cutebox.set_max_tree_order(max_tree_order)
  cutebox.set_max_tree_nparts(max_tree_nparts)
  cutebox.set_do_CCF(do_CCF)
  cutebox.set_n_grid_side(n_grid_side)

  # Check if parameters are good
  err = cutebox.verify_parameters()

"""
 Free up memory of a C catalog
"""
def freeCatalog(catalog):
  cutebox.free_catalog(catalog)

"""
 Use CUTE to read a catalog and store it in C format
 If paramfile = None we assume the parameters have already been set in CUTE by set_CUTEbox_parameters
"""
def readCatalog(filename,input_format,box_size):
  cutebox.set_input_format(input_format)
  cutebox.set_box_size(box_size)
  return cutebox.read_Catalog(filename)

"""
  Create a CUTE tracer catalog in C format from numpy arrays of X, Y, Z positions
"""
def createCatalogFromNumpy(x, y, z):
  ok = True
  if (type(x)     is not np.ndarray): ok = False
  if (type(y)     is not np.ndarray): ok = False
  if (type(z)     is not np.ndarray): ok = False
  if( not ok ):
    print("Error: all input needs to be numpy double arrays")
    return None
  return cutebox.create_catalog_from_numpy(x, y, z)

"""
 Print CUTE parameters
"""
def print_param():
  cutebox.print_parameters()
