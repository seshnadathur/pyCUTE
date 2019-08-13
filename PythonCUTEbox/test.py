import pycutebox
import numpy as np

'''
 Simple test of the pyCUTEbox wrapper
 test_type = 1 : set parameters by hand and run with do_CCF = 0 and no randoms
 test_type = 2 : set parameters from a paramfile and run with do_CFF = 0 and no randoms 
 test_type = 3 : supply CUTE with a galaxy catalog from python and run with do_CFF = 0 and no randoms
 test_type = 4 : supply all catalogs from python and run with do_CFF = 0 and with randoms
 test_type = 5 : supply all catalogs from python and run with do_CFF = 1 and with randoms
'''
def run_a_test(test_type):
  
  if(test_type == 1):
    
    #############################################################################
    
    # Set CUTE parameters directly
    pycutebox.set_CUTEbox_parameters(
        data_filename="test/data.dat",
        data_filename2="notdefault",
        random_filename="notdefault",
        use_randoms=0,
        reuse_randoms=0,
        num_lines="all",
        input_format=0,
        output_filename="test/corr.dat",
        corr_type=1,
        box_size=500.0,
        use_tree=0,
        max_tree_order=6,
        max_tree_nparts=100,
        do_CCF=0,
        n_grid_side=256,
        use_pm=0)
    
    # Run CUTEbox with parameters already set
    x, corr, paircounts = pycutebox.runCUTEbox()
    
    print x
    print corr
    print paircounts
  
    #############################################################################
  
  elif(test_type == 2):
    
    #############################################################################
    
    # Run CUTEbox directly from a parameterfile
    x, corr, paircounts = pycutebox.runCUTEbox(paramfile = "test/param_test.txt")
    
    print x
    print corr
    print paircounts
    
    #############################################################################
  
  elif(test_type == 3):
    
    #############################################################################
    
    # Read a galaxy catalog from file and supply it to CUTEbox
    galaxy_catalog = pycutebox.readCatalog(filename = "test/data.dat", input_format = 0, box_size = 500.0)
    
    # Set CUTE parameters directly
    pycutebox.set_CUTEbox_parameters(
        data_filename="not_relevant_as_we_supply_it",
        data_filename2="notdefault",
        random_filename="notdefault",
        use_randoms=0,
        reuse_randoms=0,
        num_lines="all",
        input_format=0,
        output_filename="test/corr.dat",
        corr_type=1,
        box_size=500.0,
        use_tree=0,
        max_tree_order=6,
        max_tree_nparts=100,
        do_CCF=0,
        n_grid_side=256,
        use_pm=0)
    
    # Run CUTE box with preread data
    x, corr, paircounts = pycutebox.runCUTEbox(galaxy_catalog = galaxy_catalog)
    
    print x
    print corr
    print paircounts
  
    # Free up data
    pycutebox.freeCatalog(galaxy_catalog)
    
    #############################################################################
  
  elif(test_type == 4):
    
    #############################################################################
    
    # Read all the data from file and supply it to CUTEbox
    galaxy_catalog  = pycutebox.readCatalog(filename = "test/data.dat",  input_format = 0, box_size = 500.0)
    galaxy_catalog2 = pycutebox.readCatalog(filename = "test/data2.dat", input_format = 0, box_size = 500.0)
    random_catalog  = pycutebox.readCatalog(filename = "test/rand.dat",  input_format = 0, box_size = 500.0)
    
    # Set CUTE parameters directly
    pycutebox.set_CUTEbox_parameters(
        data_filename="not_relevant_as_we_supply_it",
        data_filename2="not_relevant_as_we_supply_it",
        random_filename="not_relevant_as_we_supply_it",
        use_randoms=1,
        reuse_randoms=0,
        num_lines="all",
        input_format=0,
        output_filename="test/corr.dat",
        corr_type=1,
        box_size=500.0,
        use_tree=0,
        max_tree_order=6,
        max_tree_nparts=100,
        do_CCF=0,
        n_grid_side=256,
        use_pm=0)
    
    # Run CUTE box with randoms
    x, corr, paircounts = pycutebox.runCUTEbox(galaxy_catalog = galaxy_catalog, galaxy_catalog2 = galaxy_catalog2, random_catalog = random_catalog)
    
    print x
    print corr
    print paircounts
    
    # Free up data
    pycutebox.freeCatalog(galaxy_catalog)
    pycutebox.freeCatalog(galaxy_catalog2)
    pycutebox.freeCatalog(random_catalog)
    
    #############################################################################
  
  elif(test_type == 5):
    
    #############################################################################
    
    # Read all the data from file and supply it to CUTEbox
    galaxy_catalog  = pycutebox.readCatalog(filename = "test/data.dat",  input_format = 0, box_size = 500.0)
    galaxy_catalog2 = pycutebox.readCatalog(filename = "test/data2.dat", input_format = 0, box_size = 500.0)
    random_catalog  = pycutebox.readCatalog(filename = "test/rand.dat",  input_format = 0, box_size = 500.0)
    
    # Set CUTE parameters directly
    pycutebox.set_CUTEbox_parameters(
        data_filename="not_relevant_as_we_supply_it",
        data_filename2="not_relevant_as_we_supply_it",
        random_filename="not_relevant_as_we_supply_it",
        use_randoms=1,
        reuse_randoms=0,
        num_lines="all",
        input_format=0,
        output_filename="test/corr.dat",
        corr_type=1,
        box_size=500.0,
        use_tree=0,
        max_tree_order=6,
        max_tree_nparts=100,
        do_CCF=1,
        n_grid_side=256,
        use_pm=0)
    
    # Run CUTE box with randoms and CCF
    x, corr, paircounts = pycutebox.runCUTEbox(galaxy_catalog = galaxy_catalog, galaxy_catalog2 = galaxy_catalog2, random_catalog = random_catalog)
  
    print x
    print corr
    print paircounts
    
    # Free up data
    pycutebox.freeCatalog(galaxy_catalog)
    pycutebox.freeCatalog(galaxy_catalog2)
    pycutebox.freeCatalog(random_catalog)
  
    #############################################################################


#############################################################################
# Test different ways of using the script
#############################################################################

run_a_test(1)

run_a_test(2)

run_a_test(3)

run_a_test(4)

run_a_test(5)


