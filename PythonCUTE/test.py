import pycute 
import numpy as np

'''
# Simple test. Run from a parameterfile
paramfile="test/param.ini"
x, corr, DD, DR, RR = pycute.runCUTE(paramfile = paramfile)
print "r    = ", x
print "corr = ", corr
print "DD   = ", DD
print "DR   = ", DR
print "RR   = ", RR

# Read random catalog in python and call CUTE with this choice
paramfile="fiducial_N.ini"
random_filename = "/mnt/lustre/nadathur/BOSS_DR12_voidRSD/tracer_files_for_CUTE/zobov_voids/reconVoids-Rcut-randoms-DR12-CMASS-N-x50.txt"
random_catalog   = pycute.readCatalog(paramfile,  random_filename)
x, y, corr, D1D2, D1R2, D2R1, R1R2 = pycute.runCUTE(paramfile = paramfile, random_catalog = random_catalog)
print x
print y
print corr
print D1D2
print D1R2
print D2R1
print R1R2

# Print CUTE parameters
pycute.print_param()

# Set CUTE parameters from file
pycute.set_CUTE_parameters(paramfile="test/param.ini")
'''

# Set CUTE parameters directly
pycute.set_CUTE_parameters(
    data_filename="notdefault",
    data_filename2="notdefault",
    random_filename="notdefault",
    random_filename2="notdefault",
    reuse_randoms=0,
    num_lines="all",
    input_format=2,
    mask_filename="notdefault",
    z_dist_filename="notdefault",
    output_filename="notdefault",
    corr_estimator="LS",
    corr_type="monopole",
    np_rand_fact=10,
    omega_M=0.3,
    omega_L=0.7,
    w=-1.0,
    radial_aperture=1.0,
    dim1_max=120,
    dim2_max=1.,
    dim3_max=0.73,
    dim3_min=0.4,
    dim1_nbin=30,
    dim2_nbin=80,
    dim3_nbin=1,
    log_bin=0,
    n_logint=10,
    use_pm=1,
    n_pix_sph=2048)

# Finalize (in case of MPI)
pycute.Finalize()
