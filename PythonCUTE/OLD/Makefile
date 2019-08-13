####################################################
###          User definable stuff
####################################################

#Behavior options
USE_MPI = no                      #Use MPI parallelization? Set to "yes" or "no"
DEFINEOPTIONS = -D_VERBOSE        #Comment out for more explicit output
DEFINEOPTIONS += -D_HAVE_OMP      #Comment this out if you don't have the OpenMP headers
DEFINEOPTIONS += -D_WITH_WEIGHTS  #Comment this out if you want weighed galaxies
#DEFINEOPTIONS += -D_DEBUG        #Only used for debugging
#DEFINEOPTIONS += -D_TRUE_ACOS    #Comment this out if you don't want to use an approximate arccosine

# For linking to python
PYTHONLIB = -L/opt/apps/python/2.7.8/gcc-4.4.7/lib -lpython2.7

#GSL options
GSL_INC = /opt/local/include
GSL_LIB = /opt/local/lib
LGSL = -L$(GSL_LIB) -lgsl -lgslcblas

# DEFINES for the OpenMP version
DEFINEFLAGSCPU = $(DEFINEOPTIONS)

# COMPILER AND OPTIONS
ifeq ($(strip $(USE_MPI)),yes)
DEFINEFLAGSCPU += -D_HAVE_MPI
COMPCPU = mpicc
else
COMPCPU = gcc
endif

OPTCPU = -Wall -O3 -fopenmp $(DEFINEFLAGSCPU)

#INCLUDES AND LIBRARIES
INCLUDECOM = -I./src -I$(GSL_INC)
LIBCPU = $(LGSL) -lm

#.o FILES
#CUTE
DEF = src/define.o
COM = src/common.o
COSMO = src/cosmo.o
CORR = src/correlator.o
BOX2D = src/boxes2D.o
BOX3D = src/boxes3D.o
IO = src/io.o
MAIN = src/main.c
OFILES = $(DEF) $(COM) $(COSMO) $(CORR) $(BOX2D) $(BOX3D) $(IO) $(MAIN)

#FINAL GOAL
EXE = CUTE

#RULES
default : $(EXE) 
#RULE TO MAKE .o's FROM .c's
$(DEF) : src/define.c
	$(COMPCPU) $(OPTCPU) -fPIC -c $< -o $@ $(INCLUDECOM)
$(COM) : src/common.c
	$(COMPCPU) $(OPTCPU) -fPIC -c $< -o $@ $(INCLUDECOM)
$(COSMO) : src/cosmo.c
	$(COMPCPU) $(OPTCPU) -fPIC -c $< -o $@ $(INCLUDECOM)
$(CORR) : src/correlator.c
	$(COMPCPU) $(OPTCPU) -fPIC -c $< -o $@ $(INCLUDECOM)
$(IO) : src/io.c
	$(COMPCPU) $(OPTCPU) -fPIC -c $< -o $@ $(INCLUDECOM)
$(BOX2D) : src/boxes2D.c
	$(COMPCPU) $(OPTCPU) -fPIC -c $< -o $@ $(INCLUDECOM)
$(BOX3D) : src/boxes3D.c
	$(COMPCPU) $(OPTCPU) -fPIC -c $< -o $@ $(INCLUDECOM)

#RULES TO MAKE THE FINAL EXECUTABLES
$(EXE) : $(OFILES)
	$(COMPCPU) -fPIC -shared $(OPTCPU) $(OFILES) swigfile_wrap.o -o _CUTEPythonWrapper.so $(INCLUDECOM) $(LIBCPU) $(PYTHONLIB) 

#CLEANING RULES
clean :
	rm -f ./src/*.o

cleaner :
	rm -f ./src/*.o ./src/*~ *~ $(EXE) 
