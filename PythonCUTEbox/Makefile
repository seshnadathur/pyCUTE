####################################################
###          User definable stuff
#radial binning
PYTHON_LIBRARY = yes
NB_R = 50
R_MAX = 120
#Logarithmic binning
N_LOGINT = 9
#Angular bins (only used for 3d r-mu correlations)
NB_mu = 100
# Python paths
PYTHONINC  = -I/opt/apps/pkgs/anaconda3/2019.03/intel64/include/python3.7m
PYTHONINC += -I/opt/apps/pkgs/anaconda3/2019.03/intel64/lib/python3.7/site-packages/numpy/core/include$
PYTHONLIB  = -L/opt/apps/pkgs/anaconda3/2019.03/intel64/lib -lpython3.7m
#Behavior options
DEFINEOPTIONS = -D_VERBOSE -D_LONGIDS
DEFINEOPTIONS += -D_HAVE_OMP #Comment this out if you don't have the OpenMP headers
#DEFINEOPTIONS += -D_DEBUG
#DEFINEOPTIONS += -D_LOGBIN
### End of user-definable stuff
####################################################

#Transform to inverse scales
I_R_MAX = $(shell echo "scale=5;1./$(R_MAX)" | bc)
LOG_R_MAX = $(shell echo "scale=9;l($(R_MAX))/l(10)" | bc -l)

# DEFINES for the OpenMP version
DEFINEFLAGSCPU = $(DEFINEOPTIONS)
DEFINEFLAGSCPU += -DNB_R=$(NB_R) -DI_R_MAX=$(I_R_MAX) -DLOG_R_MAX=$(LOG_R_MAX)
DEFINEFLAGSCPU += -DN_LOGINT=$(N_LOGINT)

# COMPILER AND OPTIONS
COMPCPU = gcc
OPTCPU = -Wall -O3 -fopenmp $(DEFINEFLAGSCPU)
ifeq ($(strip $(PYTHON_LIBRARY)),yes)
DEFINEFLAGSCPU += -D_CUTE_AS_PYTHON_MODULE
OPTCPU += -fPIC
endif

#INCLUDES AND LIBRARIES
INCLUDECOM = -I./src
LIBCPU = -lm

#.o FILES
#CUTE
DEF = src/define.o
COM = src/common.o
PYCUTE = src/pycute.o
PM = src/pm.o
TREE = src/tree.o
IO = src/io.o
NEIGH = src/neighbors.o
CORR = src/correlator.o
MAIN = src/main.c
OFILES = $(DEF) $(COM) $(PYCUTE) $(PM) $(TREE) $(IO) $(NEIGH) $(CORR) $(MAIN)

#FINAL GOAL
EXE = CUTE_box

#RULES
default : $(EXE)
#RULE TO MAKE .o's FROM .c's
$(DEF) : src/define.c
	$(COMPCPU) $(OPTCPU) -c $< -o $@ $(INCLUDECOM)
$(COM) : src/common.c
	$(COMPCPU) $(OPTCPU) -c $< -o $@ $(INCLUDECOM)
$(PYCUTE) : src/pycute.c
	$(COMPCPU) $(OPTCPU) -c $< -o $@ $(INCLUDECOM)
$(CORR) : src/correlator.c
	$(COMPCPU) $(OPTCPU) -c $< -o $@ $(INCLUDECOM)
$(TREE) : src/tree.c
	$(COMPCPU) $(OPTCPU) -c $< -o $@ $(INCLUDECOM)
$(IO) : src/io.c
	$(COMPCPU) $(OPTCPU) -c $< -o $@ $(INCLUDECOM)
$(NEIGH) : src/neighbors.c
	$(COMPCPU) $(OPTCPU) -c $< -o $@ $(INCLUDECOM)
$(PM) : src/pm.c
	$(COMPCPU) $(OPTCPU) -c $< -o $@ $(INCLUDECOM)

#RULES TO MAKE THE FINAL EXECUTABLES
$(EXE) : $(OFILES)
ifeq ($(strip $(PYTHON_LIBRARY)),yes)
	swig -python $(DEFINEFLAGSCPU) CUTEboxPython.i
	$(COMPCPU) $(OPTCPU) -c CUTEboxPython_wrap.c $(PYTHONINC)
	$(COMPCPU) $(OPTCPU) -shared $(OFILES) CUTEboxPython_wrap.o -o _CUTEboxPython.so $(INCLUDECOM) $(LIBCPU) $(PYTHONLIB)
else
	$(COMPCPU) $(OPTCPU) $(OFILES) -o $(EXE) $(INCLUDECOM) $(LIBCPU)
endif

#CLEANING RULES
clean :
	rm -f ./src/*.o CUTEboxPython.py CUTEboxPython.pyc CUTEboxPython_wrap* _CUTEboxPython.so pycutebox.pyc

cleaner :
	rm -f ./src/*.o ./src/*~ *~ CUTE_box
