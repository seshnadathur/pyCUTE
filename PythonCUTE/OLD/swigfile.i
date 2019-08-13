%module CUTEPythonWrapper
%{
  #include "src/define.h"
  #include "src/common.h"
  extern int runCUTE(Catalog *random_catalog);
  extern Catalog *read_random_catalog(char *fname);
  extern void free_random_catalog();
  extern void read_run_params(char *paramfile);
%}

#include "src/define.h"
#include "src/common.h"
extern int runCUTE(Catalog *random_catalog);
extern Catalog *read_random_catalog(char *fname);
extern void free_random_catalog();
extern void read_run_params(char *paramfile);

