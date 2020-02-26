%module CUTEboxPython
%{
  #define SWIG_FILE_WITH_INIT
  #include "src/define.h"
  #include "src/common.h"
  extern int runCUTEbox(Catalog *galaxy_catalog, Catalog *galaxy_catalog2, Catalog *random_catalog, Result *result, int verbose);
  extern Catalog *read_Catalog(char *fname);
  extern void free_catalog(Catalog *cat);
  extern void read_run_params(char *paramfile);
  extern Result *make_empty_result_struct();
  extern void free_result_struct(Result *res);
  extern int get_corr_type();
  extern int get_do_CCF();
  extern int get_use_randoms();

  extern Catalog *create_catalog_from_numpy(int nx, double *x, int ny, double *y, int nz, double *z);

  extern int verify_parameters();
  extern void print_parameters();

  extern void set_data_filename(char *s);
  extern void set_data_filename2(char *s);
  extern void set_random_filename(char *s);
  extern void set_input_format(int i);
  extern void set_output_filename(char *s);
  extern void set_corr_type(int i);
  extern void set_use_pm(int i);
  extern void set_reuse_randoms(int i);
  extern void set_num_lines(char *s);
  extern void set_box_size(double x);
  extern void set_use_tree(int i);
  extern void set_max_tree_order(int i);
  extern void set_max_tree_nparts(int i);
  extern void set_use_randoms(int i);
  extern void set_do_CCF(int i);
  extern void set_n_grid_side(int i);

  struct Catalog{
  #ifdef _LONGIDS
    long np;
  #else
    int np;
  #endif
    double *pos;
  };

  struct Result {
    int nx, ny, nz;
    double *x, *y, *z, *corr,
         *D1D1, *D1D2, *D1R1, *D1R2,
         *D2D2, *D2R1, *D2R2,
         *R1R1, *R1R2,
         *R2R2;
  };
%}

%include "numpy.i"
%init %{
import_array();
%}
%apply (int DIM1, double* INPLACE_ARRAY1) {(int n0, double *a0)};
%apply (int DIM1, double* IN_ARRAY1) {(int nx, double *x), (int ny, double *y), (int nz, double *z)};

#include "src/define.h"
#include "src/common.h"
extern int runCUTEbox(Catalog *galaxy_catalog, Catalog *galaxy_catalog2, Catalog *random_catalog, Result *result, int verbose);
extern Catalog *read_Catalog(char *fname);
extern void free_catalog(Catalog *cat);
extern void read_run_params(char *paramfile);
extern Result *make_empty_result_struct();
extern void free_result_struct(Result *res);
extern int get_corr_type();
extern int get_do_CCF();
extern int get_use_randoms();

extern Catalog *create_catalog_from_numpy(int nx, double *x, int ny, double *y, int nz, double *z);

extern int verify_parameters();
extern void print_parameters();

void set_data_filename(char *s);
void set_data_filename2(char *s);
void set_random_filename(char *s);
void set_input_format(int i);
void set_output_filename(char *s);
void set_corr_type(int i);
void set_use_pm(int i);
void set_reuse_randoms(int i);
void set_num_lines(char *s);
void set_box_size(double x);
void set_use_tree(int i);
void set_max_tree_order(int i);
void set_max_tree_nparts(int i);
void set_use_randoms(int i);
void set_do_CCF(int i);
void set_n_grid_side(int i);

struct Catalog{
#ifdef _LONGIDS
  long np;
#else
  int np;
#endif
  double *pos;
};

%extend Catalog{
  int get_np(){
    return $self->np;
  }
  double get_pos(int i) {
    return $self->pos[i];
  }
}

%typemap(newfree) Catalog * {
  free_catalog($1);
}

struct Result {
  int nx, ny, nz;
  double *x, *y, *z, *corr,
         *D1D1, *D1D2, *D1R1, *D1R2,
         *D2D2, *D2R1, *D2R2,
         *R1R1, *R1R2,
         *R2R2;
};

%extend Result{
  int get_nx(){
    return $self->nx;
  }
  int get_ny(){
    return $self->ny;
  }
  int get_nz(){
    return $self->nz;
  }
  double get_x(int i) {
    if(self->nx > 0)
      return $self->x[i];
    return 0.0;
  }
  double get_y(int i) {
    if(self->ny > 0)
      return $self->y[i];
    return 0.0;
  }
  double get_z(int i) {
    if(self->nz > 0)
      return $self->z[i];
    return 0.0;
  }
  double get_corr(int i) {
    return $self->corr[i];
  }
  double get_D1D1(int i) {
    return $self->D1D1[i];
  }
  double get_D1D2(int i) {
    return $self->D1D2[i];
  }
  double get_D1R1(int i) {
    return $self->D1R1[i];
  }
  double get_D1R2(int i) {
    return $self->D1R2[i];
  }
  double get_D2D2(int i) {
    return $self->D2D2[i];
  }
  double get_D2R1(int i) {
    return $self->D2R1[i];
  }
  double get_D2R2(int i) {
    return $self->D2R2[i];
  }
  double get_R1R1(int i) {
    return $self->R1R1[i];
  }
  double get_R1R2(int i) {
    return $self->R1R2[i];
  }
  double get_R2R2(int i) {
    return $self->R2R2[i];
  }
  ~Result(){
    free($self->x);
    free($self->y);
    free($self->z);
    free($self->corr);
    free($self->D1D1);
    free($self->D1D2);
    free($self->D1R1);
    free($self->D1R2);
    free($self->D2D2);
    free($self->D2R1);
    free($self->D2R2);
    free($self->R1R1);
    free($self->R1R2);
    free($self->R2R2);
  }
}

%typemap(newobject) make_empty_result_struct;
%typemap(newfree) Result * {
  free_result_struct($1);
}
