%module CUTEPython
%{
  #include "src/define.h"
  #include "src/common.h"
  extern int runCUTE(Catalog *random_catalog, Result* corr_func);
  extern Catalog *read_random_catalog(char *paramfile, char *fname);
  extern void free_Catalog(Catalog *cat);
  extern void read_run_params(char *paramfile);
  extern Result *make_empty_result_struct();
  extern void free_result_struct(Result *res);
  extern int get_corr_type();
  extern void finalize_mpi();

  extern void initialize_binner();
  extern int verify_parameters();
  extern void set_data_filename(char *s);
  extern void set_data_filename_2(char *s);
  extern void set_random_filename(char *s);
  extern void set_random_filename_2(char *s);
  extern void set_RR_filename(char *s);
  extern void set_input_format(int i);
  extern void set_output_filename(char *s);
  extern void set_corr_type(char *s);
  extern void set_omega_M(double x);
  extern void set_omega_L(double x);
  extern void set_w(double x);
  extern void set_radial_aperture(double x);
  extern void set_dim1_max(double x);
  extern void set_dim1_min_logbin(double x);
  extern void set_dim2_max(double x);
  extern void set_dim3_max(double x);
  extern void set_dim3_min(double x);
  extern void set_dim1_nbin(int i);
  extern void set_dim2_nbin(int i);
  extern void set_dim3_nbin(int i);
  extern void set_log_bin(int i);
  extern void set_use_pm(int i);
  extern void set_n_pix_sph(int i);

  struct Catalog{
    int np;
    double *red,*cth,*phi;
  #ifdef _WITH_WEIGHTS
    double *weight;
  #endif
    np_t sum_w, sum_w2;
  };

  struct Result {
    int nbins;
    double *x, *corr, *DD, *DR, *RD, *RR;
  };
%}

#include "src/define.h"
#include "src/common.h"
extern int runCUTE(Catalog *random_catalog, Result* corr_func);
extern Catalog *read_random_catalog(char *paramfile, char *fname);
extern void free_Catalog(Catalog *cat);
extern void read_run_params(char *paramfile);
extern Result *make_empty_result_struct();
extern void free_result_struct(Result *res);
extern int get_corr_type();
extern void finalize_mpi();

extern void initialize_binner();
extern int verify_parameters();
extern void set_data_filename(char *s);
extern void set_data_filename_2(char *s);
extern void set_random_filename(char *s);
extern void set_random_filename_2(char *s);
extern void set_RR_filename(char *s);
extern void set_input_format(int i);
extern void set_output_filename(char *s);
extern void set_corr_type(char *s);
extern void set_omega_M(double x);
extern void set_omega_L(double x);
extern void set_w(double x);
extern void set_radial_aperture(double x);
extern void set_dim1_max(double x);
extern void set_dim1_min_logbin(double x);
extern void set_dim2_max(double x);
extern void set_dim3_max(double x);
extern void set_dim3_min(double x);
extern void set_dim1_nbin(int i);
extern void set_dim2_nbin(int i);
extern void set_dim3_nbin(int i);
extern void set_log_bin(int i);
extern void set_use_pm(int i);
extern void set_n_pix_sph(int i);

struct Catalog{
  int np;
  double *red,*cth,*phi;
#ifdef _WITH_WEIGHTS
  double *weight;
#endif
  np_t sum_w, sum_w2;
};

%extend Catalog{
  int get_np(){
    return $self->np;
  }
  double get_red(int i) {
    return $self->red[i];
  }
  double get_cth(int i) {
    return $self->cth[i];
  }
  double get_phi(int i) {
    return $self->phi[i];
  }
#ifdef _WITH_WEIGHTS
  double get_weight(int i) {
    return $self->weight[i];
  }
#endif
  np_t get_sum_w(){
    return $self->sum_w;
  }
  np_t get_sum_w2(){
    return $self->sum_w2;
  }
}

%typemap(newfree) Catalog * {
  free_Catalog($1);
}

struct Result {
  int nbins;
  double *x, *corr, *DD, *DR, *RD, *RR;
};

%extend Result{
  int get_nbins(){
    return $self->nbins;
  }
  double get_x(int i) {
    return $self->x[i];
  }
  double get_corr(int i) {
    return $self->corr[i];
  }
  double get_DD(int i) {
    return $self->DD[i];
  }
  double get_DR(int i) {
    return $self->DR[i];
  }
  double get_RD(int i) {
    return $self->RD[i];
  }
  double get_RR(int i) {
    return $self->RR[i];
  }
  ~Result(){
    free($self->x);
    free($self->corr);
    free($self->DD);
    free($self->RD);
    free($self->DR);
    free($self->RR);
  }
}

%typemap(newobject) make_empty_result_struct;
%typemap(newfree) Result * {
  free_result_struct($1);
}
