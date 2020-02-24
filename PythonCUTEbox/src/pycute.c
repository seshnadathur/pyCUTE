#ifdef _CUTE_AS_PYTHON_MODULE
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include "define.h"
#include "common.h"

Result *global_result = NULL;
Catalog *global_galaxy_catalog = NULL;
Catalog *global_galaxy_catalog2 = NULL;
Catalog *global_random_catalog = NULL;

int get_corr_type(){
  return corr_type;
}

int get_do_CCF(){
  return do_CCF;
}

int get_use_randoms(){
  return use_randoms;
}

Result *make_empty_result_struct(){
  int n_bins_all = 0, nx = 0, ny = 0, nz = 0;
  if(corr_type==1){
    n_bins_all=NB_R;
    nx = NB_R;
  } else if(corr_type==2){
    n_bins_all=NB_R*NB_R;
    nx = NB_R;
    ny = NB_R;
  } else if(corr_type==3){
    n_bins_all=NB_R*NB_mu;
    nx = NB_R;
    ny = NB_mu;
  }

  Result *res = malloc(sizeof(Result));
  res->nx   = nx;
  res->ny   = ny;
  res->nz   = nz;
  res->x    = malloc(sizeof(double)*nx);
  res->y    = malloc(sizeof(double)*ny);
  res->z    = malloc(sizeof(double)*nz);
  res->corr = malloc(sizeof(double)*n_bins_all);
  res->D1D1 = malloc(sizeof(double)*n_bins_all);
  res->D1D2 = malloc(sizeof(double)*n_bins_all);
  res->D1R1 = malloc(sizeof(double)*n_bins_all);
  res->D1R2 = malloc(sizeof(double)*n_bins_all);
  res->D2D2 = malloc(sizeof(double)*n_bins_all);
  res->D2R1 = malloc(sizeof(double)*n_bins_all);
  res->D2R2 = malloc(sizeof(double)*n_bins_all);
  res->R1R1 = malloc(sizeof(double)*n_bins_all);
  res->R1R2 = malloc(sizeof(double)*n_bins_all);
  res->R2R2 = malloc(sizeof(double)*n_bins_all);
  return res;
}

void set_result_3d(Result *res, int i, int j, int k, int ind, double x, double y, double z, double corr,
    double D1D1, double D1D2, double D1R1, double D1R2,
    double D2D2, double D2R1, double D2R2,
    double R1R1, double R1R2,
    double R2R2){
  if(res != NULL){
    if(i >= 0) res->x[i] = x;
    if(j >= 0) res->y[j] = y;
    if(k >= 0) res->z[k] = z;
    res->corr[ind]  = corr;
    res->D1D1[ind]  = D1D1;
    res->D1D2[ind]  = D1D2;
    res->D1R1[ind]  = D1R1;
    res->D1R2[ind]  = D1R2;
    res->D2D2[ind]  = D2D2;
    res->D2R1[ind]  = D2R1;
    res->R1R1[ind]  = R1R1;
    res->R1R2[ind]  = R1R2;
    res->R2R2[ind]  = R2R2;
  }
}

void set_result_2d(Result *res, int i, int j, int ind, double x, double y, double corr,
    double D1D1, double D1D2, double D1R1, double D1R2,
    double D2D2, double D2R1, double D2R2,
    double R1R1, double R1R2,
    double R2R2){

  set_result_3d(res, i, j, -1, ind, x, y, 0.0, corr,
      D1D1, D1D2, D1R1, D1R2,
      D2D2, D2R1, D2R2,
      R1R1, R1R2,
      R2R2);
}

void set_result(Result *res, int i, double x, double corr,
    double D1D1, double D1D2, double D1R1, double D1R2,
    double D2D2, double D2R1, double D2R2,
    double R1R1, double R1R2,
    double R2R2){

  set_result_3d(res, i, -1, -1, i, x, 0.0, 0.0, corr,
      D1D1, D1D2, D1R1, D1R2,
      D2D2, D2R1, D2R2,
      R1R1, R1R2,
      R2R2);
}

void free_result_struct(Result *res){
  if(res != NULL){
    free(res->x);
    free(res->y);
    free(res->z);
    free(res->corr);
    free(res->D1D1);
    free(res->D1D2);
    free(res->D1R1);
    free(res->D1R2);
    free(res->D2D2);
    free(res->D2R1);
    free(res->D2R2);
    free(res->R1R1);
    free(res->R1R2);
    free(res->R2R2);
    free(res);
  }
}

Catalog *create_catalog_from_numpy(int n, double *x, int n1, double *y, int n2, double *z){
  if(! ((n == n1) && (n1 == n2))){
    printf("Error: create_catalog_from_numpy inconsistent sizes of the arrays [%i %i %i]\n", n, n1, n2);
    return NULL;
  }
  Catalog *cat = malloc(sizeof(Catalog));
  cat->np = n;
  cat->pos=(double *)malloc(3*cat->np*sizeof(double));
  int i;
  for(i = 0; i < n; i++){
    cat->pos[3*i] = x[i];
    cat->pos[3*i+1] = y[i];
    cat->pos[3*i+2] = z[i];
  }
  return cat;
}

#endif
