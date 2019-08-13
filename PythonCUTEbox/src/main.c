///////////////////////////////////////////////////////////////////////
//                                                                   //
//   Copyright 2012 David Alonso                                     //
//                                                                   //
//                                                                   //
// This file is part of CUTE.                                        //
//                                                                   //
// CUTE is free software: you can redistribute it and/or modify it   //
// under the terms of the GNU General Public License as published by //
// the Free Software Foundation, either version 3 of the License, or //
// (at your option) any later version.                               //
//                                                                   //
// CUTE is distributed in the hope that it will be useful, but       //
// WITHOUT ANY WARRANTY; without even the implied warranty of        //
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU //
// General Public License for more details.                          //
//                                                                   //
// You should have received a copy of the GNU General Public License //
// along with CUTE.  If not, see <http://www.gnu.org/licenses/>.     //
//                                                                   //
///////////////////////////////////////////////////////////////////////

/*********************************************************************/
//                               Main                                //
/*********************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include "define.h"
#include "common.h"
#include "correlator.h"
#include "io.h"
#include "pm.h"
#include "neighbors.h"
#include "tree.h"

void write_CF(char *fname,double *corr,double *ercorr,
    unsigned long long *DD)
{
  //////
  // Writes correlation function to file fname
  FILE *fo;
  int ii;

  fo=fopen(fname,"w");
  if(fo==NULL) {
    char oname[64]="output_CUTE.dat";
    fprintf(stderr,"CUTE: Error opening output file %s",fname);
    fprintf(stderr,", using ./output_CUTE.dat");
    fo=fopen(oname,"w");
    if(fo==NULL) error_open_file(oname);
  }

  for(ii=0;ii<NB_R;ii++) {
    double rr;
    int ind = ii; (void)ind;
#ifdef _LOGBIN
    rr=pow(10,((ii+0.5)-NB_R)/N_LOGINT+LOG_R_MAX);
#else //_LOGBIN
    rr=(ii+0.5)/(NB_R*I_R_MAX);
#endif //_LOGBIN
    fprintf(fo,"%lE %lE %lE %llu \n",
        rr,corr[ii],ercorr[ii],DD[ii]);
#ifdef _CUTE_AS_PYTHON_MODULE
    set_result(global_result, ii, rr, corr[ind],
        (double)DD[ind], 0.0, 0.0, 0.0,
        0.0,             0.0, 0.0,
        0.0,             0.0,
        0.0);
#endif
  }
  fclose(fo);

  printf("\n");
}

void write_CF_w_rand(char *fname,double *corr,double *ercorr,
    unsigned long long *DD,unsigned long long *DR,unsigned long long *RR)
{
  //////
  // Writes correlation function to file fname
  FILE *fo;
  int ii;

  fo=fopen(fname,"w");
  if(fo==NULL) {
    char oname[64]="output_CUTE.dat";
    fprintf(stderr,"CUTE: Error opening output file %s",fname);
    fprintf(stderr,", using ./output_CUTE.dat");
    fo=fopen(oname,"w");
    if(fo==NULL) error_open_file(oname);
  }

  for(ii=0;ii<NB_R;ii++) {
    double rr;
    int ind = ii; (void)ind;
#ifdef _LOGBIN
    rr=pow(10,((ii+0.5)-NB_R)/N_LOGINT+LOG_R_MAX);
#else //_LOGBIN
    rr=(ii+0.5)/(NB_R*I_R_MAX);
#endif //_LOGBIN
    fprintf(fo,"%lE %lE %lE %llu %llu %llu \n",
        rr,corr[ii],ercorr[ii],DD[ii],DR[ii],RR[ii]);
#ifdef _CUTE_AS_PYTHON_MODULE
    set_result(global_result, ii, rr, corr[ind],
        (double)DD[ind], 0.0, (double)DR[ind], 0.0,
        0.0,             0.0, 0.0,
        (double)RR[ind], 0.0,
        0.0);
#endif
  }
  fclose(fo);

  printf("\n");
}

void write_CCF(char *fname,double *corr,double *ercorr,unsigned long long *DD,
    unsigned long long *D1R,unsigned long long *D2R,unsigned long long *RR)
{
  //////
  // Writes correlation function to file fname
  FILE *fo;
  int ii;

  fo=fopen(fname,"w");
  if(fo==NULL) {
    char oname[64]="output_CUTE.dat";
    fprintf(stderr,"CUTE: Error opening output file %s",fname);
    fprintf(stderr,", using ./output_CUTE.dat");
    fo=fopen(oname,"w");
    if(fo==NULL) error_open_file(oname);
  }

  for(ii=0;ii<NB_R;ii++) {
    double rr;
    int ind = ii; (void)ind;
#ifdef _LOGBIN
    rr=pow(10,((ii+0.5)-NB_R)/N_LOGINT+LOG_R_MAX);
#else //_LOGBIN
    rr=(ii+0.5)/(NB_R*I_R_MAX);
#endif //_LOGBIN
    fprintf(fo,"%lE %lE %lE %llu %llu %llu %llu \n",
        rr,corr[ii],ercorr[ii],DD[ii],D1R[ii],D2R[ii],RR[ii]);
#ifdef _CUTE_AS_PYTHON_MODULE
    set_result(global_result, ii, rr, corr[ind],
        0.0,             (double)DD[ind], (double)D1R[ind], 0.0,
        0.0,             (double)D2R[ind], 0.0,
        (double)RR[ind], 0.0,
        0.0);
#endif
  }
  fclose(fo);

  printf("\n");
}

void write_3d_CF_w_rand(char *fname,double *corr,double *ercorr,unsigned long long *DD,
    unsigned long long *DR,unsigned long long *RR)
{
  ///////
  // Writes 3d correlation function to file fname
  FILE *fo;
  int ii,jj;

  fo=fopen(fname,"w");
  if(fo==NULL) {
    char oname[64]="output_CUTE.dat";
    fprintf(stderr,"CUTE: Error opening output file %s",fname);
    fprintf(stderr,", using ./output_CUTE.dat");
    fo=fopen(oname,"w");
    if(fo==NULL) error_open_file(oname);
  }

  if(corr_type==1) {
    fprintf(fo,"# r xi(r) sigma_xi DD DR RR\n");
    for(ii=0;ii<NB_R;ii++) {
      double rr;
      int ind = ii; (void)ind;
#ifdef _LOGBIN
      rr=pow(10,((ii+0.5)-NB_R)/N_LOGINT+LOG_R_MAX);
#else //_LOGBIN
      rr=(ii+0.5)/(NB_R*I_R_MAX);
#endif //_LOGBIN
      fprintf(fo,"%lE %lE %lE %llu %llu %llu \n",
          rr,corr[ii],ercorr[ii],DD[ii],DR[ii],RR[ii]);
#ifdef _CUTE_AS_PYTHON_MODULE
      set_result(global_result, ii, rr, corr[ind],
          (double)DD[ind], 0.0, (double)DR[ind], 0.0,
          0.0,             0.0, 0.0,
          (double)RR[ind], 0.0,
          0.0);
#endif
    }
  }
  else if(corr_type==2) {
    fprintf(fo,"# r_t r_l xi(r_l,r_t) DD DR RR\n");
    for(ii=0;ii<NB_R;ii++) {
      for(jj=0;jj<NB_R;jj++) {
        double rl,rt;
        int ind=ii+NB_R*jj; (void)ind;
#ifdef _LOGBIN
        rl=pow(10,((ii+0.5)-NB_R)/N_LOGINT+LOG_R_MAX);
        rt=pow(10,((jj+0.5)-NB_R)/N_LOGINT+LOG_R_MAX);
#else //_LOGBIN
        rl=(ii+0.5)/(NB_R*I_R_MAX);
        rt=(jj+0.5)/(NB_R*I_R_MAX);
#endif //_LOGBIN
        fprintf(fo,"%lE %lE %lE %llu %llu %llu \n",
            rl,rt,corr[ii+NB_R*jj],DD[ii+NB_R*jj],DR[ii+NB_R*jj],RR[ii+NB_R*jj]);
#ifdef _CUTE_AS_PYTHON_MODULE
        set_result_2d(global_result, ii, jj, ind, rl, rt, corr[ind],
            (double)DD[ind], 0.0, (double)DR[ind], 0.0,
            0.0,             0.0, 0.0,
            (double)RR[ind], 0.0,
            0.0);
#endif
      }
    }
  }
  else if(corr_type==3) {
    fprintf(fo,"# r mu xi(r,mu) DD DR RR\n");
    for(ii=0;ii<NB_R;ii++) {
      for(jj=0;jj<NB_mu;jj++) {
        double r,mu;
        int ind = jj+NB_R*ii; (void)ind;
        r=(ii+0.5)/(NB_R*I_R_MAX);
        mu=(jj+0.5)/(NB_mu);
        fprintf(fo,"%lE %lE %lE %llu %llu %llu \n",
            r,mu,corr[jj+NB_R*ii],DD[jj+NB_R*ii],DR[jj+NB_R*ii],RR[jj+NB_R*ii]);
#ifdef _CUTE_AS_PYTHON_MODULE
        set_result_2d(global_result, ii, jj, ind, r, mu, corr[ind],
            (double)DD[ind], 0.0, (double)DR[ind], 0.0,
            0.0,             0.0, 0.0,
            (double)RR[ind], 0.0,
            0.0);
#endif
      }
    }
  }
}

void write_3d_CCF_w_rand(char *fname,double *corr,double *ercorr,unsigned long long *DD,
    unsigned long long *D1R,unsigned long long *D2R,unsigned long long *RR)
{
  ///////
  // Writes 3d correlation function to file fname
  FILE *fo;
  int ii,jj;

  fo=fopen(fname,"w");
  if(fo==NULL) {
    char oname[64]="output_CUTE.dat";
    fprintf(stderr,"CUTE: Error opening output file %s",fname);
    fprintf(stderr,", using ./output_CUTE.dat");
    fo=fopen(oname,"w");
    if(fo==NULL) error_open_file(oname);
  }

  if(corr_type==1) {
    fprintf(fo,"# r xi(r) sigma_xi D1D2 D1R D2R RR\n");
    for(ii=0;ii<NB_R;ii++) {
      double rr;
      int ind = ii; (void)ind;
#ifdef _LOGBIN
      rr=pow(10,((ii+0.5)-NB_R)/N_LOGINT+LOG_R_MAX);
#else //_LOGBIN
      rr=(ii+0.5)/(NB_R*I_R_MAX);
#endif //_LOGBIN
      fprintf(fo,"%lE %lE %lE %llu %llu %llu %llu \n",
          rr,corr[ii],ercorr[ii],DD[ii],D1R[ii],D2R[ii],RR[ii]);
#ifdef _CUTE_AS_PYTHON_MODULE
      set_result(global_result, ii, rr, corr[ind],
          0.0,                   (double)DD[ind],  (double)D1R[ind], 0.0,
          0.0,                   (double)D2R[ind], 0.0,
          (double)RR[ind],       0.0,
          0.0);
#endif
    }
  }
  else if(corr_type==2) {
    fprintf(fo,"# r_t r_l xi(r_l,r_t) D1D2 D1R D2R RR\n");
    for(ii=0;ii<NB_R;ii++) {
      for(jj=0;jj<NB_R;jj++) {
        double rl,rt;
        int ind = ii+NB_R*jj; (void)ind;
#ifdef _LOGBIN
        rt=pow(10,((ii+0.5)-NB_R)/N_LOGINT+LOG_R_MAX);
        rl=pow(10,((jj+0.5)-NB_R)/N_LOGINT+LOG_R_MAX);
#else //_LOGBIN
        rt=(ii+0.5)/(NB_R*I_R_MAX);
        rl=(jj+0.5)/(NB_R*I_R_MAX);
#endif //_LOGBIN
        fprintf(fo,"%lE %lE %lE %llu %llu %llu %llu \n",
            rt,rl,corr[ii+NB_R*jj],DD[ii+NB_R*jj],D1R[ii+NB_R*jj],D2R[ii+NB_R*jj],RR[ii+NB_R*jj]);
#ifdef _CUTE_AS_PYTHON_MODULE
        set_result_2d(global_result, ii, jj, ind, rt, rl, corr[ind],
            0.0,             (double)DD[ind],  (double)D1R[ind], 0.0,
            0.0,             (double)D2R[ind], 0.0,
            (double)RR[ind], 0.0,
            0.0);
#endif
      }
    }
  }
  else if(corr_type==3) {
    fprintf(fo,"# r mu xi(r,mu)D1D2 D1R D2R RR\n");
    for(ii=0;ii<NB_R;ii++) {
      for(jj=0;jj<NB_mu;jj++) {
        double r,mu;
        int ind = jj+NB_mu*ii; (void)ind;
        r=(ii+0.5)/(NB_R*I_R_MAX);
        mu=(jj+0.5)/(NB_mu);
        fprintf(fo,"%lE %lE %lE %llu %llu %llu %llu \n",
            r,mu,corr[jj+NB_mu*ii],DD[jj+NB_mu*ii],D1R[jj+NB_mu*ii],D2R[jj+NB_mu*ii],RR[jj+NB_mu*ii]);
#ifdef _CUTE_AS_PYTHON_MODULE
        set_result_2d(global_result, ii, jj, ind, r, mu, corr[ind],
            0.0,             (double)DD[ind],  (double)D1R[ind], 0.0,
            0.0,             (double)D2R[ind], 0.0,
            (double)RR[ind], 0.0,
            0.0);
#endif
      }
    }
  }
}

void write_3d_CF(char *fname,double *corr,double *ercorr,unsigned long long *DD)
{
  ///////
  // Writes 3d correlation function to file fname
  FILE *fo;
  int ii,jj;

  fo=fopen(fname,"w");
  if(fo==NULL) {
    char oname[64]="output_CUTE.dat";
    fprintf(stderr,"CUTE: Error opening output file %s",fname);
    fprintf(stderr,", using ./output_CUTE.dat");
    fo=fopen(oname,"w");
    if(fo==NULL) error_open_file(oname);
  }

  if(corr_type==1) {
    fprintf(fo,"# r xi(r) sigma_xi DD\n");
    for(ii=0;ii<NB_R;ii++) {
      double rr;
      int ind = ii; (void)ind;
#ifdef _LOGBIN
      rr=pow(10,((ii+0.5)-NB_R)/N_LOGINT+LOG_R_MAX);
#else //_LOGBIN
      rr=(ii+0.5)/(NB_R*I_R_MAX);
#endif //_LOGBIN
      fprintf(fo,"%lE %lE %lE %llu \n",
          rr,corr[ii],ercorr[ii],DD[ii]);
#ifdef _CUTE_AS_PYTHON_MODULE
      set_result(global_result, ii, rr, corr[ind],
          (double)DD[ind], 0.0, 0.0, 0.0,
          0.0,             0.0, 0.0,
          0.0,             0.0,
          0.0);
#endif
    }
  }
  else if(corr_type==2) {
    fprintf(fo,"# r_t r_l xi(r_l,r_t) DD\n");
    for(ii=0;ii<NB_R;ii++) {
      for(jj=0;jj<NB_R;jj++) {
        double rl,rt;
        int ind = ii+NB_R*jj; (void)ind;
#ifdef _LOGBIN
        rl=pow(10,((ii+0.5)-NB_R)/N_LOGINT+LOG_R_MAX);
        rt=pow(10,((jj+0.5)-NB_R)/N_LOGINT+LOG_R_MAX);
#else //_LOGBIN
        rl=(ii+0.5)/(NB_R*I_R_MAX);
        rt=(jj+0.5)/(NB_R*I_R_MAX);
#endif //_LOGBIN
        fprintf(fo,"%lE %lE %lE %llu \n",
            rl,rt,corr[ii+NB_R*jj],DD[ii+NB_R*jj]);
#ifdef _CUTE_AS_PYTHON_MODULE
        set_result_2d(global_result, ii, jj, ind, rl, rt, corr[ind],
            (double)DD[ind], 0.0, 0.0, 0.0,
            0.0,             0.0, 0.0,
            0.0,             0.0,
            0.0);
#endif
      }
    }
  }
  else if(corr_type==3) {
    fprintf(fo,"# r mu xi(r,mu) DD\n");
    for(ii=0;ii<NB_R;ii++) {
      for(jj=0;jj<NB_mu;jj++) {
        double r,mu;
        int ind = ii+NB_R*jj; (void)ind;
        r=(ii+0.5)/(NB_R*I_R_MAX);
        mu=(jj+0.5)/(NB_mu);
        fprintf(fo,"%lE %lE %lE %llu \n",
            r,mu,corr[jj+NB_R*ii],DD[jj+NB_R*ii]);
#ifdef _CUTE_AS_PYTHON_MODULE
        set_result_2d(global_result, ii, jj, ind, r, mu, corr[ind],
            (double)DD[ind], 0.0, 0.0, 0.0,
            0.0,             0.0, 0.0,
            0.0,             0.0,
            0.0);
#endif
      }
    }
  }
}

void write_3d_CCF(char *fname,double *corr,double *ercorr,unsigned long long *DD)
{
  ///////
  // Writes 3d correlation function to file fname
  FILE *fo;
  int ii,jj;

  fo=fopen(fname,"w");
  if(fo==NULL) {
    char oname[64]="output_CUTE.dat";
    fprintf(stderr,"CUTE: Error opening output file %s",fname);
    fprintf(stderr,", using ./output_CUTE.dat");
    fo=fopen(oname,"w");
    if(fo==NULL) error_open_file(oname);
  }

  if(corr_type==1) {
    fprintf(fo,"# r xi(r) sigma_xi D1D2 \n");
    for(ii=0;ii<NB_R;ii++) {
      double rr;
      int ind = ii; (void)ind;
#ifdef _LOGBIN
      rr=pow(10,((ii+0.5)-NB_R)/N_LOGINT+LOG_R_MAX);
#else //_LOGBIN
      rr=(ii+0.5)/(NB_R*I_R_MAX);
#endif //_LOGBIN
      fprintf(fo,"%lE %lE %lE %llu \n",
          rr,corr[ii],ercorr[ii],DD[ii]);
#ifdef _CUTE_AS_PYTHON_MODULE
      set_result(global_result, ii, rr, corr[ind],
          0.0,             (double)DD[ind], 0.0, 0.0,
          0.0,             0.0, 0.0,
          0.0,             0.0,
          0.0);
#endif
    }
  }
  else if(corr_type==2) {
    fprintf(fo,"# r_t r_l xi(r_l,r_t) D1D2\n");
    for(ii=0;ii<NB_R;ii++) {
      for(jj=0;jj<NB_R;jj++) {
        double rl,rt;
        int ind = ii+NB_R*jj; (void)ind;
#ifdef _LOGBIN
        rl=pow(10,((ii+0.5)-NB_R)/N_LOGINT+LOG_R_MAX);
        rt=pow(10,((jj+0.5)-NB_R)/N_LOGINT+LOG_R_MAX);
#else //_LOGBIN
        rl=(ii+0.5)/(NB_R*I_R_MAX);
        rt=(jj+0.5)/(NB_R*I_R_MAX);
#endif //_LOGBIN
        fprintf(fo,"%lE %lE %lE %llu \n",
            rl,rt,corr[ii+NB_R*jj],DD[ii+NB_R*jj]);
#ifdef _CUTE_AS_PYTHON_MODULE
        set_result_2d(global_result, ii, jj, ind, rl, rt, corr[ind],
            (double)DD[ind], 0.0, 0.0, 0.0,
            0.0,             0.0, 0.0,
            0.0,             0.0,
            0.0);
#endif
      }
    }
  }
  else if(corr_type==3) {
    fprintf(fo,"# r mu xi(r,mu) D1D2\n");
    for(ii=0;ii<NB_R;ii++) {
      for(jj=0;jj<NB_mu;jj++) {
        double r,mu;
        int ind = jj+NB_mu*ii; (void)ind;
        r=(ii+0.5)/(NB_R*I_R_MAX);
        mu=(jj+0.5)/(NB_mu);
        fprintf(fo,"%lE %lE %lE %llu \n",
            r,mu,corr[jj+NB_mu*ii],DD[jj+NB_mu*ii]);
#ifdef _CUTE_AS_PYTHON_MODULE
        set_result_2d(global_result, ii, jj, ind, r, mu, corr[ind],
            (double)DD[ind], 0.0, 0.0, 0.0,
            0.0,             0.0, 0.0,
            0.0,             0.0,
            0.0);
#endif
      }
    }
  }
}

void run_monopole_corr_bf(void)
{
  //////
  // Main routine for monopole using brute-force
  lint n_dat;
  Catalog *cat_dat;
  unsigned long long DD[NB_R];
  double corr[NB_R],ercorr[NB_R];

  timer(4);

#ifdef _VERBOSE
  printf("*** Correlation function parameters: \n");
  printf(" - Range: %.3lf < r < %.3lf (Mpc/h)\n",0.,1./(I_R_MAX));
  printf(" - #bins: %d\n",NB_R);
#ifdef _LOGBIN
  printf(" - Using logarithmic binning with %d bins per decade \n",N_LOGINT);
#else
  printf(" - Resolution: Dr = %.3lf (Mpc/h)\n",1./(I_R_MAX*NB_R));
#endif
  printf(" - Using a brute-force approach \n");
  printf("\n");
#endif

  //Read data
#ifdef _CUTE_AS_PYTHON_MODULE
  if(global_galaxy_catalog == NULL){
    cat_dat=read_catalog(fnameData,&n_dat);
  } else {
    cat_dat=global_galaxy_catalog;
    n_dat=global_galaxy_catalog->np;
  }
#else
  cat_dat=read_catalog(fnameData,&n_dat);
#endif

#ifdef _DEBUG
  write_cat(cat_dat,"debug_DatCat.dat");
#endif

  printf("*** Correlating\n");
  timer(0);
  corr_mono_box_bf(cat_dat->np,cat_dat->pos,DD);
  timer(1);
  printf("\n");

  printf("*** Writing output \n");
  make_CF(DD,n_dat,corr,ercorr);
  write_CF(fnameOut,corr,ercorr,DD);

  printf("*** Cleaning up \n");
#ifdef _CUTE_AS_PYTHON_MODULE
  if(global_galaxy_catalog == NULL)
#endif
    free_catalog(cat_dat);
  printf("\n");

  timer(5);
}

void run_monopole_corr_neighbors(void)
{
  //////
  // Main routine for monopole using neighbor boxes, no randoms
  lint n_dat;
  int nside;
  Catalog *cat_dat;
  NeighborBox *boxes;
  unsigned long long DD[NB_R];
  double corr[NB_R],ercorr[NB_R];

  timer(4);

#ifdef _VERBOSE
  printf("*** Correlation function parameters: \n");
  printf(" - Range: %.3lf < r < %.3lf (Mpc/h)\n",0.,1./(I_R_MAX));
  printf(" - #bins: %d\n",NB_R);
#ifdef _LOGBIN
  printf(" - Using logarithmic binning with %d bins per decade \n",N_LOGINT);
#else //_LOGBIN
  printf(" - Resolution: Dr = %.3lf (Mpc/h)\n",1./(I_R_MAX*NB_R));
#endif //_LOGBIN
  printf("\n");
#endif //_LOGBIN

  //Read data
#ifdef _CUTE_AS_PYTHON_MODULE
  if(global_galaxy_catalog == NULL){
    cat_dat=read_catalog(fnameData,&n_dat);
  } else {
    cat_dat=global_galaxy_catalog;
    n_dat=global_galaxy_catalog->np;
  }
#else
  cat_dat=read_catalog(fnameData,&n_dat);
#endif
  nside=optimal_nside(l_box,1./I_R_MAX,cat_dat->np);
  boxes=catalog_to_boxes(nside,*cat_dat);

#ifdef _DEBUG
  write_cat(cat_dat,"debug_DatCat.dat");
#endif

  printf("*** Correlating\n");
  timer(0);
  corr_mono_boxes(nside,boxes,DD);
  make_3d_CF(DD,n_dat,corr,ercorr);
  //corr_mono_box_neighbors(nside,boxes,cat_dat.np,cat_dat.pos,DD);
  //make_CF_double(DD,n_dat,corr,ercorr);
  timer(1);
  printf("\n");

  printf("*** Writing output \n");
  write_3d_CF(fnameOut,corr,ercorr,DD);

  printf("*** Cleaning up \n");
#ifdef _CUTE_AS_PYTHON_MODULE
  if(global_galaxy_catalog == NULL)
#endif
    free_catalog(cat_dat);
  free_boxes(nside,boxes);
  printf("\n");

  timer(5);
}

void run_monopole_CCF(int use_randoms,int reuse_randoms)
{
  //////
  // Routine for monopole cross-correlation using neighbor boxes, using randoms to account for partial box
  lint n_dat1,n_dat2,n_rand;
  int nside,nside1,nside2,i;
  Catalog *cat_dat1, *cat_dat2, *rand_dat;
  NeighborBox *data1_boxes, *data2_boxes, *rand_boxes;
  unsigned long long D1D2[NB_R], D1R[NB_R], D2R[NB_R], RR[NB_R];
  double corr[NB_R],ercorr[NB_R];

  timer(4);

#ifdef _VERBOSE
  printf("*** Correlation function parameters: \n");
  printf(" - Range: %.3lf < r < %.3lf (Mpc/h)\n",0.,1./(I_R_MAX));
  printf(" - #bins: %d\n",NB_R);
#ifdef _LOGBIN
  printf(" - Using logarithmic binning with %d bins per decade \n",N_LOGINT);
#else //_LOGBIN
  printf(" - Resolution: Dr = %.3lf (Mpc/h)\n",1./(I_R_MAX*NB_R));
#endif //_LOGBIN
  printf("\n");
#endif //_LOGBIN

  //Read data
#ifdef _CUTE_AS_PYTHON_MODULE
  if(global_galaxy_catalog == NULL){
    cat_dat1=read_catalog(fnameData,&n_dat1);
  } else {
    cat_dat1=global_galaxy_catalog;
    n_dat1=global_galaxy_catalog->np;
  }
#else
  cat_dat1=read_catalog(fnameData,&n_dat1);
#endif
  nside1=optimal_nside(l_box,1./I_R_MAX,cat_dat1->np);

#ifdef _CUTE_AS_PYTHON_MODULE
  if(global_galaxy_catalog2 == NULL){
    cat_dat2=read_catalog(fnameData2,&n_dat2);
  } else {
    cat_dat2=global_galaxy_catalog2;
    n_dat2=global_galaxy_catalog2->np;
  }
#else
  cat_dat2=read_catalog(fnameData2,&n_dat2);
#endif
  nside2=optimal_nside(l_box,1./I_R_MAX,cat_dat2->np);

  nside=MAX(nside1,nside2);
  data1_boxes=catalog_to_boxes(nside,*cat_dat1);
  data2_boxes=catalog_to_boxes(nside,*cat_dat2);

  if(use_randoms) {
    //Read randoms data
#ifdef _CUTE_AS_PYTHON_MODULE
    if(global_random_catalog == NULL){
      rand_dat=read_catalog(fnameRand,&n_rand);
    } else {
      rand_dat=global_random_catalog;
      n_rand=global_random_catalog->np;
    }
#else
    rand_dat=read_catalog(fnameRand,&n_rand);
#endif
    rand_boxes=catalog_to_boxes(nside,*rand_dat);
  }

#ifdef _DEBUG
  write_cat(cat_dat1,"debug_Dat1Cat.dat");
  write_cat(cat_dat1,"debug_Dat2Cat.dat");
  if(use_randoms) write_cat(rand_dat,"debug_RandCat.dat");
#endif

  printf("*** Correlating\n");
  timer(0);
  crosscorr_mono_box_neighbors(nside,data1_boxes,data2_boxes,D1D2);
  if(use_randoms) {
    crosscorr_mono_box_neighbors(nside,data1_boxes,rand_boxes,D1R);
    if(reuse_randoms==2) {
      for(i=0;i<NB_R;i++) {
        D2R[i]=0;
        RR[i]=0;
        corr[i]=0;
        ercorr[i]=0;
      }
      timer(1);
      printf("\n");
      printf("*** Writing output without D2R and RR - adjust output files with pre-calculated values \n");
    }
    else if(reuse_randoms==1) {
      crosscorr_mono_box_neighbors(nside,data2_boxes,rand_boxes,D2R);
      for(i=0;i<NB_R;i++) {
        RR[i]=0;
        corr[i]=0;
        ercorr[i]=0;
      }
      timer(1);
      printf("\n");
      printf("*** Writing output without RR - adjust output files with pre-calculated values \n");
    }
    else {
      crosscorr_mono_box_neighbors(nside,data2_boxes,rand_boxes,D2R);
      corr_mono_boxes(nside,rand_boxes,RR);
      timer(1);
      printf("\n");
      printf("*** Writing output \n");
      make_3d_CCF_w_rand(D1D2,D1R,D2R,RR,n_dat1,n_dat2,n_rand,corr,ercorr);
    }
    write_3d_CCF_w_rand(fnameOut,corr,ercorr,D1D2,D1R,D2R,RR);
  }
  else {    // ie, not using randoms
    timer(1);
    printf("\n");
    printf("*** Writing output \n");
    make_3d_CCF(D1D2,n_dat1,n_dat2,corr,ercorr);
    write_3d_CCF(fnameOut,corr,ercorr,D1D2);
  }

  printf("*** Cleaning up \n");

#ifdef _CUTE_AS_PYTHON_MODULE
  if(global_galaxy_catalog == NULL)
#endif
    free_catalog(cat_dat1);

#ifdef _CUTE_AS_PYTHON_MODULE
  if(global_galaxy_catalog2 == NULL)
#endif
    free_catalog(cat_dat2);

  free_boxes(nside,data1_boxes);
  free_boxes(nside,data2_boxes);
  if(use_randoms) {
#ifdef _CUTE_AS_PYTHON_MODULE
    if(global_random_catalog == NULL)
#endif
      free_catalog(rand_dat);

    free_boxes(nside,rand_boxes);
  }
  printf("\n");

  timer(5);
}

void run_monopole_corr_neighbors_w_rand(void)
{
  //////
  // Routine for monopole using neighbor boxes, using randoms to account for partial box
  lint n_dat,n_rand;
  int nside;
  Catalog *cat_dat, *rand_dat;
  NeighborBox *data_boxes, *rand_boxes;
  unsigned long long DD[NB_R], DR[NB_R], RR[NB_R];
  double corr[NB_R],ercorr[NB_R];

  timer(4);

#ifdef _VERBOSE
  printf("*** Correlation function parameters: \n");
  printf(" - Range: %.3lf < r < %.3lf (Mpc/h)\n",0.,1./(I_R_MAX));
  printf(" - #bins: %d\n",NB_R);
#ifdef _LOGBIN
  printf(" - Using logarithmic binning with %d bins per decade \n",N_LOGINT);
#else //_LOGBIN
  printf(" - Resolution: Dr = %.3lf (Mpc/h)\n",1./(I_R_MAX*NB_R));
#endif //_LOGBIN
  printf("\n");
#endif //_VERBOSE

  //Read data
#ifdef _CUTE_AS_PYTHON_MODULE
  if(global_galaxy_catalog == NULL){
    cat_dat=read_catalog(fnameData,&n_dat);
  } else {
    cat_dat=global_galaxy_catalog;
    n_dat=global_galaxy_catalog->np;
  }
#else
  cat_dat=read_catalog(fnameData,&n_dat);
#endif
  nside=optimal_nside(l_box,1./I_R_MAX,cat_dat->np);
  data_boxes=catalog_to_boxes(nside,*cat_dat);

  //Read randoms data
#ifdef _CUTE_AS_PYTHON_MODULE
  if(global_random_catalog == NULL){
    rand_dat=read_catalog(fnameRand,&n_rand);
  } else {
    rand_dat=global_random_catalog;
    n_rand=global_random_catalog->np;
  }
#else
  rand_dat=read_catalog(fnameRand,&n_rand);
#endif
  rand_boxes=catalog_to_boxes(nside,*rand_dat);

#ifdef _DEBUG
  write_cat(cat_dat,"debug_DatCat.dat");
  write_cat(rand_dat,"debug_RandCat.dat");
#endif

  printf("*** Correlating\n");
  timer(0);
  corr_mono_boxes(nside,data_boxes,DD);
  corr_mono_boxes(nside,rand_boxes,RR);
  crosscorr_mono_box_neighbors(nside,data_boxes,rand_boxes,DR);
  timer(1);
  printf("\n");

  printf("*** Writing output \n");
  make_3d_CF_w_rand(DD,DR,RR,n_dat,n_rand,corr,ercorr);
  write_3d_CF_w_rand(fnameOut,corr,ercorr,DD,DR,RR);
  //make_CF(DD,n_dat,corr,ercorr);
  //write_CF(fnameOut,corr,ercorr,DD);

  printf("*** Cleaning up \n");
#ifdef _CUTE_AS_PYTHON_MODULE
  if(global_galaxy_catalog == NULL)
#endif
    free_catalog(cat_dat);
#ifdef _CUTE_AS_PYTHON_MODULE
  if(global_random_catalog == NULL)
#endif
    free_catalog(rand_dat);
  free_boxes(nside,data_boxes);
  free_boxes(nside,rand_boxes);
  printf("\n");

  timer(5);
}

void run_3d_ps_auto_corr_boxes(int use_randoms)
{
  //////
  // Runs xi(pi,sigma) in brute-force mode
  lint n_dat,n_rand;
  int nside;
  Catalog *cat_dat, *cat_rand;
  NeighborBox *data_boxes, *rand_boxes;
  unsigned long long DD[NB_R*NB_R], DR[NB_R*NB_R], RR[NB_R*NB_R];
  double corr[NB_R*NB_R],ercorr[NB_R*NB_R];

  timer(4);

#ifdef _VERBOSE
  printf("*** 3D correlation function (pi,sigma): \n");
  printf(" - Range: (%.3lf,%.3lf) < (pi,sigma) < (%.3lf,%.3lf) Mpc/h\n",
      0.,0.,1./I_R_MAX,1./I_R_MAX);
  printf(" - #bins: (%d,%d)\n",NB_R,NB_R);
#ifdef _LOGBIN
  printf(" - Using logarithmic binning with %d bins per decade (log binning not recommended for (pi,sigma)!)\n",N_LOGINT);
#else //_LOGBIN
  printf(" - Resolution: Dr = %.3lf (Mpc/h)\n",1./(I_R_MAX*NB_R));
#endif //_LOGBIN
  printf(" - Using a brute-force approach \n");
  printf("\n");
#endif // _VERBOSE

  //Read data
#ifdef _CUTE_AS_PYTHON_MODULE
  if(global_galaxy_catalog == NULL){
    cat_dat=read_catalog(fnameData,&n_dat);
  } else {
    cat_dat=global_galaxy_catalog;
    n_dat=global_galaxy_catalog->np;
  }
#else
  cat_dat=read_catalog(fnameData,&n_dat);
#endif
  nside=optimal_nside(l_box,1./I_R_MAX,cat_dat->np);
  data_boxes=catalog_to_boxes(nside,*cat_dat);

  if(use_randoms) {
    //Read randoms data
#ifdef _CUTE_AS_PYTHON_MODULE
    if(global_random_catalog == NULL){
      cat_rand=read_catalog(fnameRand,&n_rand);
    } else {
      cat_rand=global_random_catalog;
      n_rand=global_random_catalog->np;
    }
#else
    cat_rand=read_catalog(fnameRand,&n_rand);
#endif
    rand_boxes=catalog_to_boxes(nside,*cat_rand);
  }

#ifdef _DEBUG
  write_cat(cat_dat,"debug_DatCat.dat");
  if(use_randoms) write_cat(cat_rand,"debug_RandCat.dat");
#endif

  printf("*** Correlating \n");
  printf(" - Auto-correlating data \n");
  timer(0);
  auto_3d_ps_boxes(nside,data_boxes,DD);
  if(use_randoms) {
    timer(2);
    printf(" - Auto-correlating random \n");
    auto_3d_ps_boxes(nside,rand_boxes,RR);
    timer(2);
    printf(" - Cross-correlating \n");
    cross_3d_ps_boxes(nside,data_boxes,rand_boxes,DR);
    timer(1);

    printf("*** Writing output\n");
    make_3d_CF_w_rand(DD,DR,RR,n_dat,n_rand,corr,ercorr);
    write_3d_CF_w_rand(fnameOut,corr,ercorr,DD,DR,RR);
  }
  else {    // ie, not using randoms
    timer(1);
    printf("*** Writing output\n");
    make_3d_CF(DD,n_dat,corr,ercorr);
    write_3d_CF(fnameOut,corr,ercorr,DD);
  }

  printf("*** Cleaning up\n");
#ifdef _CUTE_AS_PYTHON_MODULE
  if(global_galaxy_catalog == NULL)
#endif
    free_catalog(cat_dat);
  free_boxes(nside,data_boxes);
  if(use_randoms) {
#ifdef _CUTE_AS_PYTHON_MODULE
    if(global_random_catalog == NULL)
#endif
      free_catalog(cat_rand);
    free_boxes(nside,rand_boxes);
  }
  printf("\n");

  timer(5);
}

void run_3d_ps_cross_corr_boxes(int use_randoms,int reuse_randoms)
{
  //////
  // Runs xi(pi,sigma) cross-correlation in brute-force mode
  lint n_dat1,n_dat2,n_rand;
  int nside,nside1,nside2,i;
  Catalog *cat_dat1,*cat_dat2,*cat_rand;
  NeighborBox *data1_boxes, *data2_boxes, *rand_boxes;
  unsigned long long D1D2[NB_R*NB_R], D1R[NB_R*NB_R], D2R[NB_R*NB_R], RR[NB_R*NB_R];
  double corr[NB_R*NB_R],ercorr[NB_R*NB_R];

  timer(4);

#ifdef _VERBOSE
  printf("*** 3D cross-correlation function (pi,sigma): \n");
  printf(" - Range: (%.3lf,%.3lf) < (pi,sigma) < (%.3lf,%.3lf) Mpc/h\n",
      0.,0.,1./I_R_MAX,1./I_R_MAX);
  printf(" - #bins: (%d,%d)\n",NB_R,NB_R);
#ifdef _LOGBIN
  printf(" - Using logarithmic binning with %d bins per decade (log binning not recommended for (pi,sigma)!)\n",N_LOGINT);
#else //_LOGBIN
  printf(" - Resolution: Dr = %.3lf (Mpc/h)\n",1./(I_R_MAX*NB_R));
#endif //_LOGBIN
  printf(" - Using a brute-force approach \n");
  printf("\n");
#endif // _VERBOSE

  //Read data
#ifdef _CUTE_AS_PYTHON_MODULE
  if(global_galaxy_catalog == NULL){
    cat_dat1=read_catalog(fnameData,&n_dat1);
  } else {
    cat_dat1=global_galaxy_catalog;
    n_dat1=global_galaxy_catalog->np;
  }
#else
  cat_dat1=read_catalog(fnameData,&n_dat1);
#endif
  nside1=optimal_nside(l_box,1./I_R_MAX,cat_dat1->np);

#ifdef _CUTE_AS_PYTHON_MODULE
  if(global_galaxy_catalog2 == NULL){
    cat_dat2=read_catalog(fnameData2,&n_dat2);
  } else {
    cat_dat2=global_galaxy_catalog2;
    n_dat2=global_galaxy_catalog2->np;
  }
#else
  cat_dat2=read_catalog(fnameData2,&n_dat2);
#endif
  nside2=optimal_nside(l_box,1./I_R_MAX,cat_dat2->np);

  nside=MAX(nside1,nside2);
  data1_boxes=catalog_to_boxes(nside,*cat_dat1);
  data2_boxes=catalog_to_boxes(nside,*cat_dat2);

  if(use_randoms) {
    //Read randoms data
#ifdef _CUTE_AS_PYTHON_MODULE
    if(global_random_catalog == NULL){
      cat_rand=read_catalog(fnameRand,&n_rand);
    } else {
      cat_rand=global_random_catalog;
      n_rand=global_random_catalog->np;
    }
#else
    cat_rand=read_catalog(fnameRand,&n_rand);
#endif
    rand_boxes=catalog_to_boxes(nside,*cat_rand);
  }

#ifdef _DEBUG
  write_cat(cat_dat1,"debug_Dat1Cat.dat");
  write_cat(cat_dat1,"debug_Dat2Cat.dat");
  if(use_randoms) write_cat(cat_rand,"debug_RandCat.dat");
#endif

  printf("*** Correlating\n");
  timer(0);
  cross_3d_ps_boxes(nside,data1_boxes,data2_boxes,D1D2);
  if(use_randoms) {
    cross_3d_ps_boxes(nside,data1_boxes,rand_boxes,D1R);
    if(reuse_randoms==2) {
      for(i=0;i<NB_R*NB_R;i++) {    // will reuse pre-calculated D2R and RR to save time!
        D2R[i]=0;
        RR[i]=0;
        corr[i]=0;
        ercorr[i]=0;
      }
      timer(1);
      printf("\n");
      printf("*** Writing output without D2R and RR - adjust output files with pre-calculated values\n");
    }
    else if(reuse_randoms==1) {
      cross_3d_ps_boxes(nside,data2_boxes,rand_boxes,D2R);
      for(i=0;i<NB_R*NB_R;i++) {    // will reuse pre-calculated D2R and RR to save time!
        RR[i]=0;
        corr[i]=0;
        ercorr[i]=0;
      }
      timer(1);
      printf("\n");
      printf("*** Writing output without RR - adjust output files with pre-calculated values\n");
    }  else {
      cross_3d_ps_boxes(nside,data2_boxes,rand_boxes,D2R);
      auto_3d_ps_boxes(nside,rand_boxes,RR);
      timer(1);
      printf("\n");
      printf("*** Writing output\n");
      make_3d_CCF_w_rand(D1D2,D1R,D2R,RR,n_dat1,n_dat2,n_rand,corr,ercorr);
    }
    write_3d_CCF_w_rand(fnameOut,corr,ercorr,D1D2,D1R,D2R,RR);
  }
  else {    // ie, not using randoms
    timer(1);
    printf("\n");
    printf("*** Writing output\n");
    make_3d_CCF(D1D2,n_dat1,n_dat2,corr,ercorr);
    write_3d_CCF(fnameOut,corr,ercorr,D1D2);
  }

  printf("*** Cleaning up\n");
#ifdef _CUTE_AS_PYTHON_MODULE
  if(global_galaxy_catalog == NULL)
#endif
    free_catalog(cat_dat1);
#ifdef _CUTE_AS_PYTHON_MODULE
  if(global_galaxy_catalog2 == NULL)
#endif
    free_catalog(cat_dat2);
  free_boxes(nside,data1_boxes);
  free_boxes(nside,data2_boxes);
  if(use_randoms) {
#ifdef _CUTE_AS_PYTHON_MODULE
    if(global_random_catalog == NULL)
#endif
      free_catalog(cat_rand);
    free_boxes(nside,rand_boxes);
  }
  printf("\n");

  timer(5);
}

void run_3d_rmu_auto_corr_boxes(int use_randoms)
{
  //////
  // Runs xi(r,mu) in brute-force mode
  lint n_dat,n_rand;
  int nside;
  Catalog *cat_dat, *cat_rand;
  NeighborBox *data_boxes, *rand_boxes;
  unsigned long long DD[NB_R*NB_mu], DR[NB_R*NB_mu], RR[NB_R*NB_mu];
  double corr[NB_R*NB_mu],ercorr[NB_R*NB_mu];

  timer(4);

#ifdef _VERBOSE
  printf("*** 3D correlation function (r,mu): \n");
  printf(" - Range: (%.3lf,%.3lf) < (r,mu) < (%.3lf,%.3lf) Mpc/h\n",
      0.,0.,1./I_R_MAX,1.);
  printf(" - #bins: (%d,%d)\n",NB_R,NB_mu);
#ifdef _LOGBIN
  printf(" - Using logarithmic r-binning with %d bins per decade (log binning not recommended for (r,mu)!)\n",N_LOGINT);
#else //_LOGBIN
  printf(" - Resolution: Dr = %.3lf (Mpc/h)\n",1./(I_R_MAX*NB_R));
#endif //_LOGBIN
  printf(" - Using a brute-force approach \n");
  printf("\n");
#endif // _VERBOSE

  //Read data
#ifdef _CUTE_AS_PYTHON_MODULE
  if(global_galaxy_catalog == NULL){
    cat_dat=read_catalog(fnameData,&n_dat);
  } else {
    cat_dat=global_galaxy_catalog;
    n_dat=global_galaxy_catalog->np;
  }
#else
  cat_dat=read_catalog(fnameData,&n_dat);
#endif
  nside=optimal_nside(l_box,1./I_R_MAX,cat_dat->np);
  data_boxes=catalog_to_boxes(nside,*cat_dat);

  if(use_randoms) {
    //Read randoms data
#ifdef _CUTE_AS_PYTHON_MODULE
    if(global_random_catalog == NULL){
      cat_rand=read_catalog(fnameRand,&n_rand);
    } else {
      cat_rand=global_random_catalog;
      n_rand=global_random_catalog->np;
    }
#else
    cat_rand=read_catalog(fnameRand,&n_rand);
#endif
    rand_boxes=catalog_to_boxes(nside,*cat_rand);
  }

#ifdef _DEBUG
  write_cat(cat_dat,"debug_DatCat.dat");
  if(use_randoms) write_cat(cat_rand,"debug_RandCat.dat");
#endif

  printf("*** Correlating \n");
  printf(" - Auto-correlating data \n");
  timer(0);
  auto_3d_rmu_boxes(nside,data_boxes,DD);
  if(use_randoms) {
    timer(2);
    printf(" - Auto-correlating random \n");
    auto_3d_rmu_boxes(nside,rand_boxes,RR);
    timer(2);
    printf(" - Cross-correlating \n");
    cross_3d_rmu_boxes(nside,data_boxes,rand_boxes,DR);
    timer(1);

    printf("*** Writing output\n");
    make_3d_CF_w_rand(DD,DR,RR,n_dat,n_rand,corr,ercorr);
    write_3d_CF_w_rand(fnameOut,corr,ercorr,DD,DR,RR);
  }
  else {
    timer(1);
    printf("*** Writing output\n");
    make_3d_CF(DD,n_dat,corr,ercorr);
    write_3d_CF(fnameOut,corr,ercorr,DD);
  }

  printf("*** Cleaning up\n");
#ifdef _CUTE_AS_PYTHON_MODULE
  if(global_galaxy_catalog == NULL)
#endif
    free_catalog(cat_dat);
  free_boxes(nside,data_boxes);
  if(use_randoms) {
#ifdef _CUTE_AS_PYTHON_MODULE
    if(global_random_catalog == NULL)
#endif
      free_catalog(cat_rand);
    free_boxes(nside,rand_boxes);
  }
  printf("\n");

  timer(5);
}

void run_3d_rmu_cross_corr_boxes(int use_randoms, int reuse_randoms)
{
  //////
  // Runs xi(r,mu) cross-correlation in brute-force mode
  lint n_dat1,n_dat2,n_rand;
  int nside,nside1,nside2,i;
  Catalog *cat_dat1, *cat_dat2, *cat_rand;
  NeighborBox *data1_boxes, *data2_boxes, *rand_boxes;
  unsigned long long D1D2[NB_R*NB_mu], D1R[NB_R*NB_mu], D2R[NB_R*NB_mu], RR[NB_R*NB_mu];
  double corr[NB_R*NB_mu],ercorr[NB_R*NB_mu];

  timer(4);

#ifdef _VERBOSE
  printf("*** 3D cross-correlation function (r,mu): \n");
  printf(" - Range: (%.3lf,%.3lf) < (r,mu) < (%.3lf,%.3lf) Mpc/h\n",
      0.,0.,1./I_R_MAX,1.);
  printf(" - #bins: (%d,%d)\n",NB_R,NB_mu);
#ifdef _LOGBIN
  printf(" - Using logarithmic r-binning with %d bins per decade (log binning not recommended for (r,mu)!)\n",N_LOGINT);
#else //_LOGBIN
  printf(" - Resolution: Dr = %.3lf (Mpc/h)\n",1./(I_R_MAX*NB_R));
#endif //_LOGBIN
  printf(" - Using a brute-force approach \n");
  printf("\n");
#endif // _VERBOSE

  //Read data
#ifdef _CUTE_AS_PYTHON_MODULE
  if(global_galaxy_catalog == NULL){
    cat_dat1=read_catalog(fnameData,&n_dat1);
  } else {
    cat_dat1=global_galaxy_catalog;
    n_dat1=global_galaxy_catalog->np;
  }
#else
  cat_dat1=read_catalog(fnameData,&n_dat1);
#endif
  nside1=optimal_nside(l_box,1./I_R_MAX,cat_dat1->np);

#ifdef _CUTE_AS_PYTHON_MODULE
  if(global_galaxy_catalog2 == NULL){
    cat_dat2=read_catalog(fnameData2,&n_dat2);
  } else {
    cat_dat2=global_galaxy_catalog2;
    n_dat2=global_galaxy_catalog2->np;
  }
#else
  cat_dat2=read_catalog(fnameData2,&n_dat2);
#endif
  nside2=optimal_nside(l_box,1./I_R_MAX,cat_dat2->np);
  nside=MAX(nside1,nside2);
  data1_boxes=catalog_to_boxes(nside,*cat_dat1);
  data2_boxes=catalog_to_boxes(nside,*cat_dat2);

  if(use_randoms) {
    //Read randoms data
#ifdef _CUTE_AS_PYTHON_MODULE
    if(global_random_catalog == NULL){
      cat_rand=read_catalog(fnameRand,&n_rand);
    } else {
      cat_rand=global_random_catalog;
      n_rand=global_random_catalog->np;
    }
#else
    cat_rand=read_catalog(fnameRand,&n_rand);
#endif
    rand_boxes=catalog_to_boxes(nside,*cat_rand);
  }

#ifdef _DEBUG
  write_cat(cat_dat1,"debug_Dat1Cat.dat");
  write_cat(cat_dat1,"debug_Dat2Cat.dat");
  if(use_randoms) write_cat(cat_rand,"debug_RandCat.dat");
#endif

  printf("*** Correlating\n");
  timer(0);
  cross_3d_rmu_boxes(nside,data1_boxes,data2_boxes,D1D2);
  if(use_randoms) {
    cross_3d_rmu_boxes(nside,data1_boxes,rand_boxes,D1R);
    if(reuse_randoms==2) {
      for(i=0;i<NB_R*NB_mu;i++) {    // will reuse pre-calculated D2R and RR to save time!
        D2R[i]=0;
        RR[i]=0;
        corr[i]=0;
        ercorr[i]=0;
      }
      timer(1);
      printf("\n");
      printf("*** Writing output without D2R and RR - adjust output files with pre-calculated values \n");
    }
    else if(reuse_randoms==1) {
      cross_3d_rmu_boxes(nside,data2_boxes,rand_boxes,D2R);
      for(i=0;i<NB_R*NB_mu;i++) {    // will reuse pre-calculated D2R and RR to save time!
        RR[i]=0;
        corr[i]=0;
        ercorr[i]=0;
      }
      timer(1);
      printf("\n");
      printf("*** Writing output without RR - adjust output files with pre-calculated values \n");
    }
    else {
      cross_3d_rmu_boxes(nside,data2_boxes,rand_boxes,D2R);
      auto_3d_rmu_boxes(nside,rand_boxes,RR);
      timer(1);
      printf("\n");
      printf("*** Writing output\n");
      make_3d_CCF_w_rand(D1D2,D1R,D2R,RR,n_dat1,n_dat2,n_rand,corr,ercorr);
    }

    write_3d_CCF_w_rand(fnameOut,corr,ercorr,D1D2,D1R,D2R,RR);
  }
  else {    // ie, not using randoms
    timer(1);
    printf("\n");
    printf("*** Writing output\n");
    make_3d_CCF(D1D2,n_dat1,n_dat2,corr,ercorr);
    write_3d_CCF(fnameOut,corr,ercorr,D1D2);
  }

  printf("*** Cleaning up\n");
#ifdef _CUTE_AS_PYTHON_MODULE
  if(global_galaxy_catalog == NULL)
#endif
    free_catalog(cat_dat1);
#ifdef _CUTE_AS_PYTHON_MODULE
  if(global_galaxy_catalog2 == NULL)
#endif
    free_catalog(cat_dat2);
  free_boxes(nside,data1_boxes);
  free_boxes(nside,data2_boxes);
  if(use_randoms) {
#ifdef _CUTE_AS_PYTHON_MODULE
    if(global_random_catalog == NULL)
#endif
      free_catalog(cat_rand);
    free_boxes(nside,rand_boxes);
  }
  printf("\n");

  timer(5);
}

void run_monopole_corr_tree(void)
{
  //////
  // Main routine for monopole using tree
  lint n_dat;
  Catalog *cat_dat;
  branch *tree;
  unsigned long long DD[NB_R];
  double corr[NB_R],ercorr[NB_R];

  timer(4);

#ifdef _VERBOSE
  printf("*** Monopole: \n");
  printf(" - Range: %.3lf < r < %.3lf (Mpc/h)\n",0.,1./(I_R_MAX));
  printf(" - #bins: %d\n",NB_R);
#ifdef _LOGBIN
  printf(" - Using logarithmic binning with %d bins per decade \n",N_LOGINT);
#else
  printf(" - Resolution: Dr = %.3lf (Mpc/h)\n",1./(I_R_MAX*NB_R));
#endif
  printf(" - Using a tree algorithm \n");
  printf("\n");
#endif

  //Read data
  cat_dat=read_catalog(fnameData,&n_dat);
  tree=mk_tree(*cat_dat);
  printf("\n");

#ifdef _DEBUG
  write_cat(cat_dat,"debug_DatCat.dat");
  write_tree(tree,"debug_DatTree.dat");
  compute_tree_stats(tree,"debug_TreeStats.dat");
#endif //_DEBUG

  printf("*** Correlating\n");
  timer(0);
  corr_mono_box_tree(cat_dat->np,cat_dat->pos,tree,DD);
  timer(1);
  printf("\n");

  printf("*** Writing output \n");
  make_CF(DD,n_dat,corr,ercorr);
  write_CF(fnameOut,corr,ercorr,DD);

  printf("*** Cleaning up \n");
#ifdef _CUTE_AS_PYTHON_MODULE
  if(global_galaxy_catalog == NULL)
#endif
    free_catalog(cat_dat);
  free_branch(tree);

  timer(5);
}

void run_monopole_corr_pm(void)
{
  //////
  // Main routine for monopole using pm

  double *grid,*new_grid;
  int new_n_grid;
  //float *grid;
  /*lint n_dat;
    Catalog cat_dat;*/
  unsigned long long DD[NB_R];
  double corr[NB_R],ercorr[NB_R];
  timer(4);

#ifdef _VERBOSE
  printf("*** Correlation function parameters: \n");
  printf(" - Range: %.3lf < r < %.3lf (Mpc/h)\n",0.,1./(I_R_MAX));
  printf(" - #bins: %d\n",NB_R);
#ifdef _LOGBIN
  printf(" - Using logarithmic binning with %d bins per decade \n",N_LOGINT);
#else
  printf(" - Resolution: Dr = %.3lf (Mpc/h)\n",1./(I_R_MAX*NB_R));
#endif
  printf(" - Using a PM approach\n");
  printf("\n");
#endif


  /*  //Read data
      cat_dat=read_catalog(fnameData,&n_dat);
      printf("*** Calculating PM grid \n");
      grid=pos_2_tsc(cat_dat);
      printf("\n");
#ifdef _DEBUG
write_cat(cat_dat,"debug_DatCat.dat");
write_grid(grid,"debug_DatGrid.dat");
#endif
free_catalog(cat_dat);
   */
  grid = read_grid();
  new_n_grid = 146; // (int)(l_box*I_R_MAX*NB_R);
  if(new_n_grid<n_grid) {
    printf("  Appropriate grid size: %d\n",new_n_grid);
    new_grid = resize_grid(grid,new_n_grid);

    n_grid = new_n_grid;
    grid = new_grid;
  }
#ifdef _DEBUG
  write_grid(grid,"debug_DatGrid.dat");
#endif

  printf("*** Correlating\n");
  timer(0);
  corr_mono_box_pm(grid,corr,ercorr,DD);
  timer(1);
  printf("\n");

  write_CF(fnameOut,corr,ercorr,DD);

  printf("*** Cleaning up \n");
  free(grid);
  printf("\n");

  timer(5);

}

void run_monopole_corr_p3m(void)
{
  //////
}

#ifdef _CUTE_AS_PYTHON_MODULE

int runCUTEbox(Catalog *galaxy_catalog, Catalog *galaxy_catalog2, Catalog *random_catalog, Result *result, int verbose){
  int ii;
  cute_verbose = verbose;

#else


  int main(int argc,char **argv){
    //////
    // Main routine
    int ii;
    char fnameIn[128];
    if(argc!=2) {
      printf("Usage ./CUTE_box <input file>\n");
      exit(1);
    }
    sprintf(fnameIn,"%s",argv[1]);

#endif

    setbuf(stdout, NULL);

    printf("\n");
    printf("-----------------------------------------------------------\n");
    printf("|| CUTE - Correlation Utilities and Two-point Estimation ||\n");
    printf("-----------------------------------------------------------\n\n");

#ifdef _CUTE_AS_PYTHON_MODULE
  global_galaxy_catalog  = galaxy_catalog;
  global_galaxy_catalog2 = galaxy_catalog2;
  global_random_catalog  = random_catalog;
  
  // Initialize
  if(galaxy_catalog != NULL){
    printf("Using external data catalog with np = %d\n", 
        (int)galaxy_catalog->np);
  }
  if(galaxy_catalog2 != NULL){
    printf("Using second external data catalog with np = %d\n", 
        (int)galaxy_catalog2->np);
  }
  if(random_catalog != NULL){
    printf("Using external random catalog with np = %d\n", 
        (int)random_catalog->np);
  }

  global_result = result;
#endif

    //Initialize random number generator
#ifdef _DEBUG
    srand(1234);
#else
    srand(time(NULL));
#endif
#ifdef _VERBOSE
    printf("Initializing random number generator\n");
    printf("First random number : %d \n",rand());
#endif

#ifdef _VERBOSE
    //Calculate number of threads
    ii=0;
#pragma omp parallel
    {
#pragma omp atomic
      ii++;
    }
    printf("Using %d threads \n",ii);
#endif

    printf("\n");

#ifndef _CUTE_AS_PYTHON_MODULE  
    read_run_params(fnameIn);
#endif

    if(corr_type==1) {
      if(do_CCF)
        run_monopole_CCF(use_randoms,reuse_randoms);
      else if(use_pm==1)
        run_monopole_corr_pm();
      else if(use_pm==2)
        run_monopole_corr_p3m();
      else {
        if(use_tree)
          run_monopole_corr_tree();
        else if(use_randoms){
          run_monopole_corr_neighbors_w_rand();
        }
        else{
          run_monopole_corr_neighbors();
        }
      }
    }
    else if(corr_type==2) {
      if(do_CCF)
        run_3d_ps_cross_corr_boxes(use_randoms,reuse_randoms);
      else
        run_3d_ps_auto_corr_boxes(use_randoms);
    }
    else if(corr_type==3) {
      if(do_CCF)
        run_3d_rmu_cross_corr_boxes(use_randoms,reuse_randoms);
      else {
        run_3d_rmu_auto_corr_boxes(use_randoms);
      }
    }

    printf("             Done !!!             \n\n");

    return 0;
  }
