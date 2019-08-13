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
//          Read-write routines for the periodic-box mode            //
/*********************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "define.h"
#include "common.h"

#ifdef _CUTE_AS_PYTHON_MODULE
static int param_errors = 0;
#endif

typedef struct {
  int npart[6];
  double mass[6];
  double time;
  double redshift;
  int flag_sfr;
  int flag_feedback;
  int npartTotal[6];
  int flag_cooling;
  int num_files;
  double BoxSize;
  double Omega0;
  double OmegaLambda;
  double HubbleParam;
  char fill[256-6*4-6*8-2*8-2*4-6*4-2*4-4*8];
  // fills to 256 Bytes
} gad_header;

typedef struct {
  char label[4];
  int size;
} gad_title;

static double wrap_double(double x)
{
  //////
  // Returns x mod(l_box)
  if(x<0) 
    return wrap_double(x+l_box);
  else if(x>=l_box) 
    return wrap_double(x-l_box);
  else 
    return x;
}

static size_t my_fread(void *p,size_t size,
    size_t nmemb,FILE *stream)
{
  //////
  // Self-checked binary reading routine
  size_t nread;

  if((nread=fread(p,size,nmemb,stream))!=nmemb) {
    fprintf(stderr,"CUTE: error reading binary file \n");
    exit(1);
  }
  return nread;
}

void make_CF_double(unsigned long long DD[],int nD,
    double corr[],double ercorr[])
{
  //////
  // Creates correlation function and poisson errors
  // from pair counts DD
  double *edd;
  double rho_av=nD/(l_box*l_box*l_box);
  int ii;

  edd=(double *)malloc(sizeof(double)*NB_R);
  if(edd==NULL)
    error_mem_out();

#ifndef _LOGBIN
  DD[0]-=nD; //Subtract diagonal (self pair counts)
#endif //_LOGBIN

  for(ii=0;ii<NB_R;ii++)
    edd[ii]=1./sqrt((double)DD[ii]);

  for(ii=0;ii<NB_R;ii++) {
    if(DD[ii]==0) {
      corr[ii]=0;
      ercorr[ii]=0;
    }
    else {
      double r0,r1,vr,rho_r;
#ifdef _LOGBIN
      r0=pow(10.,(((double)ii-NB_R)/N_LOGINT)+LOG_R_MAX);
      r1=pow(10.,(((double)ii+1-NB_R)/N_LOGINT)+LOG_R_MAX);
#else //_LOGBIN
      r0=ii/(I_R_MAX*NB_R);
      r1=(ii+1)/(I_R_MAX*NB_R);
#endif //_LOGBIN
      vr=4*M_PI*(r1*r1*r1-r0*r0*r0)/3;
      rho_r=DD[ii]/(nD*vr);
      corr[ii]=rho_r/rho_av-1;
      ercorr[ii]=(1+corr[ii])*edd[ii];
    }
  }

  free(edd);
}

void make_CF(unsigned long long DD[],int nD,
    double corr[],double ercorr[])
{
  //////
  // Creates correlation function and poisson errors
  // from pair counts DD
  double *edd;
  double rho_av=nD/(l_box*l_box*l_box);
  int ii;

  edd=(double *)malloc(sizeof(double)*NB_R);
  if(edd==NULL)
    error_mem_out();

  for(ii=0;ii<NB_R;ii++)
    edd[ii]=1./sqrt((double)DD[ii]*2); // multiply by 2 because DD counted each pair only once

  for(ii=0;ii<NB_R;ii++) {
    if(DD[ii]==0) {
      corr[ii]=0;
      ercorr[ii]=0;
    }
    else {
      double r0,r1,vr,rho_r;
#ifdef _LOGBIN
      r0=pow(10.,(((double)ii-NB_R)/N_LOGINT)+LOG_R_MAX);
      r1=pow(10.,(((double)ii+1-NB_R)/N_LOGINT)+LOG_R_MAX);
#else //_LOGBIN
      r0=ii/(I_R_MAX*NB_R);
      r1=(ii+1)/(I_R_MAX*NB_R);
#endif //_LOGBIN
      vr=4*M_PI*(r1*r1*r1-r0*r0*r0)/3;
      rho_r=2.*DD[ii]/(nD*vr);	// multiply by 2 because DD counted each pair only once
      corr[ii]=rho_r/rho_av-1;
      ercorr[ii]=(1+corr[ii])*edd[ii];
    }
  }

  free(edd);
}

void make_3d_CCF_w_rand(unsigned long long D1D2[],unsigned long long D1R[],unsigned long long D2R[],
    unsigned long long RR[],int nD1,int nD2,int nR,double corr[],double ercorr[])
{
  //////
  // Creates cross-correlation function
  // from pair counts D1D2, D1R, D2R and RR
  // (Error estimate may be nonsensical, don't trust it)
  double *ddd,*dd1r,*dd2r,*drr;
  double norm_dd=1.0*((double)nD1)*nD2; 	
  double norm_d1r=1.0*((double)nD1)*nR;		
  double norm_d2r=1.0*((double)nD2)*nR;		
  double norm_rr=0.5*((double)nR*nR-nR);	
  double fact;
  int ii,jj,index;

  if(corr_type==1) {
    ddd=(double *)malloc(sizeof(double)*NB_R);
    if(ddd==NULL)
      error_mem_out();
    dd1r=(double *)malloc(sizeof(double)*NB_R);
    if(dd1r==NULL)
      error_mem_out();
    dd2r=(double *)malloc(sizeof(double)*NB_R);
    if(dd2r==NULL)
      error_mem_out();
    drr=(double *)malloc(sizeof(double)*NB_R);
    if(drr==NULL)
      error_mem_out();

    for(ii=0;ii<NB_R;ii++) {
      ddd[ii]=(double)(D1D2[ii]/norm_dd);
      dd1r[ii]=(double)(D1R[ii]/norm_d1r);
      dd2r[ii]=(double)(D2R[ii]/norm_d2r);
      drr[ii]=(double)(RR[ii]/norm_rr);
      if(drr[ii]==0) {
        corr[ii]=0;
        ercorr[ii]=0;
      }
      else {
        corr[ii]=ddd[ii]/drr[ii]-dd1r[ii]/drr[ii]-dd2r[ii]/drr[ii]+1.0;
        fact = norm_dd/norm_rr*RR[ii]*(1+corr[ii]) + 4.*(drr[ii]*norm_dd*(1+corr[ii]))*(drr[ii]*norm_dd*(1+corr[ii]))/sqrt(nD1*nD2);
        ercorr[ii]=sqrt(3.)*(norm_rr/(norm_dd*RR[ii]))*sqrt(fact);
      }
    }
    free(ddd);
    free(dd1r);
    free(dd2r);
    free(drr);
  }
  else if(corr_type==2) {
    ddd=(double *)malloc(sizeof(double)*NB_R*NB_R);
    if(ddd==NULL)
      error_mem_out();
    dd1r=(double *)malloc(sizeof(double)*NB_R*NB_R);
    if(dd1r==NULL)
      error_mem_out();
    dd2r=(double *)malloc(sizeof(double)*NB_R*NB_R);
    if(dd2r==NULL)
      error_mem_out();
    drr=(double *)malloc(sizeof(double)*NB_R*NB_R);
    if(drr==NULL)
      error_mem_out();

    for(ii=0;ii<NB_R;ii++) {
      for(jj=0;jj<NB_R;jj++) {
        index = jj+ii*NB_R;
        ddd[index]=(double)(D1D2[index]/norm_dd);
        dd1r[index]=(double)(D1R[index]/norm_d1r);
        dd2r[index]=(double)(D2R[index]/norm_d2r);
        drr[index]=(double)(RR[index]/norm_rr);
        if(drr[index]==0) {
          corr[index]=0;
          ercorr[index]=0;
        }
        else {
          corr[index]=ddd[index]/drr[index]-dd1r[index]/drr[index]-dd2r[index]/drr[index]+1.0;
          fact = norm_dd/norm_rr*RR[index]*(1+corr[index]) + 
            4.*(drr[index]*norm_dd*(1+corr[index]))*(drr[index]*norm_dd*(1+corr[index]))/sqrt(nD1*nD2);
          ercorr[index]=sqrt(3.)*(norm_rr/(norm_dd*RR[index]))*sqrt(fact);
        }
      }
    }
    free(ddd);
    free(dd1r);
    free(dd2r);
    free(drr);
  }
  else if(corr_type==3) {
    ddd=(double *)malloc(sizeof(double)*NB_R*NB_mu);
    if(ddd==NULL)
      error_mem_out();
    dd1r=(double *)malloc(sizeof(double)*NB_R*NB_mu);
    if(dd1r==NULL)
      error_mem_out();
    dd2r=(double *)malloc(sizeof(double)*NB_R*NB_mu);
    if(dd2r==NULL)
      error_mem_out();
    drr=(double *)malloc(sizeof(double)*NB_R*NB_mu);
    if(drr==NULL)
      error_mem_out();

    for(ii=0;ii<NB_R;ii++) {
      for(jj=0;jj<NB_mu;jj++) {
        index = jj+ii*NB_mu;
        ddd[index]=(double)(D1D2[index]/norm_dd);
        dd1r[index]=(double)(D1R[index]/norm_d1r);
        dd2r[index]=(double)(D2R[index]/norm_d2r);
        drr[index]=(double)(RR[index]/norm_rr);
        if(drr[index]==0) {
          corr[index]=0;
          ercorr[index]=0;
        }
        else {
          corr[index]=ddd[index]/drr[index]-dd1r[index]/drr[index]-dd2r[index]/drr[index]+1.0;
          fact = norm_dd/norm_rr*RR[index]*(1+corr[index]) + 
            4.*(drr[index]*norm_dd*(1+corr[index]))*(drr[index]*norm_dd*(1+corr[index]))/sqrt(nD1*nD2);
          ercorr[index]=sqrt(3.)*(norm_rr/(norm_dd*RR[index]))*sqrt(fact);
        }
      }
    }
    free(ddd);
    free(dd1r);
    free(dd2r);
    free(drr);
  }
}

void make_3d_CF_w_rand(unsigned long long DD[],unsigned long long DR[],
    unsigned long long RR[],int nD,int nR,double corr[],double ercorr[])
{
  //////
  // Creates correlation function and poisson errors
  // from pair counts DD, DR and RR
  double *ddd,*ddr,*drr;
  double norm_dd=0.5*((double)nD*nD-nD); 	
  double norm_rr=0.5*((double)nR*nR-nR);	
  double norm_dr=1.0*((double)nD)*nR;		
  double fact;
  int ii,jj,index;

  if(corr_type==1) {
    ddd=(double *)malloc(sizeof(double)*NB_R);
    if(ddd==NULL)
      error_mem_out();
    ddr=(double *)malloc(sizeof(double)*NB_R);
    if(ddr==NULL)
      error_mem_out();
    drr=(double *)malloc(sizeof(double)*NB_R);
    if(drr==NULL)
      error_mem_out();

    for(ii=0;ii<NB_R;ii++) {
      ddd[ii]=(double)(DD[ii]/norm_dd);
      ddr[ii]=(double)(DR[ii]/norm_dr);
      drr[ii]=(double)(RR[ii]/norm_rr);
      if(drr[ii]==0) {
        corr[ii]=0;
        ercorr[ii]=0;
      }
      else {
        corr[ii]=ddd[ii]/drr[ii]-2*ddr[ii]/drr[ii]+1.0;
        fact = norm_dd/norm_rr*RR[ii]*(1+corr[ii]) + 4.*(drr[ii]*norm_dd*(1+corr[ii]))*(drr[ii]*norm_dd*(1+corr[ii]))/nD;
        ercorr[ii]=sqrt(3.)*(norm_rr/(norm_dd*RR[ii]))*sqrt(fact);
      }
    }
    free(ddd);
    free(ddr);
    free(drr);
  }
  else if(corr_type==2) {
    ddd=(double *)malloc(sizeof(double)*NB_R*NB_R);
    if(ddd==NULL)
      error_mem_out();
    ddr=(double *)malloc(sizeof(double)*NB_R*NB_R);
    if(ddr==NULL)
      error_mem_out();
    drr=(double *)malloc(sizeof(double)*NB_R*NB_R);
    if(drr==NULL)
      error_mem_out();

    for(ii=0;ii<NB_R;ii++) {
      for(jj=0;jj<NB_R;jj++) {
        index = ii+jj*NB_R;
        ddd[index]=(double)(DD[index]/norm_dd);
        ddr[index]=(double)(DR[index]/norm_dr);
        drr[index]=(double)(RR[index]/norm_rr);
        if(drr[index]==0) {
          corr[index]=0;
          ercorr[index]=0;
        }
        else {
          corr[index]=ddd[index]/drr[index]-2*ddr[index]/drr[index]+1.0;
          fact = norm_dd/norm_rr*RR[index]*(1+corr[index]) + 
            4.*(drr[index]*norm_dd*(1+corr[index]))*(drr[index]*norm_dd*(1+corr[index]))/nD;
          ercorr[index]=sqrt(3.)*(norm_rr/(norm_dd*RR[index]))*sqrt(fact);
        }
      }
    }
    free(ddd);
    free(ddr);
    free(drr);
  }
  else if(corr_type==3) {
    ddd=(double *)malloc(sizeof(double)*NB_R*NB_mu);
    if(ddd==NULL)
      error_mem_out();
    ddr=(double *)malloc(sizeof(double)*NB_R*NB_mu);
    if(ddr==NULL)
      error_mem_out();
    drr=(double *)malloc(sizeof(double)*NB_R*NB_mu);
    if(drr==NULL)
      error_mem_out();

    for(ii=0;ii<NB_R;ii++) {
      for(jj=0;jj<NB_mu;jj++) {
        index = jj+ii*NB_mu;
        ddd[index]=(double)(DD[index]/norm_dd);
        ddr[index]=(double)(DR[index]/norm_dr);
        drr[index]=(double)(RR[index]/norm_rr);
        if(drr[index]==0) {
          corr[index]=0;
          ercorr[index]=0;
        }
        else {
          corr[index]=ddd[index]/drr[index]-2*ddr[index]/drr[index]+1.0;
          fact = norm_dd/norm_rr*RR[index]*(1+corr[index]) + 
            4.*(drr[index]*norm_dd*(1+corr[index]))*(drr[index]*norm_dd*(1+corr[index]))/nD;
          ercorr[index]=sqrt(3.)*(norm_rr/(norm_dd*RR[index]))*sqrt(fact);
        }
      }
    }
    free(ddd);
    free(ddr);
    free(drr);
  } 
}

void make_3d_CCF(unsigned long long D1D2[],int nD1,int nD2,double corr[],double ercorr[])
{
  //////
  // Creates cross-correlation function from pair counts D1D2
  // assuming homogeneous selection function (no randoms)
  // (Error estimate may be nonsensical, don't trust it)
  double *ddd, *vols;
  double norm_dd=1.0*((double)nD1)*nD2;
  double inv_box_vol=1.0/(l_box*l_box*l_box);
  double r0,r1;
  int ii,jj,index;

  if(corr_type==1) {
    ddd=(double *)malloc(sizeof(double)*NB_R);
    if(ddd==NULL)
      error_mem_out();
    vols=(double *)malloc(sizeof(double)*NB_R);
    if(vols==NULL)
      error_mem_out();

    for(ii=0;ii<NB_R;ii++) {
      ddd[ii]=(double)(D1D2[ii]/norm_dd);
#ifdef _LOGBIN
      r0=pow(10.,(((double)ii-NB_R)/N_LOGINT)+LOG_R_MAX);
      r1=pow(10.,(((double)ii+1-NB_R)/N_LOGINT)+LOG_R_MAX);
#else 
      r0=ii/(I_R_MAX*NB_R);
      r1=(ii+1)/(I_R_MAX*NB_R);
#endif 
      vols[ii]=4*M_PI*(r1*r1*r1-r0*r0*r0)/3;
      corr[ii]=ddd[ii]/(vols[ii]*inv_box_vol)-1.0;
      ercorr[ii]=(1+corr[ii])/sqrt(D1D2[ii]);
    }
    free(ddd);
    free(vols);
  }
  else if(corr_type==2) {
    ddd=(double *)malloc(sizeof(double)*NB_R*NB_R);
    if(ddd==NULL)
      error_mem_out();
    vols=(double *)malloc(sizeof(double)*NB_R*NB_R);
    if(vols==NULL)
      error_mem_out();

    for(ii=0;ii<NB_R;ii++) {
      for(jj=0;jj<NB_R;jj++) {
        index = jj+ii*NB_R;
#ifdef _LOGBIN
        r0=pow(10.,(((double)jj-NB_R)/N_LOGINT)+LOG_R_MAX);
        r1=pow(10.,(((double)jj+1-NB_R)/N_LOGINT)+LOG_R_MAX);
#else 
        r0=jj/(I_R_MAX*NB_R);
        r1=(jj+1)/(I_R_MAX*NB_R);
#endif 
        vols[index]=2*M_PI*(r1+r0)*(r1-r0)*(r1-r0);
        ddd[index]=(double)(D1D2[index]/norm_dd);
        corr[index]=ddd[index]/(vols[index]*inv_box_vol)-1.0;
        ercorr[index]=(1+corr[index])/sqrt(D1D2[index]);
      }
    }
    free(ddd);
    free(vols);
  }
  else if(corr_type==3) {
    ddd=(double *)malloc(sizeof(double)*NB_R*NB_mu);
    if(ddd==NULL)
      error_mem_out();
    vols=(double *)malloc(sizeof(double)*NB_R*NB_mu);
    if(vols==NULL)
      error_mem_out();

    for(ii=0;ii<NB_R;ii++) {
      for(jj=0;jj<NB_mu;jj++) {
        index = jj+ii*NB_mu;
#ifdef _LOGBIN
        r0=pow(10.,(((double)ii-NB_R)/N_LOGINT)+LOG_R_MAX);
        r1=pow(10.,(((double)ii+1-NB_R)/N_LOGINT)+LOG_R_MAX);
#else 
        r0=ii/(I_R_MAX*NB_R);
        r1=(ii+1)/(I_R_MAX*NB_R);
#endif
        //vols[index]=M_PI*(r1+r0)*(r1+r0)*(r1-r0)/NB_mu;
        vols[index]=4.*M_PI*(r1*r1*r1-r0*r0*r0)/(3.0*NB_mu);
        ddd[index]=(double)(D1D2[index]/norm_dd);
        corr[index]=ddd[index]/(vols[index]*inv_box_vol)-1.0;
        ercorr[index]=(1+corr[index])/sqrt(D1D2[index]);
      }
    }
    free(ddd);
    free(vols);
  }
}

void make_3d_CF(unsigned long long DD[],int nD,double corr[],double ercorr[])
{
  //////
  // Creates correlation function and poisson errors
  // from pair counts DD, DR and RR
  double *ddd,*vols;
  double inv_box_vol=1.0/(l_box*l_box*l_box);
  int ii,jj,index;
  double r0,r1;

  fprintf(stderr,"DD[0]=%llu \n",DD[0]);

  if(corr_type==1) {
    ddd=(double *)malloc(sizeof(double)*NB_R);
    if(ddd==NULL)
      error_mem_out();
    vols=(double *)malloc(sizeof(double)*NB_R);
    if(vols==NULL)
      error_mem_out();

    for(ii=0;ii<NB_R;ii++) {
      ddd[ii]=(double)(DD[ii]/(1.0*((double)nD)*nD));
#ifdef _LOGBIN
      r0=pow(10.,(((double)ii-NB_R)/N_LOGINT)+LOG_R_MAX);
      r1=pow(10.,(((double)ii+1-NB_R)/N_LOGINT)+LOG_R_MAX);
#else //_LOGBIN
      r0=ii/(I_R_MAX*NB_R);
      r1=(ii+1)/(I_R_MAX*NB_R);
#endif //_LOGBIN
      vols[ii]=4*M_PI*(r1*r1*r1-r0*r0*r0)/3;
      corr[ii]=2*ddd[ii]/(vols[ii]*inv_box_vol)-1.0;
      ercorr[ii]=(1+corr[ii])/sqrt(DD[ii]);
    }
    free(ddd);
    free(vols);
  }
  else if(corr_type==2) {
    ddd=(double *)malloc(sizeof(double)*NB_R*NB_R);
    if(ddd==NULL)
      error_mem_out();
    vols=(double *)malloc(sizeof(double)*NB_R*NB_R);
    if(vols==NULL)
      error_mem_out();

    for(ii=0;ii<NB_R;ii++) {
      for(jj=0;jj<NB_R;jj++) {
        index = ii+jj*NB_R;
        ddd[index]=(double)(DD[index]/(1.0*((double)nD)*nD));
#ifdef _LOGBIN
        r0=pow(10.,(((double)jj-NB_R)/N_LOGINT)+LOG_R_MAX);
        r1=pow(10.,(((double)jj+1-NB_R)/N_LOGINT)+LOG_R_MAX);
#else //_LOGBIN
        r0=jj/(I_R_MAX*NB_R);
        r1=(jj+1)/(I_R_MAX*NB_R);
#endif //_LOGBIN
        vols[index]=2*M_PI*(r1+r0)*(r1-r0)*(r1-r0);
        corr[index]=2*ddd[index]/(vols[index]*inv_box_vol)-1.0;
        ercorr[index]=(1+corr[index])/sqrt(DD[index]);
      }
    }
    free(ddd);
    free(vols);
  }
  else if(corr_type==3) {
    ddd=(double *)malloc(sizeof(double)*NB_R*NB_mu);
    if(ddd==NULL)
      error_mem_out();
    vols=(double *)malloc(sizeof(double)*NB_R*NB_mu);
    if(vols==NULL)
      error_mem_out();

    for(ii=0;ii<NB_R;ii++) {
      for(jj=0;jj<NB_mu;jj++) {
        index = jj+ii*NB_mu;
        ddd[index]=(double)(DD[index]/(1.0*((double)nD)*nD));
#ifdef _LOGBIN
        r0=pow(10.,(((double)ii-NB_R)/N_LOGINT)+LOG_R_MAX);
        r1=pow(10.,(((double)ii+1-NB_R)/N_LOGINT)+LOG_R_MAX);
#else //_LOGBIN
        r0=ii/(I_R_MAX*NB_R);
        r1=(ii+1)/(I_R_MAX*NB_R);
#endif //_LOGBIN
        //vols[index]=M_PI*(r1+r0)*(r1+r0)*(r1-r0)/NB_mu;
        vols[index]=4.*M_PI*(r1*r1*r1-r0*r0*r0)/(3.0*NB_mu);
        corr[index]=2*ddd[index]/(vols[index]*inv_box_vol)-1.0;
        ercorr[index]=(1+corr[index])/sqrt(DD[index]);
      }
    }
    free(ddd);
    free(vols);
  }
}

static void check_params(void)
{
  //////
  // Check all parameters are there and are sensible

  //vital files
  if(!strcmp(fnameData,"default")) {
    fprintf(stderr,"CUTE: Data catalog was not provided \n");
#ifndef _CUTE_AS_PYTHON_MODULE   
    exit(1);
#else
    param_errors++;
#endif
  }
  if(!strcmp(fnameOut,"default")) {
    fprintf(stderr,"CUTE: Output filename was not provided \n");
#ifndef _CUTE_AS_PYTHON_MODULE   
    exit(1);
#else
    param_errors++;
#endif
  }
  if(l_box<0) {
    fprintf(stderr,"CUTE: Box size was not provided \n");
#ifndef _CUTE_AS_PYTHON_MODULE   
    exit(1);
#else
    param_errors++;
#endif
  }
  //input format
  if((input_format!=0)&&(input_format!=1)&&(input_format!=2)) {
    fprintf(stderr,"CUTE: wrong input format. Using standard ASCII file \n");
    input_format=0;
  }
  if((input_format>0)&&(n_objects!=-1)) {
    fprintf(stderr,"CUTE: can't select #objects for GADGET input format.");
    fprintf(stderr," Reading all objects \n");
    n_objects=-1;
  }
  //numbers of objects from data and random
  if((n_objects!=-1)&&(n_objects<0)) {
    fprintf(stderr,"CUTE: Wrong wumber of lines from the data file. ");
    fprintf(stderr,"CUTE: All objects in the catalog will be read.");
    n_objects=-1;
  }
  if((use_pm)&&(use_tree)) {
    fprintf(stderr,"CUTE: Tree and PM algorithms CAN'T be combined.");
    fprintf(stderr," Using tree \n");
    use_pm=0;
    use_tree=1;
  }

  if(use_tree<0) {
    fprintf(stderr,"CUTE: No Tree option was provided \n");
#ifndef _CUTE_AS_PYTHON_MODULE   
    exit(1);
#else
    param_errors++;
#endif
  }
  else if(use_tree) {
    if(max_tree_order<0) {
      fprintf(stderr,"CUTE: Maximum tree order was not provided. Using 7 \n");
      max_tree_order=7;
    }
    if(max_tree_nparts<0) {
      fprintf(stderr,"CUTE: Maximum #particles per leaf was not provided.");
      fprintf(stderr," Using 10\n");
      max_tree_nparts=10;
    }
#ifndef _LOGBIN //Check resolution
    double binsize=1/(I_R_MAX*NB_R);
    double leafsize=l_box/pow(2,max_tree_order);

    if(binsize<=leafsize) {
      fprintf(stderr,"CUTE: Warning! binsize is smaller than cell size (%.3lf < %.3lf). ",
          binsize,leafsize);
      fprintf(stderr," Tree may be suboptimal \n");
    }
#endif //_LOGBIN
  }

  if(!((corr_type==1)||(corr_type==2)||(corr_type==3))) {
    fprintf(stderr,"CUTE: Invalid correlation type selection in input corr_type\n");
    fprintf(stderr," Valid values are: \n\tcorr_type=1 (monopole, xi(r))\n\tcorr_type=2 (3d, xi(sigma,pi))\n\tcorr_type=3 (3d, xi(r,mu)) \n");
#ifndef _CUTE_AS_PYTHON_MODULE   
    exit(1);
#else
    param_errors++;
#endif
  }

  if(do_CCF) {
    if(!strcmp(fnameData2,"default")) {
      fprintf(stderr,"CUTE: CCF requested but second data catalog was not provided \n");
#ifndef _CUTE_AS_PYTHON_MODULE   
      exit(1);
#else
      param_errors++;
#endif
    }
    if(!strcmp(fnameRand,"default")) {
      fprintf(stderr,"CUTE: CCF requires random catalog, but this was not provided \n");
#ifndef _CUTE_AS_PYTHON_MODULE   
      exit(1);
#else
      param_errors++;
#endif
    }
  }
  if((use_randoms)&&(!strcmp(fnameRand,"default"))) {
    fprintf(stderr,"CUTE: use_randoms selected but random catalog was not provided, proceeding without randoms \n");
    use_randoms=0;
  }

  if((use_randoms)&&(reuse_randoms)&&(do_CCF)) {
    fprintf(stderr,"CUTE: only reporting D1D2 and D1R pair counts! CCF, D2R and RR will be set to zero\n");
  }

  if((use_pm!=0)&&(use_pm!=1)) {
    fprintf(stderr,"CUTE: No PM option was provided \n");
#ifndef _CUTE_AS_PYTHON_MODULE   
    exit(1);
#else
    param_errors++;
#endif
  }
  if(use_pm) {
    if(n_grid<=0) {
      fprintf(stderr,"CUTE: No PM grid size was provided \n");
#ifndef _CUTE_AS_PYTHON_MODULE   
      exit(1);
#else
      param_errors++;
#endif
    }

    //Check resolution
    double cellsize=l_box/n_grid;
#ifdef _LOGBIN
    fprintf(stderr,"CUTE: Warning! Logarithmic binning with PM algorithm. ");
    fprintf(stderr,"Scales below cellsize (%.3lf) won't be correctly calculated \n",cellsize);
#else //_LOGBIN
    double binsize=1/(I_R_MAX*NB_R);
    if(binsize<=cellsize) {
      fprintf(stderr,"CUTE: Warning! binsize is smaller than cell size (%.3lf < %.3lf). ",
          binsize,cellsize);
      fprintf(stderr," Using PM is not recommended \n");
    }
#endif //_LOGBIN
  }
}

#ifdef _CUTE_AS_PYTHON_MODULE
void print_parameters(){    
  printf("\n===================================\n");
  printf("CUTE Parameters: \n");
  printf("===================================\n");
  printf(" data_filename    = %s\n", fnameData);
  printf(" data_filename2   = %s\n", fnameData2);
  printf(" random_filename  = %s\n", fnameRand);
  printf(" use_randoms      = %i\n", use_randoms);
  printf(" reuse_randoms    = %i\n", reuse_randoms);
  printf(" num_lines        = %i\n", (int)n_objects);
  printf(" input_format     = %i\n", input_format);
  printf(" output_filename  = %s\n", fnameOut);
  printf(" corr_type        = %i\n", corr_type);
  printf(" use_pm           = %i\n", use_pm);
  printf(" use_tree         = %i\n", use_tree);
  printf(" box_size         = %f\n", l_box);
  printf(" max_tree_order   = %i\n", max_tree_order);
  printf(" max_tree_nparts  = %i\n", max_tree_nparts);
  printf(" do_CCF           = %i\n", do_CCF);
  printf(" n_grid_side      = %i\n", n_grid);
  printf("===================================\n\n");
}
#endif

void read_run_params(char *fname)
{
  //////
  // Reads and checks the parameter file
  FILE *fi;
  int n_lin,ii;

  printf("*** Reading run parameters \n");
  //Read parameters from file
  fi=fopen(fname,"r");
  if(fi==NULL) error_open_file(fname);
  n_lin=linecount(fi);
  rewind(fi);
  for(ii=0;ii<n_lin;ii++) {
    char s0[512],s1[64],s2[256];
    if(fgets(s0,sizeof(s0),fi)==NULL)
      error_read_line(fname,ii+1);
    if((s0[0]=='#')||(s0[0]=='\n')) continue;
    int sr=sscanf(s0,"%s %s",s1,s2);
    if(sr!=2)
      error_read_line(fname,ii+1);

    if(!strcmp(s1,"data_filename="))
      sprintf(fnameData,"%s",s2);
    else if(!strcmp(s1,"data2_filename="))
      sprintf(fnameData2,"%s",s2);
    else if(!strcmp(s1,"random_filename="))
      sprintf(fnameRand,"%s",s2);
    else if(!strcmp(s1,"num_lines=")) {
      if(!strcmp(s2,"all"))
        n_objects=-1;
      else
        n_objects=atoi(s2);
    }
    else if(!strcmp(s1,"input_format="))
      input_format=atoi(s2);
    else if(!strcmp(s1,"output_filename="))
      sprintf(fnameOut,"%s",s2);
    else if(!strcmp(s1,"box_size=")) {
      l_box=atof(s2);
      l_box_half=0.5*l_box;
    }
    else if(!strcmp(s1,"use_tree="))
      use_tree=atoi(s2);
    else if(!strcmp(s1,"max_tree_order="))
      max_tree_order=atoi(s2);
    else if(!strcmp(s1,"max_tree_nparts="))
      max_tree_nparts=atoi(s2);
    else if(!strcmp(s1,"use_pm="))
      use_pm=atoi(s2);
    else if(!strcmp(s1,"use_randoms="))
      use_randoms=atoi(s2);
    else if(!strcmp(s1,"reuse_randoms="))
      reuse_randoms=atoi(s2);
    else if(!strcmp(s1,"do_CCF="))
      do_CCF=atoi(s2);
    else if(!strcmp(s1,"corr_type="))
      corr_type=atoi(s2);
    else if(!strcmp(s1,"n_grid_side="))
      n_grid=atoi(s2);
    else
      fprintf(stderr,"CUTE: Unknown parameter %s\n",s1);
  }
  fclose(fi);

  check_params();

  printf("\n");
}

static void gad_check_block(int b1,int b2)
{
  //////
  // Checks that a block from a snapshot is
  // consistent from its begin/end values
  if(b1!=b2) {
    fprintf(stderr,"CUTE: Corrupted block!\n");
    exit(1);
  }
}

static int gad_seek_block(FILE *snap,char name[])
{
  //////
  // Seeks block from title
  gad_title tit;
  int block1,block2;

  rewind(snap);

  while(1>0) {
    if(!(fread(&block1,sizeof(int),1,snap))||
        feof(snap)||ferror(snap)) {
      fprintf(stderr,"CUTE: Block %s not found!!\n",name);
      exit(1);
    }
    my_fread(&tit,sizeof(gad_title),1,snap);
    my_fread(&block2,sizeof(int),1,snap);
    gad_check_block(block1,block2);
    if(strncmp(tit.label,name,3)!=0)
      fseek(snap,tit.size,SEEK_CUR);
    else
      break;
  }

  return 0;
}

static Catalog *read_ascii(char *fname,lint *np)
{ 
  //////
  // Reads catalog from ascii file with
  // default format.
  Catalog *cat = (Catalog *)malloc(sizeof(Catalog));
  lint ii,n_lin;
  FILE *fd;

  //Open file and count lines
  fd=fopen(fname,"r");
  if(fd==NULL) error_open_file(fname);
  if(n_objects==-1) 
    n_lin=linecount(fd);
  else
    n_lin=n_objects;
  rewind(fd);

#ifdef _VERBOSE
  printf("  %ld objects will be read \n",(long)n_lin);
#endif

  //Allocate catalog memory
  *np=n_lin;
  cat->np=n_lin;
  cat->pos=(double *)malloc(3*cat->np*sizeof(double));
  if(cat->pos==NULL)
    error_mem_out();

  rewind(fd);
  //Read galaxies in mask
  for(ii=0;ii<n_lin;ii++) {
    double xx,yy,zz;
    char s0[1024];
    int sr;
    if(fgets(s0,sizeof(s0),fd)==NULL)
      error_read_line(fname,ii+1);
    sr=sscanf(s0,"%lf %lf %lf",&xx,&yy,&zz);
    if(sr!=3)
      error_read_line(fname,ii+1);
    cat->pos[3*ii]=xx;
    cat->pos[3*ii+1]=yy;
    cat->pos[3*ii+2]=zz;
  }
  fclose(fd);

  return cat;
}

static int check_num_files(char *prefix)
{
  FILE *fil;

  fil=fopen(prefix,"r");
  if(fil!=NULL) {
    fclose(fil);
    return 1;
  }
  else {
    int nfils=0;
    while(nfils>=0) {
      char fname[256];
      sprintf(fname,"%s.%d",prefix,nfils);
      fil=fopen(fname,"r");
      if(fil!=NULL) {
        fclose(fil);
        nfils++;
      }
      else {
        if(nfils==0) {
          fprintf(stderr,"CUTE: can't find file %s or %s.x\n",prefix,prefix);
          return -1;
        }
        else if(nfils==1) {
          fprintf(stderr,"CUTE: only file %s found. Weird.\n",fname);
          return -1;
        }
        else {
          return nfils;
        }
      }
    }
  }

  fprintf(stderr,"CUTE: this shouldn't have happened \n");
  return -1;
}

static Catalog *read_snapshot_single(char *fname,lint *np,int input)
{
  //////
  // Creates catalog from a single snapshot file
  Catalog *cat = (Catalog *)malloc(sizeof(Catalog));
  lint ii;
  gad_header head;
  int block1,block2;

  FILE *snap=fopen(fname,"r");
  if(snap==NULL) error_open_file(fname);

  //Read header
  if(input==2)
    gad_seek_block(snap,"HEAD");
  my_fread(&block1,sizeof(int),1,snap);
  my_fread(&head,sizeof(gad_header),1,snap);
  my_fread(&block2,sizeof(int),1,snap);
  gad_check_block(block1,block2);

  if(head.num_files!=1) {
    fprintf(stderr,"CUTE: Multi-file input not expected \n");
    exit(1);
  }

#ifdef _VERBOSE
  printf("  The cosmological model is:\n");
  printf("   - Omega_M = %.3lf\n",head.Omega0);
  printf("   - Omega_L = %.3lf\n",head.OmegaLambda);
  printf("   - h = %.3lf\n",head.HubbleParam);
  printf("  This file contains: \n");
  for(ii=0;ii<6;ii++) {
    printf("   - %d particles of type %d with mass",
        head.npart[ii],(int)ii);
    printf(" %.3lE (%d in total)\n",
        head.mass[ii],head.npartTotal[ii]);
  }
  printf("  The box size is %.3lf\n",head.BoxSize);
  printf("  Redshift z = %.3lf \n",head.redshift);
#endif //_VERBOSE

  l_box=head.BoxSize;
  l_box_half=l_box*0.5;
  cat->np=0;
  for(ii=0;ii<6;ii++) {
    cat->np+=head.npart[ii];
    if(head.npart[ii]!=head.npartTotal[ii]) {
      fprintf(stderr,"CUTE: error reading snapshot \n");
      exit(1);
    }
  }
  *np=cat->np;

  cat->pos=(double *)malloc(3*cat->np*sizeof(double));
  if(cat->pos==NULL)
    error_mem_out();

  if(input==2)
    gad_seek_block(snap,"POS");
  my_fread(&block1,sizeof(int),1,snap);
  for(ii=0;ii<cat->np;ii++) {
    float pos[3];
    my_fread(pos,sizeof(float),3,snap);
    cat->pos[3*ii]=(double)(pos[0]);
    cat->pos[3*ii+1]=(double)(pos[1]);
    cat->pos[3*ii+2]=(double)(pos[2]);
  }
  my_fread(&block2,sizeof(int),1,snap);
  gad_check_block(block1,block2);
  fclose(snap);

  return cat;
}

static Catalog *read_gadget(char *prefix,lint *np,int input)
{
  Catalog *cat = (Catalog *)malloc(sizeof(Catalog));
  int nfils=check_num_files(prefix);
  if(nfils<=0) exit(1);

#ifdef _VERBOSE
  printf("  Reading from GADGET snapshot format \n");
#endif //_VERBOSE

  if(nfils==1) {
    printf("  Reading single snapshot file\n");
    cat=read_snapshot_single(prefix,np,input);
    return cat;
  }
  else {
    lint ii;
    char fname[256];
    gad_header head;
    int block1,block2;
    FILE *snap;

    printf("  Reading %d snapshot files \n",nfils);
    //    fprintf(stderr,"CUTE: multi-file input not supported \n");
    //    exit(1);

    sprintf(fname,"%s.0",prefix);
    snap=fopen(fname,"r");
    if(snap==NULL) error_open_file(fname);

    //Read header
    if(input==2)
      gad_seek_block(snap,"HEAD");
    my_fread(&block1,sizeof(int),1,snap);
    my_fread(&head,sizeof(gad_header),1,snap);
    my_fread(&block2,sizeof(int),1,snap);
    gad_check_block(block1,block2);

    if(head.num_files!=nfils) {
      fprintf(stderr,
          "CUTE: Header and existing files do not match %d != %d.\n",
          nfils,head.num_files);
      fprintf(stderr,"      There may be some files missing\n");
      exit(1);
    }

#ifdef _VERBOSE
    printf("  The cosmological model is:\n");
    printf("   - Omega_M = %.3lf\n",head.Omega0);
    printf("   - Omega_L = %.3lf\n",head.OmegaLambda);
    printf("   - h = %.3lf\n",head.HubbleParam);
    printf("  This file contains: \n");
    for(ii=0;ii<6;ii++) {
      printf("   - %d particles of type %d with mass %.3lE\n",
          head.npartTotal[ii],(int)ii,head.mass[ii]);
    }
    printf("  The box size is %.3lf\n",head.BoxSize);
    printf("  Redshift z = %.3lf \n",head.redshift);
#endif //_VERBOSE

    l_box=head.BoxSize;
    l_box_half=l_box*0.5;
    cat->np=0;
    for(ii=0;ii<6;ii++)
      cat->np+=head.npartTotal[ii];
    *np=cat->np;

    cat->pos=(double *)malloc(3*cat->np*sizeof(double));
    if(cat->pos==NULL)
      error_mem_out();
    fclose(snap);

    lint np_read=0;
    for(ii=0;ii<nfils;ii++) {
      lint np_new;
      lint jj;

      sprintf(fname,"%s.%d",prefix,(int)ii);
      snap=fopen(fname,"r");
      if(snap==NULL) error_open_file(fname);
#ifdef _VERBOSE
      printf("  Reading file  %s \n",fname);
#endif //_VERBOSE

      //Read header
      if(input==2)
        gad_seek_block(snap,"HEAD");
      my_fread(&block1,sizeof(int),1,snap);
      my_fread(&head,sizeof(gad_header),1,snap);
      my_fread(&block2,sizeof(int),1,snap);
      gad_check_block(block1,block2);

      np_new=0;
      for(jj=0;jj<6;jj++)
        np_new+=head.npart[jj];
      printf("  %ld parts in file %ld \n",(long)np_new,(long)ii);

      if(np_read+np_new>cat->np) {
        fprintf(stderr,
            "CUTE: files seem to contain too many particles\n");
        fprintf(stderr,"      file %s, %ld > %ld \n",
            fname,(long)(np_read+np_new),(long)(cat->np));
        exit(1);
      }

      if(input==2)
        gad_seek_block(snap,"POS");
      my_fread(&block1,sizeof(int),1,snap);
      for(jj=np_read;jj<np_read+np_new;jj++) {
        float pos[3];
        my_fread(pos,sizeof(float),3,snap);
        cat->pos[3*jj]=(double)(pos[0]);
        cat->pos[3*jj+1]=(double)(pos[1]);
        cat->pos[3*jj+2]=(double)(pos[2]);
      }
      my_fread(&block2,sizeof(int),1,snap);
      gad_check_block(block1,block2);
      fclose(snap);

      np_read+=np_new;
    }
    if(np_read!=cat->np) {
      fprintf(stderr,
          "CUTE: #particles read disagrees with header: %ld != %ld\n",
          (long)np_read,(long)(cat->np));
      exit(1);
    }

    return cat;
  }
}

Catalog *read_catalog(char *fname,lint *np)
{
  //////
  // Creates catalog from file fname
  lint ii;
  double x_mean=0,y_mean=0,z_mean=0;
  Catalog *cat = (Catalog *)malloc(sizeof(Catalog));

  printf("*** Reading catalog ");
#ifdef _VERBOSE
  printf("from file %s",fname);
#endif
  printf("\n");

  if(input_format)
    cat=read_gadget(fname,np,input_format);
  else
    cat=read_ascii(fname,np);

  //Correct particles out of bounds and calculate CoM
  for(ii=0;ii<cat->np;ii++) {
    double xx,yy,zz;
    xx=cat->pos[3*ii];
    yy=cat->pos[3*ii+1];
    zz=cat->pos[3*ii+2];
    if((xx<0)||(xx>=l_box)) xx=wrap_double(xx);
    if((yy<0)||(yy>=l_box)) yy=wrap_double(yy);
    if((zz<0)||(zz>=l_box)) zz=wrap_double(yy);
    cat->pos[3*ii]=xx;
    cat->pos[3*ii+1]=yy;
    cat->pos[3*ii+2]=zz;
    x_mean+=xx/cat->np;
    y_mean+=yy/cat->np;
    z_mean+=zz/cat->np;
  }

#ifdef _VERBOSE
  printf("  The center of mass is (%.3lf,%.3lf,%.3lf) \n",
      x_mean,y_mean,z_mean);
#endif //_VERBOSE

  printf("\n");
  return cat;
}


#ifdef _CUTE_AS_PYTHON_MODULE

int verify_parameters(){
  param_errors = 0;
  check_params();
  printf("Checking CUTE parameters. Total error count: %i\n",param_errors);
  return param_errors;
}
Catalog *read_Catalog(char *fname){
  lint n;
  return read_catalog(fname, &n);
}
void set_data_filename(char *s){
  sprintf(fnameData,"%s",s);
}
void set_data_filename2(char *s){
  sprintf(fnameData2,"%s",s);
}
void set_random_filename(char *s){
  sprintf(fnameRand,"%s",s);
}
void set_input_format(int i){
  input_format=i;
}
void set_output_filename(char *s){
  sprintf(fnameOut,"%s",s);
}
void set_corr_type(int i){
  if(i < 1 || i > 3) {
    fprintf(stderr,"CUTE: wrong corr type %i.",i);
    fprintf(stderr," Possible types are 1,2 or 3");
  }
  corr_type = i;
}
void set_use_pm(int i){
  use_pm=i;
}
void set_reuse_randoms(int i){
  reuse_randoms = i;
}
void set_num_lines(char *s){
  if(!strcmp(s,"all"))
    n_objects=-1;
  else
    n_objects=atoi(s);
}
void set_box_size(double x){
  l_box=x;
  l_box_half=0.5*l_box;
}
void set_use_tree(int i){
  use_tree = i;
}
void set_max_tree_order(int i){
  max_tree_order = i;
}
void set_max_tree_nparts(int i){
  max_tree_nparts = i;
}
void set_use_randoms(int i){
  use_randoms = i;
}
void set_do_CCF(int i){
  do_CCF = i;
}
void set_n_grid_side(int i){
  n_grid = i;
}
#endif
