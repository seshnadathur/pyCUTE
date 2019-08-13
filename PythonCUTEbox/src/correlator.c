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
//                      Correlators with OpenMP                      //
/*********************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "define.h"
#include "common.h"

//#define _DO_BATCHES
#ifdef _DO_BATCHES
#define COUNT_LIM 1000000000
#endif //_DO_BATCHES

static const double I_DR=I_R_MAX*NB_R;
static const double R2_MAX=1./(I_R_MAX*I_R_MAX);

void corr_mono_box_bf(lint np,double *pos,
    unsigned long long hh[])
{
  //////
  // Correlator for monopole in the periodic-box case
  // by brute-force
  int i;
  for(i=0;i<NB_R;i++)
    hh[i]=0; //Clear shared histogram

#pragma omp parallel default(none)		\
  shared(pos,np,hh,l_box,l_box_half)
  {
    lint ii;
    unsigned long long hthread[NB_R]; //Histogram filled by each thread

    for(ii=0;ii<NB_R;ii++)
      hthread[ii]=0; //Clear private histogram

#pragma omp for nowait schedule(dynamic)
    for(ii=0;ii<np;ii++) {
      lint jj;
      double *pos1=&(pos[3*ii]);
      for(jj=0;jj<np;jj++) {
        double xr[3];
        double r2;
        int ir;
        xr[0]=ABS(pos1[0]-pos[3*jj]);
        xr[1]=ABS(pos1[1]-pos[3*jj+1]);
        xr[2]=ABS(pos1[2]-pos[3*jj+2]);
        if(xr[0]>l_box_half) xr[0]=l_box-xr[0];
        if(xr[1]>l_box_half) xr[1]=l_box-xr[1];
        if(xr[2]>l_box_half) xr[2]=l_box-xr[2];
        r2=xr[0]*xr[0]+xr[1]*xr[1]+xr[2]*xr[2]; //Relative distance squared
        if(r2>R2_MAX) continue;
#ifdef _LOGBIN
        if(r2>0) {
          ir=(int)(N_LOGINT*(0.5*log10(r2)-LOG_R_MAX)+NB_R);
          if((ir<NB_R)&&(ir>=0))
            (hthread[ir])++;
        }
#else //_LOGBIN
        ir=(int)(sqrt(r2)*I_DR);
        if(ir<NB_R) //Check bound
          (hthread[ir])++;
#endif //_LOGBIN
      }
    } // end pragma omp for
#pragma omp critical
    {
      for(ii=0;ii<NB_R;ii++) //Check bound
        hh[ii]+=hthread[ii]; //Add private histograms to shared one
    } // end pragma omp critical
  } // end pragma omp parallel
}

void corr_mono_box_pm(double *grid,double corr[],double ercorr[],
    unsigned long long hh[])
{
  //////
  // Correlator for monopole in the periodic-box case
  // by using a particle mesh
  int *ibin_box;
  double agrid=l_box/n_grid;
  double agrid2=agrid*agrid;
  double r_max=MIN(1/I_R_MAX,l_box_half);
  lint index_max=(int)(r_max/agrid)+1;
  lint i;

  for(i=0;i<NB_R;i++) {
    hh[i]=0;
    corr[i]=0;
    ercorr[i]=0;
  }

  printf("Using a distance cube of order %ld for r_max = %.3lf \n",
      (long)index_max,r_max);
  ibin_box=(int *)malloc(index_max*index_max*index_max*sizeof(int));
  if (ibin_box==NULL) error_mem_out();

  for(i=0;i<index_max;i++) {	// calculate distances for cells within the distance cube
    lint j;
    for(j=0;j<index_max;j++) {
      lint k;
      for(k=0;k<index_max;k++) {
        lint ir2=i*i+j*j+k*k;
        double r2=agrid2*ir2;
        int ibin;
        ibin=(int)(sqrt(r2)*I_DR);
        ibin_box[k+index_max*(j+index_max*i)]=ibin;
      }
    }
  }

  for(i=-index_max+1;i<index_max;i++) {	// count number of cells in distance cube within distance bins
    lint j;
    for(j=-index_max+1;j<index_max;j++) {
      lint k;
      for(k=-index_max+1;k<index_max;k++) {
        int ibin=ibin_box[labs(k)+index_max*
          (labs(j)+index_max*labs(i))];
        if(ibin<NB_R) hh[ibin]++;
      }
    }
  }

#pragma omp parallel default(none)			\
  shared(n_grid,grid,index_max,ibin_box,corr)
  {
    lint ii;
    double corr_thr[NB_R];
    lint n_grid2=n_grid*n_grid,index_max2=index_max*index_max;
    int ngm1=n_grid-1;
#ifdef _DO_BATCHES
    double corr_batch[NB_R];
    unsigned int hh_batch[NB_R];
#endif //_DO_BATCHES

    for(ii=0;ii<NB_R;ii++) {
      corr_thr[ii]=0;
#ifdef _DO_BATCHES
      corr_batch[ii]=0;
      hh_batch[ii]=0;
#endif //_DO_BATCHES
    }

#pragma omp for nowait schedule(dynamic)
    for(ii=0;ii<n_grid2*n_grid;ii++) {
      lint i1=ii/(n_grid2);
      lint k1=ii%n_grid;
      lint j1=(ii-k1-i1*n_grid2)/n_grid;
      double d1=grid[ii];
      lint ir;
      for(ir=-index_max+1;ir<index_max;ir++) {
        lint jr;
        lint i2=(i1+ir)&ngm1;
        lint irr=labs(ir)*index_max2;
        i2*=n_grid2;
        for(jr=-index_max+1;jr<index_max;jr++) {
          lint kr;
          lint j2=(j1+jr)&ngm1;
          lint jrr=labs(jr)*index_max;
          j2*=n_grid;
          for(kr=-index_max+1;kr<index_max;kr++) {
            int ibin=ibin_box[labs(kr)+jrr+irr];
            if(ibin>=NB_R) continue;
            else {
              double d2;
              lint k2=(k1+kr)&ngm1;
              d2=d1*grid[k2+j2+i2];
#ifdef _DO_BATCHES
              corr_batch[ibin]+=d2;
              hh_batch[ibin]++;
              if(hh_batch[ibin]>COUNT_LIM) {
                corr_thr[ibin]+=corr_batch[ibin];
                hh_batch[ibin]=0;
                corr_batch[ibin]=0;
              }
#else //_DO_BATCHES
              corr_thr[ibin]+=d2;
#endif //_DO_BATCHES
            }
          }
        }
      }
    } //end pragma omp for

#ifdef _DO_BATCHES
    for(ii=0;ii<NB_R;ii++)
      corr_thr[ii]+=corr_batch[ii];
#endif //_DO_BATCHES

#pragma omp critical
    {
      for(ii=0;ii<NB_R;ii++)
        corr[ii]+=corr_thr[ii];
    } //end pragma omp critical
  } //end pragma omp parallel

  for(i=0;i<NB_R;i++) {
    hh[i]*=(n_grid*((lint)(n_grid*n_grid)));
    if(hh[i]>0) corr[i]/=hh[i];
    else corr[i]=0;
  }

  free(ibin_box);

  return;
}

static int limit_dist2_point2box(double *x,float *x_l,float *x_h,
    float *d2_l,float *d2_h)
{
  //////
  // Returns, in d2_l and d2_h, the minimum and maximum distance
  // squared from point x[3] to box defined by x_l[3] and x_m[3]
  int ii;

  for(ii=0;ii<3;ii++) {
    float d_l,d_h;
    if(x[ii]<x_l[ii]) {
      d_l=x_l[ii]-x[ii];
      d_h=x_h[ii]-x[ii];
    }
    else if(x[ii]>x_h[ii]) {
      d_l=x[ii]-x_h[ii];
      d_h=x[ii]-x_l[ii];
    }
    else {
      d_l=0;
      d_h=MAX(x_h[ii]-x[ii],x[ii]-x_l[ii]);
    }

    if(d_l>l_box_half) { //Branch entirely out of range, wrap around
      d_l=l_box-d_h;
      d_h=l_box-d_l;
    }
    else if(d_h>l_box_half) return 1; //Branch partly out or range, report

    *d2_l+=d_l*d_l;
    *d2_h+=d_h*d_h;
  }

  return 0;
}

static void bin_branch(branch *br,double *x,unsigned long long hh[])
{
  //////
  // Bins branch br into histogram hh according to distance to x[3].
  // Considers different possibilities:
  //    -Branch is empty or not
  //    -Branch/leaf is partly beyond l_box_half
  //    -Branch is completely out of range
  //    -Branch fits inside one bin
  //    -Branch/leaf spans several bins
  if(br->np) { //Bin only if branch is non-empty
    lint ii;
    float d2_l=0,d2_h=0;
    //First, calculate branch bounds and whether it is partly out of range
    if(limit_dist2_point2box(x,br->x_lo,br->x_hi,&d2_l,&d2_h)) {
      if(br->leaf) { //Leaf partly out of range. Iterate
        for(ii=0;ii<br->np;ii++) {
          double *pos=((double *)(br->sons)+3*ii);
          double r2;
          int ir;
          double xr[3];

          xr[0]=ABS(pos[0]-x[0]);
          xr[1]=ABS(pos[1]-x[1]);
          xr[2]=ABS(pos[2]-x[2]);
          if(xr[0]>l_box_half) xr[0]=l_box-xr[0];
          if(xr[1]>l_box_half) xr[1]=l_box-xr[1];
          if(xr[2]>l_box_half) xr[2]=l_box-xr[2];

          r2=xr[0]*xr[0]+xr[1]*xr[1]+xr[2]*xr[2];
          if(r2>R2_MAX) continue;
#ifdef _LOGBIN
          if(r2>0) {
            ir=(int)(N_LOGINT*(0.5*log10(r2)-LOG_R_MAX)+NB_R);
            if((ir>=0)&&(ir<NB_R))
              hh[ir]++;
          }
#else //_LOGBIN
          ir=(int)(sqrt(r2)*I_DR);
          if(ir<NB_R)
            hh[ir]++;
#endif //_LOGBIN
        }
      }
      else { //Branch partly out of range, open
        for(ii=0;ii<8;ii++)
          bin_branch(((branch **)(br->sons))[ii],x,hh);
      }
    }
    else if(d2_l>=R2_MAX) return;
    else { //If branch is completely in range or wrapped around
      int ir_l,ir_h;

      //Calculate bins for bounds
#ifdef _LOGBIN
      if(d2_l>0) {
        ir_l=(int)(N_LOGINT*(0.5*log10(d2_l)-LOG_R_MAX)+NB_R);
        ir_h=(int)(N_LOGINT*(0.5*log10(d2_h)-LOG_R_MAX)+NB_R);
      }
      else {
        ir_l=-2;
        ir_h=5;
      }
#else //_LOGBIN
      ir_l=(int)(sqrtf(d2_l)*I_DR);
      ir_h=(int)(sqrtf(d2_h)*I_DR);
#endif //_LOGBIN

      if((ir_l==ir_h)&&(ir_l>=0)) //If branch completely inside a bin
        hh[ir_l]+=br->np;
      else { //If Branch spans several bins
        if(br->leaf) { //If leaf, iterate
          for(ii=0;ii<br->np;ii++) {
            double *pos=((double *)(br->sons)+3*ii);
            double r2;
            int ir;
            double xr[3];
            xr[0]=ABS(pos[0]-x[0]);
            xr[1]=ABS(pos[1]-x[1]);
            xr[2]=ABS(pos[2]-x[2]);
            if(xr[0]>l_box_half) xr[0]=l_box-xr[0];
            if(xr[1]>l_box_half) xr[1]=l_box-xr[1];
            if(xr[2]>l_box_half) xr[2]=l_box-xr[2];

            r2=xr[0]*xr[0]+xr[1]*xr[1]+xr[2]*xr[2];
            if(r2>R2_MAX) continue;
#ifdef _LOGBIN
            if(r2>0) {
              ir=(int)(N_LOGINT*(0.5*log10(r2)-LOG_R_MAX)+NB_R);
              if((ir>=0)&&(ir<NB_R))
                hh[ir]++;
            }
#else //_LOGBIN
            ir=(int)(sqrt(r2)*I_DR);
            if(ir<NB_R)
              hh[ir]++;
#endif //_LOGBIN
          }
        }
        else { //If branch, open
          for(ii=0;ii<8;ii++)
            bin_branch(((branch **)(br->sons))[ii],x,hh);
        }
      }
    }
  }
}

void corr_mono_box_tree(lint np,double *pos,branch *tree,
    unsigned long long hh[])
{
  //////
  // Correlator for monopole in the periodic-box case
  // using a tree
  int i;

  for(i=0;i<NB_R;i++)
    hh[i]=0;

#pragma omp parallel default(none)		\
  shared(np,pos,tree,hh)
  {
    lint ii;
    unsigned long long hthread[NB_R];

    for(ii=0;ii<NB_R;ii++)
      hthread[ii]=0;

#pragma omp for nowait schedule(dynamic)
    for(ii=0;ii<np;ii++) {
      bin_branch(tree,&(pos[3*ii]),hthread);
    } //end pragma omp for

#pragma omp critical
    {
      for(ii=0;ii<NB_R;ii++)
        hh[ii]+=hthread[ii];
    } //end pragma omp critical
  } //end pragma omp parallel
}

void corr_mono_box_neighbors(int nside,NeighborBox *boxes,
    lint np,double *pos,
    unsigned long long hh[])
{
  //////
  // Correlator for monopole in the periodic-box case
  // using neighbor boxes - this original version counts each pair twice
  // and counts each particle as a pair with itself
  double agrid=l_box/nside;
  double r_max=1/I_R_MAX;
  int index_max=(int)(r_max/agrid)+1;
  int i;

  printf("  Boxes will be correlated up to %d box sizes \n",index_max);

  for(i=0;i<NB_R;i++)
    hh[i]=0; //Clear shared histogram

#pragma omp parallel default(none)			\
  shared(index_max,nside,boxes,hh,l_box,np,pos,agrid)
  {
    lint ii;
    double a2grid=agrid*agrid;
    unsigned long long hthread[NB_R]; //Histogram filled by each thread

    for(ii=0;ii<NB_R;ii++)
      hthread[ii]=0; //Clear private histogram

#pragma omp for nowait schedule(dynamic)
    for(ii=0;ii<np;ii++) {	// loop over all particles
      int ix0,iy0,iz0;
      double x0,y0,z0;
      int idz;
      x0=pos[3*ii];
      y0=pos[3*ii+1];
      z0=pos[3*ii+2];
      ix0=(int)(x0/l_box*nside); // find sub-box it lies in
      iy0=(int)(y0/l_box*nside);
      iz0=(int)(z0/l_box*nside);

      for(idz=-index_max;idz<=index_max;idz++) {
        int idy,idz_dist2;
        int iwrapz=0;
        int iz1=iz0+idz;
        if(iz1<0) {
          iz1+=nside;
          iwrapz=1;
        }
        else if(iz1>=nside) {
          iz1-=nside;
          iwrapz=1;
        }
        idz_dist2=MAX(0,abs(idz)-1);
        idz_dist2=idz_dist2*idz_dist2;
        for(idy=-index_max;idy<=index_max;idy++) {
          int idx,idy_dist2;
          int iwrapy=0;
          int iy1=iy0+idy;
          if(iy1<0) {
            iy1+=nside;
            iwrapy=1;
          }
          else if(iy1>=nside) {
            iy1-=nside;
            iwrapy=1;
          }
          idy_dist2=MAX(0,abs(idy)-1);
          idy_dist2=idy_dist2*idy_dist2;
          for(idx=-index_max;idx<=index_max;idx++) {
            int ibox,idx_dist;
            int iwrapx=0;
            int ix1=ix0+idx;
            double d2max;
            int jj;
            if(ix1<0) {
              ix1+=nside;
              iwrapx=1;
            }
            else if(ix1>=nside) {
              ix1-=nside;
              iwrapx=1;
            }
            ibox=ix1+nside*(iy1+nside*iz1); // index of second sub-box (may be same as original)
            idx_dist=MAX(0,abs(idx)-1);
            d2max=a2grid*(idx_dist*idx_dist+idy_dist2+idz_dist2); // box-to-box distance
            if(d2max>R2_MAX) continue;
            for(jj=0;jj<boxes[ibox].np;jj++) {	// loop over all particles in second box (may include original) 
              double xr[3];
              double r2;
              int ir;
              xr[0]=fabs(x0-(boxes[ibox].pos)[3*jj]);
              xr[1]=fabs(y0-(boxes[ibox].pos)[3*jj+1]);
              xr[2]=fabs(z0-(boxes[ibox].pos)[3*jj+2]);
              if(iwrapx) xr[0]=l_box-xr[0];	//check PBC due to wrapping
              if(iwrapy) xr[1]=l_box-xr[1];
              if(iwrapz) xr[2]=l_box-xr[2];
              r2=xr[0]*xr[0]+xr[1]*xr[1]+xr[2]*xr[2];
              if(r2>R2_MAX) continue;
#ifdef _LOGBIN
              if(r2>0) {
                ir=(int)(N_LOGINT*(0.5*log10(r2)-LOG_R_MAX)+NB_R);
                if((ir<NB_R)&&(ir>=0))
                  (hthread[ir])++;
              }
#else //_LOGBIN
              ir=(int)(sqrt(r2)*I_DR);
              if(ir<NB_R) //Check bound
                (hthread[ir])++;
#endif //_LOGBIN
            }
          }
        }
      }
    } // end pragma omp for

#pragma omp critical
    {
      for(ii=0;ii<NB_R;ii++) //Check bound
        hh[ii]+=hthread[ii]; //Add private histograms to shared one
    } // end pragma omp critical
  } // end pragma omp parallel
}

void corr_mono_boxes(int nside,NeighborBox *boxes,
    unsigned long long hh[])
{
  //////
  // Correlator for monopole in the periodic-box case
  // using neighbor boxes - this version counts each pair only once
  double agrid=l_box/nside;
  double r_max=1/I_R_MAX;
  int index_max=(int)(r_max/agrid)+1;
  int i;

  printf("  Boxes will be correlated up to %d box sizes \n",index_max);

  for(i=0;i<NB_R;i++)
    hh[i]=0; //Clear shared histogram

#pragma omp parallel default(none)			\
  shared(index_max,nside,boxes,hh,l_box,agrid)
  {
    lint ii,ibox;
    unsigned long long hthread[NB_R]; //Histogram filled by each thread

    for(ii=0;ii<NB_R;ii++)
      hthread[ii]=0; //Clear private histogram

#pragma omp for nowait schedule(dynamic)
    for(ibox=0;ibox<nside*nside*nside;ibox++) { //loop over sub-boxes
      int ix0,iy0,iz0;
      int idz,np_box,np_2box;
      double x0,y0,z0,xr[3],r2;
      int ir;
      lint jj,this_box;

      ix0=ibox%nside;		// coordinates of current sub-box
      iy0=(ibox%(nside*nside))/nside;
      iz0=ibox/(nside*nside);

      if(boxes[ibox].np>0) {	//the sub-box is not empty
        np_box = boxes[ibox].np;
        for(ii=0;ii<np_box;ii++)  { 	//loop over particles in the sub-box
          x0=(boxes[ibox].pos)[3*ii];				
          y0=(boxes[ibox].pos)[3*ii+1];
          z0=(boxes[ibox].pos)[3*ii+2];
          for(jj=ii+1;jj<np_box;jj++) { //loop over pairs; jj=ii+1 to start ensures each pair counted only once
            xr[0]=fabs(x0-(boxes[ibox].pos)[3*jj]);
            xr[1]=fabs(y0-(boxes[ibox].pos)[3*jj+1]);
            xr[2]=fabs(z0-(boxes[ibox].pos)[3*jj+2]);
            r2=xr[0]*xr[0]+xr[1]*xr[1]+xr[2]*xr[2];
            if(r2>R2_MAX) continue;
#ifdef _LOGBIN
            if(r2>0) {
              ir=(int)(N_LOGINT*(0.5*log10(r2)-LOG_R_MAX)+NB_R);
              if((ir<NB_R)&&(ir>=0))
                (hthread[ir])++;
            }
#else //_LOGBIN
            ir=(int)(sqrt(r2)*I_DR);
            if(ir<NB_R) //Check bound
              (hthread[ir])++;
#endif //_LOGBIN			
          }	//end for loop jj
        }	//end loop over ii
        //now look for nearby sub-boxes
        for(idz=-index_max;idz<=index_max;idz++) { 
          int idy;
          int iwrapz=0;
          int iz1=iz0+idz;
          if(iz1<0) {
            iz1+=nside;
            iwrapz=1;
          }
          else if(iz1>=nside) {
            iz1-=nside;
            iwrapz=1;
          }
          for(idy=-index_max;idy<=index_max;idy++) {
            int idx;
            int iwrapy=0;
            int iy1=iy0+idy;
            if(iy1<0) {
              iy1+=nside;
              iwrapy=1;
            }
            else if(iy1>=nside) {
              iy1-=nside;
              iwrapy=1;
            }
            for(idx=-index_max;idx<=index_max;idx++) {
              int iwrapx=0;
              int ix1=ix0+idx;
              if(ix1<0) {
                ix1+=nside;
                iwrapx=1;
              }
              else if(ix1>=nside) {
                ix1-=nside;
                iwrapx=1;
              }
              this_box=ix1+nside*(iy1+nside*iz1);		//index of nearby sub-box
              if((this_box>ibox)&&(boxes[this_box].np>0)) {	//only count pairs of sub-boxes once
                np_2box = boxes[this_box].np;
                for(ii=0;ii<np_box;ii++)  {
                  x0=(boxes[ibox].pos)[3*ii];				
                  y0=(boxes[ibox].pos)[3*ii+1];
                  z0=(boxes[ibox].pos)[3*ii+2];
                  for(jj=0;jj<np_2box;jj++) { 
                    xr[0]=fabs(x0-(boxes[this_box].pos)[3*jj]);
                    xr[1]=fabs(y0-(boxes[this_box].pos)[3*jj+1]);
                    xr[2]=fabs(z0-(boxes[this_box].pos)[3*jj+2]);
                    if(iwrapx) xr[0]=l_box-xr[0];	//check PBC due to wrapping
                    if(iwrapy) xr[1]=l_box-xr[1];
                    if(iwrapz) xr[2]=l_box-xr[2];
                    r2=xr[0]*xr[0]+xr[1]*xr[1]+xr[2]*xr[2];
                    if(r2>R2_MAX) continue;
#ifdef _LOGBIN
                    if(r2>0) {
                      ir=(int)(N_LOGINT*(0.5*log10(r2)-LOG_R_MAX)+NB_R);
                      if((ir<NB_R)&&(ir>=0))
                        (hthread[ir])++;
                    }
#else //_LOGBIN
                    ir=(int)(sqrt(r2)*I_DR);
                    if(ir<NB_R) //Check bound
                      (hthread[ir])++;
#endif //_LOGBIN
                  }	//end for loop jj
                }	//end for loop ii
              }	//end if this_box
            }	//end for idx
          }	//end for idy
        }	//end for idz
      }	//end if boxes[ibox].np>0
    } // end omp for over ibox

#pragma omp critical
    {
      for(ii=0;ii<NB_R;ii++) //Check bound
        hh[ii]+=hthread[ii]; //Add private histograms to shared one
    } // end omp critical

  } // end omp parallel

}

void crosscorr_mono_box_neighbors(int nside,NeighborBox *boxes1,
    NeighborBox *boxes2,unsigned long long hh[])
{
  //////
  // Cross-correlator for monopole in the periodic-box case
  // using neighbor boxes  
  double agrid=l_box/nside;
  double r_max=1/I_R_MAX;
  int index_max=(int)(r_max/agrid)+1;
  int i;

  printf("  Boxes will be correlated up to %d box sizes \n",index_max);

  for(i=0;i<NB_R;i++)
    hh[i]=0; //Clear shared histogram

#pragma omp parallel default(none)			\
  shared(index_max,nside,boxes1,boxes2,hh,l_box,agrid)
  {
    lint ii,ibox;
    unsigned long long hthread[NB_R]; //Histogram filled by each thread

    for(ii=0;ii<NB_R;ii++)
      hthread[ii]=0; //Clear private histogram

#pragma omp for nowait schedule(dynamic)
    for(ibox=0;ibox<nside*nside*nside;ibox++) {	//loop over sub-boxes
      int ix0,iy0,iz0;
      int idz,np1_box,np2_box;
      double x0,y0,z0,xr[3],r2;
      int ir;
      lint jj,this_box;

      ix0=ibox%nside;		// coordinates of current sub-box
      iy0=(ibox%(nside*nside))/nside;
      iz0=ibox/(nside*nside);

      for(idz=-index_max;idz<=index_max;idz++) { //loop over boxes adjacent in z-direction
        int idy;
        int iwrapz=0;
        int iz1=iz0+idz;
        if(iz1<0) {	//check for PBC wraparound
          iz1+=nside;
          iwrapz=1;
        }
        else if(iz1>=nside) {
          iz1-=nside;
          iwrapz=1;
        }
        for(idy=-index_max;idy<=index_max;idy++) { //loop over boxes adjacent in y-direction
          int idx;
          int iwrapy=0;
          int iy1=iy0+idy;
          if(iy1<0) {	//check for PBC wraparound
            iy1+=nside;
            iwrapy=1;
          }
          else if(iy1>=nside) {
            iy1-=nside;
            iwrapy=1;
          }
          for(idx=-index_max;idx<=index_max;idx++) { //loop over boxes adjacent in x-direction
            int iwrapx=0;
            int ix1=ix0+idx;
            if(ix1<0) {	//check for PBC wraparound
              ix1+=nside;
              iwrapx=1;
            }
            else if(ix1>=nside) {
              ix1-=nside;
              iwrapx=1;
            }
            this_box=ix1+nside*(iy1+nside*iz1);
            np1_box=boxes1[ibox].np;
            np2_box=boxes2[this_box].np;
            for(ii=0;ii<np1_box;ii++) { //otherwise, find distances to particle pairs in this box
              x0=(boxes1[ibox].pos)[3*ii];				
              y0=(boxes1[ibox].pos)[3*ii+1];
              z0=(boxes1[ibox].pos)[3*ii+2];
              for(jj=0;jj<np2_box;jj++) {	
                xr[0]=fabs(x0-(boxes2[this_box].pos)[3*jj]);	//calculate distance between particles
                xr[1]=fabs(y0-(boxes2[this_box].pos)[3*jj+1]);
                xr[2]=fabs(z0-(boxes2[this_box].pos)[3*jj+2]);
                if(iwrapx) xr[0]=l_box-xr[0];	//check PBC due to wrapping
                if(iwrapy) xr[1]=l_box-xr[1];
                if(iwrapz) xr[2]=l_box-xr[2];
                r2=xr[0]*xr[0]+xr[1]*xr[1]+xr[2]*xr[2];
                if(r2>R2_MAX) continue;
#ifdef _LOGBIN
                if(r2>0) {
                  ir=(int)(N_LOGINT*(0.5*log10(r2)-LOG_R_MAX)+NB_R);
                  if((ir<NB_R)&&(ir>=0))
                    (hthread[ir])++;
                }
#else //_LOGBIN
                ir=(int)(sqrt(r2)*I_DR);
                if(ir<NB_R) //Check bound
                  (hthread[ir])++;
#endif //_LOGBIN
              }	//end for loop jj
            }	//end for loop ii
          }	//end for loop idx
        }	//end for loop idy
      }	//end for loop idz
    } // end pragma omp for (ibox)

#pragma omp critical
    {
      for(ii=0;ii<NB_R;ii++) //Check bound
        hh[ii]+=hthread[ii]; //Add private histograms to shared one
    } // end pragma omp critical
  } // end pragma omp parallel
}

void auto_3d_ps_boxes(int nside,NeighborBox *boxes,
    unsigned long long hh[])
{
  //////
  // Correlator for xi(pi,sigma) in the periodic-box case
  // counts each pair only once, does not count self-pairs
  double agrid=l_box/nside;
  double r2_max=2./(I_R_MAX*I_R_MAX);
  //double rt2_max=1./(I_R_MAX*I_R_MAX);
  int index_max=(int)(sqrt(r2_max)/agrid)+1;
  int i;

  printf("  Boxes will be correlated up to %d box sizes \n",index_max);

  for(i=0;i<NB_R*NB_R;i++)
    hh[i]=0; //Clear shared histogram

#pragma omp parallel default(none)			\
  shared(index_max,nside,boxes,hh,l_box,agrid)
  {
    lint ii,ibox;
    unsigned long long hthread[NB_R*NB_R]; //Histogram filled by each thread

    for(ii=0;ii<NB_R*NB_R;ii++)
      hthread[ii]=0; //Clear private histogram

#pragma omp for nowait schedule(dynamic)
    for(ibox=0;ibox<nside*nside*nside;ibox++) { //loop over sub-boxes
      int ix0,iy0,iz0;
      int idz,np_box,np_2box;
      double x0,y0,z0,xr[3],r2;
      int irl,irt;
      lint jj,this_box;

      ix0=ibox%nside;		// coordinates of current sub-box
      iy0=(ibox%(nside*nside))/nside;
      iz0=ibox/(nside*nside);

      if(boxes[ibox].np>0) {	//the sub-box is not empty
        np_box = boxes[ibox].np;

        for(ii=0;ii<np_box;ii++)  { 	//loop over particles in the sub-box
          x0=(boxes[ibox].pos)[3*ii];
          y0=(boxes[ibox].pos)[3*ii+1];
          z0=(boxes[ibox].pos)[3*ii+2];

          for(jj=ii+1;jj<np_box;jj++) { //loop over pairs; jj=ii+1 to start ensures each pair counted only once
            xr[0]=x0-(boxes[ibox].pos)[3*jj];
            xr[1]=y0-(boxes[ibox].pos)[3*jj+1];
            xr[2]=z0-(boxes[ibox].pos)[3*jj+2];
            r2=xr[0]*xr[0]+xr[1]*xr[1]+xr[2]*xr[2];

            if(r2<R2_MAX) {
              double rl = fabs(xr[2]);	//takes the l-o-s direction to be z-axis!!
              double rt2=r2-rl*rl;
#ifdef _LOGBIN	  //here as a safety check, ideally should not use logarithmic binning for this case
              //as rl or rt can always be zero (whereas r2 should not ever be)
              if(rl>0) {
                irl=(int)(N_LOGINT*(log10(rl)-LOG_R_MAX)+NB_R);
                if((irl<NB_R)&&(irl>=0)) {
                  if((rt2>0)&&(rt2<R2_MAX)) {
                    irt=(int)(N_LOGINT*(0.5*log10(rt2)-LOG_R_MAX)+NB_R);
                    if((irt<NB_R)&&(irt>=0)) {
                      (hthread[irl+NB_R*irt])++;
                    }
                  }
                }
              }
#else //_LOGBIN
              irl=(int)(rl*I_DR);
              if((irl<NB_R)&&(irl>=0)) {
                if(rt2<R2_MAX) {
                  irt=(int)(sqrt(rt2)*I_DR);
                  if((irt<NB_R)&&(irt>=0)) {
                    (hthread[irl+NB_R*irt])++;
                  }
                }
              }
#endif //_LOGBIN
            } //endif r2<R2_MAX
          }	//end for loop jj
        }   //end loop over ii
        //now look for nearby sub-boxes
        for(idz=-index_max;idz<=index_max;idz++) { 
          int idy;
          int iwrapz=0;
          int iz1=iz0+idz;
          if(iz1<0) {
            iz1+=nside;
            iwrapz=1;
          }
          else if(iz1>=nside) {
            iz1-=nside;
            iwrapz=1;
          }
          for(idy=-index_max;idy<=index_max;idy++) {
            int idx;
            int iwrapy=0;
            int iy1=iy0+idy;
            if(iy1<0) {
              iy1+=nside;
              iwrapy=1;
            }
            else if(iy1>=nside) {
              iy1-=nside;
              iwrapy=1;
            }
            for(idx=-index_max;idx<=index_max;idx++) {
              int iwrapx=0;
              int ix1=ix0+idx;
              if(ix1<0) {
                ix1+=nside;
                iwrapx=1;
              }
              else if(ix1>=nside) {
                ix1-=nside;
                iwrapx=1;
              }
              this_box=ix1+nside*(iy1+nside*iz1);		//index of nearby sub-box
              if((this_box>ibox)&&(boxes[this_box].np>0)) {	//only count pairs of sub-boxes once
                np_2box = boxes[this_box].np;
                for(ii=0;ii<np_box;ii++)  {
                  x0=(boxes[ibox].pos)[3*ii];
                  y0=(boxes[ibox].pos)[3*ii+1];
                  z0=(boxes[ibox].pos)[3*ii+2];
                  for(jj=0;jj<np_2box;jj++) {
                    xr[0]=fabs(x0-(boxes[this_box].pos)[3*jj]);
                    xr[1]=fabs(y0-(boxes[this_box].pos)[3*jj+1]);
                    xr[2]=fabs(z0-(boxes[this_box].pos)[3*jj+2]);
                    if(iwrapx) xr[0]=l_box-xr[0];	//check PBC due to wrapping
                    if(iwrapy) xr[1]=l_box-xr[1];
                    if(iwrapz) xr[2]=l_box-xr[2];
                    r2=xr[0]*xr[0]+xr[1]*xr[1]+xr[2]*xr[2];

                    if(r2<R2_MAX) {
                      double rl = fabs(xr[2]);	//takes the l-o-s direction to be z-axis!!
                      double rt2=r2-rl*rl;
#ifdef _LOGBIN	      //here as a safety check, ideally should not use logarithmic binning for this case
                      //as rl or rt can always be zero (whereas r2 should not ever be)
                      if(rl>0) {
                        irl=(int)(N_LOGINT*(log10(rl)-LOG_R_MAX)+NB_R);
                        if((irl<NB_R)&&(irl>=0)) {
                          if((rt2>0)&&(rt2<R2_MAX)) {
                            irt=(int)(N_LOGINT*(0.5*log10(rt2)-LOG_R_MAX)+NB_R);
                            if((irt<NB_R)&&(irt>=0)) {
                              (hthread[irl+NB_R*irt])++;
                            }
                          }
                        }
                      }
#else //_LOGBIN
                      irl=(int)(rl*I_DR);
                      if((irl<NB_R)&&(irl>=0)) {
                        if(rt2<R2_MAX) {
                          irt=(int)(sqrt(rt2)*I_DR);
                          if((irt<NB_R)&&(irt>=0)) {
                            (hthread[irl+NB_R*irt])++;
                          }
                        }
                      }
#endif //_LOGBIN
                    } //endif r2<R2_MAX
                  }	//end for loop jj
                }	//end for loop ii
              }	//end if this_box
            }	//end for idx
          }	//end for idy
        }	//end for idz
      }	//end if boxes[ibox].np>0
    } // end loop over ibox

#pragma omp critical
    {
      for(ii=0;ii<NB_R*NB_R;ii++) //Check bound
        hh[ii]+=hthread[ii]; //Add private histograms to shared one
    } // end omp critical

  } // end omp parallel

}

void cross_3d_ps_boxes(int nside,NeighborBox *boxes1,
    NeighborBox *boxes2,unsigned long long hh[])
{
  //////
  // Cross-correlator for xi(pi,sigma) in the periodic-box case
  double agrid=l_box/nside;
  double r2_max=2./(I_R_MAX*I_R_MAX);
  //double rt2_max=1./(I_R_MAX*I_R_MAX);
  int index_max=(int)(sqrt(r2_max)/agrid)+1;
  int i;

  printf("  Boxes will be correlated up to %d box sizes \n",index_max);

  for(i=0;i<NB_R*NB_R;i++)
    hh[i]=0; //Clear shared histogram

#pragma omp parallel default(none)			\
  shared(index_max,nside,boxes1,boxes2,hh,l_box,agrid)
  {
    lint ii,ibox;
    unsigned long long hthread[NB_R*NB_R]; //Histogram filled by each thread

    for(ii=0;ii<NB_R*NB_R;ii++)
      hthread[ii]=0; //Clear private histogram

#pragma omp for nowait schedule(dynamic)
    for(ibox=0;ibox<nside*nside*nside;ibox++) {	//loop over sub-boxes
      int ix0,iy0,iz0;
      int idz,np1_box,np2_box;
      double x0,y0,z0,xr[3],r2;
      int irt,irl;
      lint jj,this_box;

      ix0=ibox%nside;		// coordinates of current sub-box
      iy0=(ibox%(nside*nside))/nside;
      iz0=ibox/(nside*nside);

      for(idz=-index_max;idz<=index_max;idz++) { //loop over boxes adjacent in z-direction
        int idy;
        int iwrapz=0;
        int iz1=iz0+idz;
        if(iz1<0) {	//check for PBC wraparound
          iz1+=nside;
          iwrapz=1;
        }
        else if(iz1>=nside) {
          iz1-=nside;
          iwrapz=1;
        }
        for(idy=-index_max;idy<=index_max;idy++) { //loop over boxes adjacent in y-direction
          int idx;
          int iwrapy=0;
          int iy1=iy0+idy;
          if(iy1<0) {	//check for PBC wraparound
            iy1+=nside;
            iwrapy=1;
          }
          else if(iy1>=nside) {
            iy1-=nside;
            iwrapy=1;
          }
          for(idx=-index_max;idx<=index_max;idx++) { //loop over boxes adjacent in x-direction
            int iwrapx=0;
            int ix1=ix0+idx;
            if(ix1<0) {	//check for PBC wraparound
              ix1+=nside;
              iwrapx=1;
            }
            else if(ix1>=nside) {
              ix1-=nside;
              iwrapx=1;
            }
            this_box=ix1+nside*(iy1+nside*iz1);
            np1_box=boxes1[ibox].np;
            np2_box=boxes2[this_box].np;
            for(ii=0;ii<np1_box;ii++) { //otherwise, find distances to particle pairs in this box
              x0=(boxes1[ibox].pos)[3*ii];
              y0=(boxes1[ibox].pos)[3*ii+1];
              z0=(boxes1[ibox].pos)[3*ii+2];
              for(jj=0;jj<np2_box;jj++) {
                xr[0]=fabs(x0-(boxes2[this_box].pos)[3*jj]);	//calculate distance between particles
                xr[1]=fabs(y0-(boxes2[this_box].pos)[3*jj+1]);
                xr[2]=fabs(z0-(boxes2[this_box].pos)[3*jj+2]);
                if(iwrapx) xr[0]=l_box-xr[0];	//check PBC due to wrapping
                if(iwrapy) xr[1]=l_box-xr[1];
                if(iwrapz) xr[2]=l_box-xr[2];
                r2=xr[0]*xr[0]+xr[1]*xr[1]+xr[2]*xr[2];

                if(r2<R2_MAX) {
                  double rl = fabs(xr[2]);	//takes the l-o-s direction to be z-axis!!
                  double rt2=r2-rl*rl;
#ifdef _LOGBIN	  //here as a safety check, ideally should not use logarithmic binning for this case
                  //as rl or rt can always be zero (whereas r2 should not ever be)
                  if(rl>0) {
                    irl=(int)(N_LOGINT*(log10(rl)-LOG_R_MAX)+NB_R);
                    if((irl<NB_R)&&(irl>=0)) {
                      if((rt2>0)&&(rt2<R2_MAX)) {
                        irt=(int)(N_LOGINT*(0.5*log10(rt2)-LOG_R_MAX)+NB_R);
                        if((irt<NB_R)&&(irt>=0)) {
                          (hthread[irl+NB_R*irt])++;
                        }
                      }
                    }
                  }
#else //_LOGBIN
                  irl=(int)(rl*I_DR);
                  if((irl<NB_R)&&(irl>=0)) {
                    if(rt2<R2_MAX) {
                      irt=(int)(sqrt(rt2)*I_DR);
                      if((irt<NB_R)&&(irt>=0)) {
                        (hthread[irl+NB_R*irt])++;
                      }
                    }
                  }
#endif //_LOGBIN
                } //endif r2<R2_MAX
              }	//end for loop jj
            }	//end for loop ii
          }	//end for loop idx
        }	//end for loop idy
      }	//end for loop idz
    } // end loop over ibox

#pragma omp critical
    {
      for(ii=0;ii<NB_R*NB_R;ii++) //Check bound
        hh[ii]+=hthread[ii]; //Add private histograms to shared one
    } // end pragma omp critical

  } // end pragma omp parallel

}

void auto_3d_rmu_boxes(int nside,NeighborBox *boxes,
    unsigned long long hh[])
{
  //////
  // Correlator for xi(r,mu) in the periodic-box case
  double agrid=l_box/nside;
  int index_max=(int)(sqrt(R2_MAX)/agrid)+1;
  int i;

  printf("  Boxes will be correlated up to %d box sizes \n",index_max);

  for(i=0;i<NB_R*NB_mu;i++)
    hh[i]=0; //Clear shared histogram

#pragma omp parallel default(none)			\
  shared(index_max,nside,boxes,hh,l_box,agrid)
  {
    lint ii,ibox;
    unsigned long long hthread[NB_R*NB_mu]; //Histogram filled by each thread

    for(ii=0;ii<NB_R*NB_mu;ii++)
      hthread[ii]=0; //Clear private histogram

#pragma omp for nowait schedule(dynamic)
    for(ibox=0;ibox<nside*nside*nside;ibox++) { //loop over sub-boxes
      int ix0,iy0,iz0;
      int idz,np_box,np_2box;
      double x0,y0,z0,xr[3],r2;
      int ir,imu;
      lint jj,this_box;

      ix0=ibox%nside;		// coordinates of current sub-box
      iy0=(ibox%(nside*nside))/nside;
      iz0=ibox/(nside*nside);

      if(boxes[ibox].np>0) {	//the sub-box is not empty
        np_box = boxes[ibox].np;

        for(ii=0;ii<np_box;ii++)  { 	//loop over particles in the sub-box
          x0=(boxes[ibox].pos)[3*ii];
          y0=(boxes[ibox].pos)[3*ii+1];
          z0=(boxes[ibox].pos)[3*ii+2];

          for(jj=ii+1;jj<np_box;jj++) { //loop over pairs; jj=ii+1 to start ensures each pair counted only once
            xr[0]=x0-(boxes[ibox].pos)[3*jj];
            xr[1]=y0-(boxes[ibox].pos)[3*jj+1];
            xr[2]=z0-(boxes[ibox].pos)[3*jj+2];
            r2=xr[0]*xr[0]+xr[1]*xr[1]+xr[2]*xr[2];

            if(r2<R2_MAX) {
              ir = (int)(sqrt(r2)*I_DR);
              if((ir<NB_R)&&(ir>=0)) {
                if(r2==0) imu=0;
                else {
                  double mu = fabs(xr[2])/sqrt(r2);	//takes the l-o-s direction to be z-axis!!
                  imu = (int)(mu*NB_mu);
                }
                if((imu<NB_mu)&&(imu>=0)) {
                  (hthread[imu+NB_mu*ir])++;
                }
              }
            } //endif r2<R2_MAX
          }	//end for loop jj
        }	//end loop over ii

        //now look for nearby sub-boxes
        for(idz=-index_max;idz<=index_max;idz++) {
          int idy;
          //int iwrapz=0;
          int iz1=iz0+idz;
          if(iz1<0) {
            iz1+=nside;
            //iwrapz=1;
          }
          else if(iz1>=nside) {
            iz1-=nside;
            //iwrapz=1;
          }
          for(idy=-index_max;idy<=index_max;idy++) {
            int idx;
            //int iwrapy=0;
            int iy1=iy0+idy;
            if(iy1<0) {
              iy1+=nside;
              //iwrapy=1;
            }
            else if(iy1>=nside) {
              iy1-=nside;
              //iwrapy=1;
            }
            for(idx=-index_max;idx<=index_max;idx++) {
              //int iwrapx=0;
              int ix1=ix0+idx;
              if(ix1<0) {
                ix1+=nside;
                //iwrapx=1;
              }
              else if(ix1>=nside) {
                ix1-=nside;
                //iwrapx=1;
              }
              this_box=ix1+nside*(iy1+nside*iz1);		//index of nearby sub-box
              if((this_box>ibox)&&(boxes[this_box].np>0)) {	//only count pairs of sub-boxes once
                np_2box = boxes[this_box].np;
                for(ii=0;ii<np_box;ii++)  {
                  x0=(boxes[ibox].pos)[3*ii];
                  y0=(boxes[ibox].pos)[3*ii+1];
                  z0=(boxes[ibox].pos)[3*ii+2];
                  for(jj=0;jj<np_2box;jj++) {
                    xr[0] = x0 - (boxes[this_box].pos)[3*jj];
                    xr[1] = y0 - (boxes[this_box].pos)[3*jj+1];
                    xr[2] = z0 - (boxes[this_box].pos)[3*jj+2];
                    if(fabs(xr[0])>l_box/2) xr[0] -= l_box*xr[0]/fabs(xr[0]);
                    if(fabs(xr[1])>l_box/2) xr[1] -= l_box*xr[1]/fabs(xr[1]);
                    if(fabs(xr[2])>l_box/2) xr[2] -= l_box*xr[2]/fabs(xr[2]);
                    r2=xr[0]*xr[0]+xr[1]*xr[1]+xr[2]*xr[2];

                    if(r2<R2_MAX) {
                      ir = (int)(sqrt(r2)*I_DR);
                      if((ir<NB_R)&&(ir>=0)) {
                        if(r2==0) imu=0;
                        else {
                          double mu = fabs(xr[2])/sqrt(r2);	//takes the l-o-s direction to be z-axis!!
                          imu = (int)(mu*NB_mu);
                        }
                        if((imu<NB_mu)&&(imu>=0)) {
                          (hthread[imu+NB_mu*ir])++;
                        }
                      }
                    } //endif r2<R2_MAX
                  }	//end for loop jj
                }	//end for loop ii
              }	//end if this_box
            }	//end for idx
          }	//end for idy
        }	//end for idz
      }	//end if boxes[ibox].np>0
    } // end omp for over ibox

#pragma omp critical
    {
      for(ii=0;ii<NB_R*NB_mu;ii++) //Check bound
        hh[ii]+=hthread[ii]; //Add private histograms to shared one
    } // end omp critical

  } // end omp parallel

}

void cross_3d_rmu_boxes(int nside,NeighborBox *boxes1,
    NeighborBox *boxes2,unsigned long long hh[])
{
  //////
  // Cross-correlator for xi(r,mu) in the periodic-box case
  double agrid=l_box/nside;
  //double r2_max=2./(I_R_MAX*I_R_MAX);
  int index_max=(int)(sqrt(R2_MAX)/agrid)+1;
  int i;

  printf("  Boxes will be correlated up to %d box sizes \n",index_max);

  for(i=0;i<NB_R*NB_mu;i++)
    hh[i]=0; //Clear shared histogram

#pragma omp parallel default(none)			\
  shared(index_max,nside,boxes1,boxes2,hh,l_box)
  {
    lint ii,ibox;
    unsigned long long hthread[NB_R*NB_mu]; //Histogram filled by each thread

    for(ii=0;ii<NB_R*NB_mu;ii++)
      hthread[ii]=0; //Clear private histogram

#pragma omp for nowait schedule(dynamic)
    for(ibox=0;ibox<nside*nside*nside;ibox++) {	//loop over sub-boxes
      int ix0,iy0,iz0;
      int idz,np1_box,np2_box;
      double x0,y0,z0,xr[3],r2;
      int ir,imu;
      lint jj,this_box;

      ix0=ibox%nside;		// coordinates of current sub-box
      iy0=(ibox%(nside*nside))/nside;
      iz0=ibox/(nside*nside);

      for(idz=-index_max;idz<=index_max;idz++) { //loop over boxes adjacent in z-direction
        int idy;
        //int iwrapz=0;
        int iz1=iz0+idz;
        if(iz1<0) {	//check for PBC wraparound
          iz1+=nside;
          //iwrapz=1;
        }
        else if(iz1>=nside) {
          iz1-=nside;
          //iwrapz=1;
        }
        for(idy=-index_max;idy<=index_max;idy++) { //loop over boxes adjacent in y-direction
          int idx;
          //int iwrapy=0;
          int iy1=iy0+idy;
          if(iy1<0) {	//check for PBC wraparound
            iy1+=nside;
            //iwrapy=1;
          }
          else if(iy1>=nside) {
            iy1-=nside;
            //iwrapy=1;
          }
          for(idx=-index_max;idx<=index_max;idx++) { //loop over boxes adjacent in x-direction
            //int iwrapx=0;
            int ix1=ix0+idx;
            if(ix1<0) {	//check for PBC wraparound
              ix1+=nside;
              //iwrapx=1;
            }
            else if(ix1>=nside) {
              ix1-=nside;
              //iwrapx=1;
            }
            this_box=ix1+nside*(iy1+nside*iz1);
            np1_box=boxes1[ibox].np;
            np2_box=boxes2[this_box].np;
            for(ii=0;ii<np1_box;ii++) { //otherwise, find distances to particle pairs in this box
              x0=(boxes1[ibox].pos)[3*ii];
              y0=(boxes1[ibox].pos)[3*ii+1];
              z0=(boxes1[ibox].pos)[3*ii+2];
              for(jj=0;jj<np2_box;jj++) {
                xr[0] = x0 - (boxes2[this_box].pos)[3*jj];
                xr[1] = y0 - (boxes2[this_box].pos)[3*jj+1];
                xr[2] = z0 - (boxes2[this_box].pos)[3*jj+2];
                if(fabs(xr[0])>l_box/2) xr[0] -= l_box*xr[0]/fabs(xr[0]);
                if(fabs(xr[1])>l_box/2) xr[1] -= l_box*xr[1]/fabs(xr[1]);
                if(fabs(xr[2])>l_box/2) xr[2] -= l_box*xr[2]/fabs(xr[2]);
                r2=xr[0]*xr[0]+xr[1]*xr[1]+xr[2]*xr[2];

                if(r2<R2_MAX) {
                  ir = (int)(sqrt(r2)*I_DR);
                  if((ir<NB_R)&&(ir>=0)) {
                    if(r2==0) imu=0;
                    else {
                      double mu = fabs(xr[2])/sqrt(r2);	//takes the l-o-s direction to be z-axis!!
                      imu = (int)(mu*NB_mu);
                    }
                    if((imu<NB_mu)&&(imu>=0)) {
                      (hthread[imu+NB_mu*ir])++;
                    }
                  }
                } //endif r2<R2_MAX
              }	//end for loop jj
            }	//end for loop ii
          }	//end for loop idx
        }	//end for loop idy
      }	//end for loop idz
    } // end pragma omp for (ibox)

#pragma omp critical
    {
      for(ii=0;ii<NB_R*NB_mu;ii++) //Check bound
        hh[ii]+=hthread[ii]; //Add private histograms to shared one
    } // end pragma omp critical
  } // end pragma omp parallel
}

