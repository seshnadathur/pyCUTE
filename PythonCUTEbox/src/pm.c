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
//             Routines to create the particle mesh                  //
/*********************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <complex.h>
#include "define.h"
#include "common.h"

static lint wrap_int(lint n,int ngrid)
{
  //////
  // Returns n mod(ngrid)
  if(n<0)
    return wrap_int(n+ngrid,ngrid);
  else if(n>=ngrid)
    return wrap_int(n-ngrid,ngrid);
  else
    return n;
}

double *read_grid(void)
{
  //////
  // Reads a grid[n_grid*n_grid*n_grid] with    
  // the pre-calculated density field from 
  // given file and returns it.
  printf("  Reading grid density data ...\n");
  lint n_grid_tot = n_grid*((lint)(n_grid*n_grid));
  double *grid=(double *)calloc(n_grid_tot,sizeof(double));
  if(grid==NULL) error_mem_out();
  FILE *fp;

  //printf("  n_grid = %d, n_grid_tot = %ld\n",n_grid,n_grid_tot);
  fp = fopen(fnameData,"r");
  //printf("  Opened file\n");
  fread(grid,sizeof(double),n_grid_tot,fp);
  fclose(fp);

  // the file contains overdensity delta
  // couple of checks ...
  //printf("  grid[0] = %lf\n",grid[0]);
  //printf("  grid[%d] = %lf\n",n_grid-1,grid[n_grid-1]);
  //printf("  grid[%ld] = %lf\n",n_grid_tot-1,grid[n_grid_tot-1]);

  return grid;

}

double *resize_grid(double *grid,int new_n_grid)
{
  //////
  // Takes overdensity Delta values on a grid 
  // and resizes the grid to lower resolution 
  // to save computational time later

  lint n_grid_tot = n_grid*((lint)(n_grid*n_grid));
  lint new_n_grid_tot = new_n_grid*((lint)(new_n_grid*new_n_grid));
  lint i;
  double agrid = l_box/n_grid;
  double new_agrid = l_box/new_n_grid;
  double *new_grid=(double *)calloc(new_n_grid_tot,sizeof(double));
  if(new_grid==NULL) error_mem_out();
  double ivcell=1./(agrid*agrid*agrid);

  printf("  Resizing density grid resolution from %lf to %lf ...\n",agrid,new_agrid);

  for(i=0;i<new_n_grid_tot;i++) {
    new_grid[i]=0; 
  }

///   
#pragma omp parallel default(none)		\
  shared(n_grid,n_grid_tot,ivcell,agrid,new_agrid,new_n_grid,grid,new_grid)
  {
    lint ii;
#pragma omp for nowait
    for(ii=0;ii<n_grid_tot;ii++) {	// loop over each cell of old grid
	lint indx,indy,indz;
      	lint i0x,i0y,i0z;
      	double xp,yp,zp;

      	indx = ii/(n_grid*n_grid);	// x,y,z indices of old grid cell
      	indy = ii%(n_grid*n_grid)/n_grid;
      	indz = ii%n_grid;

      	xp = (indx+0.5)*agrid;		// x,y,z coordinates
       	yp = (indy+0.5)*agrid;
      	zp = (indz+0.5)*agrid;
      
      	i0x=(lint)(xp/new_agrid);	// indices of cell on new grid
      	i0y=(lint)(yp/new_agrid);
      	i0z=(lint)(zp/new_agrid);

	new_grid[i0z+new_n_grid*(i0y+new_n_grid*i0x)] += (1+grid[ii])/ivcell;

    } //end omp for
  } //end omp parallel

#pragma omp parallel default(none)		\
  shared(new_grid,new_n_grid_tot,new_agrid)
  {
    lint jj;

#pragma omp for nowait
    for(jj=0;jj<new_n_grid_tot;jj++)
	new_grid[jj] = new_grid[jj]/(new_agrid*new_agrid*new_agrid) - 1.;
  }

  return new_grid;
///
   
}

double *pos_2_ngp(Catalog cat)
{
  //////
  // Returns a grid[n_grid*n_grid*n_grid] with    
  // the density field calculated using the NGP
  // (nearest-grid-point) technique from a set
  // of np particle positions **pos.
  printf("  Calculating NGP...\n");
  lint n_grid_tot=n_grid*((lint)(n_grid*n_grid));
  double *grid=(double *)calloc(n_grid_tot,sizeof(double));
  if(grid==NULL) error_mem_out();

#pragma omp parallel default(none)		\
  shared(cat,n_grid,n_grid_tot,grid,l_box)
  {
    lint ii;
    double agrid=l_box/n_grid;
#pragma omp for nowait
    for(ii=0;ii<cat.np;ii++) {
      lint i0x,i0y,i0z;
      double xp=CLAMP(cat.pos[3*ii],0,l_box);
      double yp=CLAMP(cat.pos[3*ii+1],0,l_box);
      double zp=CLAMP(cat.pos[3*ii+2],0,l_box);
      i0x=(lint)(xp/agrid);
      i0y=(lint)(yp/agrid);
      i0z=(lint)(zp/agrid);

#pragma omp atomic
      grid[i0x+n_grid*(i0y+n_grid*i0z)]++;
    } //end omp for
  } //end omp parallel
  
  //Substract background
#pragma omp parallel default(none)		\
  shared(grid,n_grid,n_grid_tot,cat)
  {
    lint ii;
    double n_avg=(double)(cat.np)/n_grid_tot;
    double i_n_avg=1/n_avg;
    
#pragma omp for
    for(ii=0;ii<n_grid_tot;ii++)
      grid[ii]=grid[ii]*i_n_avg-1;
  }

  return grid;
}

double *pos_2_cic(Catalog cat)
{
  //////
  // Returns a grid[n_grid*n_grid*n_grid] with    
  // the density field calculated using the CIC
  // (cloud-in-cell) technique from a set of np
  // particle positions **pos.
  printf("  Calculating CIC...\n");
  lint n_grid_tot=n_grid*((lint)(n_grid*n_grid));
  double *grid=(double *)calloc(n_grid_tot,sizeof(double));
  if(grid==NULL) error_mem_out();

#pragma omp parallel default(none)		\
  shared(cat,n_grid,n_grid_tot,grid,l_box)
  {
    lint ii;
    double agrid=l_box/n_grid;
    double ivcell=1./(agrid*agrid*agrid);
#pragma omp for nowait
    for(ii=0;ii<cat.np;ii++) {
      lint i0x,i0y,i0z;
      lint i1x,i1y,i1z;
      double a1x,a1y,a1z;
      double a0x,a0y,a0z;
      double xp=CLAMP(cat.pos[3*ii],0,l_box);
      double yp=CLAMP(cat.pos[3*ii+1],0,l_box);
      double zp=CLAMP(cat.pos[3*ii+2],0,l_box);
      i0x=(lint)(xp/agrid);
      i0y=(lint)(yp/agrid);
      i0z=(lint)(zp/agrid);

      a1x=xp-(i0x+0.5)*agrid;
      a1y=yp-(i0y+0.5)*agrid;
      a1z=zp-(i0z+0.5)*agrid;

      if(a1x<0) {
        a1x=-a1x;
        i1x=wrap_int(i0x-1,n_grid);
      }
      else
	i1x=wrap_int(i0x+1,n_grid);
      a0x=agrid-a1x;

      if(a1y<0) {
        a1y=-a1y;
        i1y=wrap_int(i0y-1,n_grid);
      }
      else
	i1y=wrap_int(i0y+1,n_grid);
      a0y=agrid-a1y;

      if(a1z<0) {
        a1z=-a1z;
        i1z=wrap_int(i0z-1,n_grid);
      }
      else
	i1z=wrap_int(i0z+1,n_grid);
      a0z=agrid-a1z;

#pragma omp atomic
      grid[i0x+n_grid*(i0y+n_grid*i0z)]+=
	a0x*a0y*a0z*ivcell;
#pragma omp atomic
      grid[i1x+n_grid*(i0y+n_grid*i0z)]+=
	a1x*a0y*a0z*ivcell;
#pragma omp atomic
      grid[i0x+n_grid*(i1y+n_grid*i0z)]+=
	a0x*a1y*a0z*ivcell;
#pragma omp atomic
      grid[i1x+n_grid*(i1y+n_grid*i0z)]+=
	a1x*a1y*a0z*ivcell;
#pragma omp atomic
      grid[i0x+n_grid*(i0y+n_grid*i1z)]+=
	a0x*a0y*a1z*ivcell;
#pragma omp atomic
      grid[i1x+n_grid*(i0y+n_grid*i1z)]+=
	a1x*a0y*a1z*ivcell;
#pragma omp atomic
      grid[i0x+n_grid*(i1y+n_grid*i1z)]+=
	a0x*a1y*a1z*ivcell;
#pragma omp atomic
      grid[i1x+n_grid*(i1y+n_grid*i1z)]+=
	a1x*a1y*a1z*ivcell;
    } //end omp for
  } //end omp parallel

  //Substract background
#pragma omp parallel default(none)		\
  shared(grid,n_grid,n_grid_tot,cat)
  {
    lint ii;
    double n_avg=(double)(cat.np)/n_grid_tot;
    double i_n_avg=1/n_avg;

#pragma omp for
    for(ii=0;ii<n_grid_tot;ii++)
      grid[ii]=grid[ii]*i_n_avg-1;
  }

  return grid;
}

double *pos_2_tsc(Catalog cat)
{
  //////
  // Returns a grid[n_grid*n_grid*n_grid] with    
  // the density field calculated using the TSC
  // (triangular-shaped-cloud) technique from a
  // set of np particle positions **pos.
  printf("  Calculating TSC...\n");
  lint n_grid_tot=n_grid*((lint)(n_grid*n_grid));
  double *grid=(double *)calloc(n_grid_tot,sizeof(double));
  if(grid==NULL) error_mem_out();

#pragma omp parallel default(none)		\
  shared(cat,n_grid,n_grid_tot,grid,l_box)
  {
    lint ii;
    double iagrid=n_grid/l_box;
    double agrid=l_box/n_grid;
#pragma omp for nowait
    for(ii=0;ii<cat.np;ii++) {
      lint i0x,i0y,i0z;
      lint ipx,ipy,ipz;
      lint imx,imy,imz;
      double a0x,a0y,a0z;
      double apx,apy,apz;
      double amx,amy,amz;
      double xp=CLAMP(cat.pos[3*ii],0,l_box);
      double yp=CLAMP(cat.pos[3*ii+1],0,l_box);
      double zp=CLAMP(cat.pos[3*ii+2],0,l_box);
      i0x=(lint)(xp/agrid);
      i0y=(lint)(yp/agrid);
      i0z=(lint)(zp/agrid);
      a0x=xp*iagrid-(i0x+0.5);
      a0y=yp*iagrid-(i0y+0.5);
      a0z=zp*iagrid-(i0z+0.5);

      amx=0.5*(0.5-a0x)*(0.5-a0x);
      amy=0.5*(0.5-a0y)*(0.5-a0y);
      amz=0.5*(0.5-a0z)*(0.5-a0z);
      apx=0.5*(0.5+a0x)*(0.5+a0x);
      apy=0.5*(0.5+a0y)*(0.5+a0y);
      apz=0.5*(0.5+a0z)*(0.5+a0z);
      a0x=0.75-a0x*a0x;
      a0y=0.75-a0y*a0y;
      a0z=0.75-a0z*a0z;

      ipx=wrap_int(i0x+1,n_grid);
      ipy=wrap_int(i0y+1,n_grid);
      ipz=wrap_int(i0z+1,n_grid);
      imx=wrap_int(i0x-1,n_grid);
      imy=wrap_int(i0y-1,n_grid);
      imz=wrap_int(i0z-1,n_grid);

#pragma omp atomic
      grid[imx+n_grid*(imy+n_grid*imz)]+=
	amx*amy*amz;
#pragma omp atomic
      grid[i0x+n_grid*(imy+n_grid*imz)]+=
	a0x*amy*amz;
#pragma omp atomic
      grid[ipx+n_grid*(imy+n_grid*imz)]+=
	apx*amy*amz;
#pragma omp atomic
      grid[imx+n_grid*(i0y+n_grid*imz)]+=
	amx*a0y*amz;
#pragma omp atomic
      grid[i0x+n_grid*(i0y+n_grid*imz)]+=
	a0x*a0y*amz;
#pragma omp atomic
      grid[ipx+n_grid*(i0y+n_grid*imz)]+=
	apx*a0y*amz;
#pragma omp atomic
      grid[imx+n_grid*(ipy+n_grid*imz)]+=
	amx*apy*amz;
#pragma omp atomic
      grid[i0x+n_grid*(ipy+n_grid*imz)]+=
	a0x*apy*amz;
#pragma omp atomic
      grid[ipx+n_grid*(ipy+n_grid*imz)]+=
	apx*apy*amz;
#pragma omp atomic
      grid[imx+n_grid*(imy+n_grid*i0z)]+=
	amx*amy*a0z;
#pragma omp atomic
      grid[i0x+n_grid*(imy+n_grid*i0z)]+=
	a0x*amy*a0z;
#pragma omp atomic
      grid[ipx+n_grid*(imy+n_grid*i0z)]+=
	apx*amy*a0z;
#pragma omp atomic
      grid[imx+n_grid*(i0y+n_grid*i0z)]+=
	amx*a0y*a0z;
#pragma omp atomic
      grid[i0x+n_grid*(i0y+n_grid*i0z)]+=
	a0x*a0y*a0z;
#pragma omp atomic
      grid[ipx+n_grid*(i0y+n_grid*i0z)]+=
	apx*a0y*a0z;
#pragma omp atomic
      grid[imx+n_grid*(ipy+n_grid*i0z)]+=
	amx*apy*a0z;
#pragma omp atomic
      grid[i0x+n_grid*(ipy+n_grid*i0z)]+=
	a0x*apy*a0z;
#pragma omp atomic
      grid[ipx+n_grid*(ipy+n_grid*i0z)]+=
	apx*apy*a0z;
#pragma omp atomic
      grid[imx+n_grid*(imy+n_grid*ipz)]+=
	amx*amy*apz;
#pragma omp atomic
      grid[i0x+n_grid*(imy+n_grid*ipz)]+=
	a0x*amy*apz;
#pragma omp atomic
      grid[ipx+n_grid*(imy+n_grid*ipz)]+=
	apx*amy*apz;
#pragma omp atomic
      grid[imx+n_grid*(i0y+n_grid*ipz)]+=
	amx*a0y*apz;
#pragma omp atomic
      grid[i0x+n_grid*(i0y+n_grid*ipz)]+=
	a0x*a0y*apz;
#pragma omp atomic
      grid[ipx+n_grid*(i0y+n_grid*ipz)]+=
	apx*a0y*apz;
#pragma omp atomic
      grid[imx+n_grid*(ipy+n_grid*ipz)]+=
	amx*apy*apz;
#pragma omp atomic
      grid[i0x+n_grid*(ipy+n_grid*ipz)]+=
	a0x*apy*apz;
#pragma omp atomic
      grid[ipx+n_grid*(ipy+n_grid*ipz)]+=
	apx*apy*apz;
    } //end omp for
  } //end omp parallel

  //Substract background
#pragma omp parallel default(none)		\
  shared(grid,n_grid,n_grid_tot,cat)
  {
    lint ii;
    double n_avg=(double)(cat.np)/n_grid_tot;
    double i_n_avg=1/n_avg;
    
#pragma omp for
    for(ii=0;ii<n_grid_tot;ii++)
      grid[ii]=grid[ii]*i_n_avg-1;
  }

  return grid;
}
