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

#ifndef _CUTE_DEFINE_
#define _CUTE_DEFINE_

#ifdef _LONGIDS
typedef long lint;
#else //_LONGIDS
typedef int lint;
#endif //_LONGIDS

extern char fnameData[256];
extern char fnameData2[256];
extern char fnameRand[256];
extern char fnameOut[256];
extern int input_format;
extern int use_tree;
extern int max_tree_order;
extern int max_tree_nparts;
extern int use_pm;
extern int use_randoms;
extern int reuse_randoms;
extern int do_CCF;
extern int corr_type;
extern lint n_objects;
extern float l_box;
extern float l_box_half;
extern int n_grid;
/*                MACROS            */
// Other possible macros
//_DEBUG, _VERBOSE, _TRUE_ACOS, _LOGBIN

#define MIN(a, b) (((a) < (b)) ? (a) : (b)) //Minimum of two numbers
#define MAX(a, b) (((a) > (b)) ? (a) : (b)) //Maximum of two numbers
#define ABS(a)   (((a) < 0) ? -(a) : (a)) //Absolute value
#define CLAMP(x, low, high)  (((x) > (high)) ? (high) : (((x) < (low)) ? (low) : (x))) //min(max(a,low),high)

#ifndef N_LOGINT
#define N_LOGINT 20 //# bins per decade for logarithmic binning
#endif

//Monopole
#ifndef I_R_MAX
#define I_R_MAX 0.005 //1/r_max
#endif
#ifndef LOG_R_MAX
#define LOG_R_MAX 2.30103 //log10(r_max)
#endif
#ifndef NB_R
#define NB_R 256 //# bins in r for monopole correlation
#endif

//2D binning
#ifndef NB_mu
#define NB_mu 100
#endif

typedef struct {
  lint np;          //#objects in the catalog
  double *pos;
} Catalog;         //Catalog (double precision)

typedef struct branch {
  float x_lo[3];
  float x_hi[3];
  char leaf;
  lint np;
  void *sons; //sons
} branch; //Tree node/branch

typedef struct {
  int np;
  double *pos;
} NeighborBox; //Neighbor box

extern int cute_verbose;

#ifdef _CUTE_AS_PYTHON_MODULE

typedef struct {
  int nx, ny, nz;
  double *x, *y, *z, *corr, 
         *D1D1, *D1D2, *D1R1, *D1R2,
         *D2D2, *D2R1, *D2R2, 
         *R1R1, *R1R2, 
         *R2R2;
} Result;

extern Result *global_result;
extern Catalog *global_galaxy_catalog;
extern Catalog *global_galaxy_catalog2;
extern Catalog *global_random_catalog;

#endif

#endif //_CUTE_DEFINE_
