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

#ifndef _CUTE_IO_BOX_
#define _CUTE_IO_BOX_

void make_CF(unsigned long long DD[],int nD,
	     double corr[],double ercorr[]);

void make_CF_double(unsigned long long DD[],int nD,
	     double corr[],double ercorr[]);

void make_3d_CCF_w_rand(unsigned long long D1D2[],unsigned long long D1R[],unsigned long long D2R[],
	     unsigned long long RR[],int nD1,int nD2,int nR,double corr[],double ercorr[]);

void make_3d_CF_w_rand(unsigned long long DD[],unsigned long long DR[],
	     unsigned long long RR[],int nD,int nR,double corr[],double ercorr[]);

void make_3d_CCF(unsigned long long D1D2[],int nD1,int nD2,double corr[],double ercorr[]);

void make_3d_CF(unsigned long long DD[],int nD,double corr[],double ercorr[]);

void read_run_params(char *fname);

Catalog *read_catalog(char *fname,lint *np);

#endif //_CUTE_IO_BOX_
