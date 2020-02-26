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

#ifndef _CUTE_COMMON_
#define _CUTE_COMMON_

void timer(int i);

lint linecount(FILE *f);

int optimal_nside(double lb,double rmax,lint np);

void free_catalog(Catalog *cat);

void error_mem_out(void);

void error_open_file(char *fname);

void error_read_line(char *fname,lint nlin);

#ifdef _DEBUG
void write_cat(Catalog cat,char *fn);

void write_grid(double *grid,char *fn);

void write_tree(branch *tree,char *fn);
#endif //_DEBUG

#ifdef _CUTE_AS_PYTHON_MODULE

Result *make_empty_result_struct();
void free_result_struct(Result *res);
void set_result(Result *res, int i, double x, double corr,
    double D1D1, double D1D2, double D1R1, double D1R2,
    double D2D2, double D2R1, double D2R2,
    double R1R1, double R1R2,
    double R2R2);
void set_result_2d(Result *res, int i, int j, int ind, double x, double y, double corr,
    double D1D1, double D1D2, double D1R1, double D1R2,
    double D2D2, double D2R1, double D2R2,
    double R1R1, double R1R2,
    double R2R2);
void set_result_3d(Result *res,
    int i, int j, int k, int ind,
    double x, double y, double z, double corr,
    double D1D1, double D1D2, double D1R1, double D1R2,
    double D2D2, double D2R1, double D2R2,
    double R1R1, double R1R2,
    double R2R2);

Catalog *create_catalog_from_numpy(int nx, double *x, int ny, double *y, int nz, double *z);

#endif

#endif //_CUTE_COMMON_
