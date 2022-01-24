/*
 *  This file is part of t8dg.
 *  t8dg is a discontinuous galerkin solver on adaptive meshes
 *  that uses the adaptive meshing library t8code to distribute
 *  the workload equally among parallel processes to achieve good scaling.
 *
 *  Copyright (C) 2020 the developers
 *
 *  t8dg is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  t8dg is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with t8dg.  If not, see <http://www.gnu.org/licenses/>.
 */

/** @file t8dg_square_matrix.h */
#ifndef SRC_T8DG_DMATRIX_H_
#define SRC_T8DG_DMATRIX_H_
#include "t8dg.h"

T8DG_EXTERN_C_BEGIN ();

typedef struct t8dg_dmatrix t8dg_dmatrix_t;

double              t8dg_dmatrix_at (const t8dg_dmatrix_t * dmatrix, const int irow, const int icolumn);

void                t8dg_dmatrix_set_at (t8dg_dmatrix_t * matrix, const int irow, const int icolumn, const double value);

void                t8dg_dmatrix_mult_sc_array (const t8dg_dmatrix_t * A, const sc_array_t * x, sc_array_t * b);

void                t8dg_dmatrix_transpose_mult_sc_array (const t8dg_dmatrix_t * A, const sc_array_t * x, sc_array_t * b);

void                t8dg_dmatrix_scale_row (t8dg_dmatrix_t * matrix, const int irow, const double alpha);

t8dg_dmatrix_t     *t8dg_dmatrix_new (int nrows, int ncolumns);

t8dg_dmatrix_t     *t8dg_dmatrix_new_zero (int nrows, int ncolumns);

t8dg_dmatrix_t     *t8dg_dmatrix_new_data (int nrows, int ncolumns, double *data);

void                t8dg_dmatrix_destroy (t8dg_dmatrix_t ** pmatrix);

typedef double      t8dg_square_3D_matrix_t[DIM3][DIM3];

void                t8dg_square3D_matrix_invert_sub_matrix (t8dg_square_3D_matrix_t matrix_invers, t8dg_square_3D_matrix_t matrix, int dim);

void                t8dg_square3D_matrix_determinant_sub_matrix (double *det, t8dg_square_3D_matrix_t matrix, int dim);

void                t8dg_square3D_matrix_scale (t8dg_square_3D_matrix_t matrix, double scaling_factor, int dim);

void                t8dg_square3D_matrix_apply_rotation_reflection_matrix_to_matrix (t8dg_square_3D_matrix_t
                                                                                     matrix_invers,
                                                                                     t8dg_square_3D_matrix_t matrix,
                                                                                     int idx_rotation_reflection, int dim);

T8DG_EXTERN_C_END ();

#endif /* SRC_T8DG_DMATRIX_H_ */
