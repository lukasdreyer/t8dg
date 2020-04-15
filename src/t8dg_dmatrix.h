/*
 * t8dg_square_3D_matrix.h
 *
 *  Created on: Mar 21, 2020
 *      Author: lukas
 */

/** @file t8dg_square_matrix.h */
#ifndef SRC_T8DG_DMATRIX_H_
#define SRC_T8DG_DMATRIX_H_
#include "t8dg.h"

typedef struct t8dg_dmatrix t8dg_dmatrix_t;

double              t8dg_dmatrix_at (t8dg_dmatrix_t * dmatrix, int irow, int icolumn);

void                t8dg_dmatrix_set_at (t8dg_dmatrix_t * matrix, int irow, int icolumn, double value);

void                t8dg_dmatrix_mult_sc_array (t8dg_dmatrix_t * A, sc_array_t * x, sc_array_t * b);

void                t8dg_dmatrix_transpose_mult_sc_array (t8dg_dmatrix_t * A, sc_array_t * x, sc_array_t * b);

t8dg_dmatrix_t     *t8dg_dmatrix_new (int nrows, int ncolumns);

void                t8dg_dmatrix_destroy (t8dg_dmatrix_t ** pmatrix);

typedef double      t8dg_square_3D_matrix_t[DIM3][DIM3];

void                t8dg_square3D_matrix_invert_sub_matrix (t8dg_square_3D_matrix_t matrix_invers, t8dg_square_3D_matrix_t matrix, int dim);

void                t8dg_square3D_matrix_determinant_sub_matrix (double *det, t8dg_square_3D_matrix_t matrix, int dim);

void                t8dg_square3D_matrix_scale (t8dg_square_3D_matrix_t matrix, double scaling_factor, int dim);

void                t8dg_square3D_matrix_apply_rotation_reflection_matrix_to_matrix (t8dg_square_3D_matrix_t
                                                                                     matrix_invers,
                                                                                     t8dg_square_3D_matrix_t matrix,
                                                                                     int idx_rotation_reflection, int dim);

#endif /* SRC_T8DG_DMATRIX_H_ */
