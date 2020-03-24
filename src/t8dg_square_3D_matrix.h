/*
 * t8dg_square_3D_matrix.h
 *
 *  Created on: Mar 21, 2020
 *      Author: lukas
 */

#ifndef SRC_T8DG_SQUARE_3D_MATRIX_H_
#define SRC_T8DG_SQUARE_3D_MATRIX_H_
#include "t8dg.h"

typedef double t8dg_square_3D_matrix_t[DIM3][DIM3];


void t8dg_square3D_matrix_invert_sub_matrix(t8dg_square_3D_matrix_t matrix_invers, t8dg_square_3D_matrix_t matrix, int dim);

void t8dg_square3D_matrix_determinant_sub_matrix(double *det, t8dg_square_3D_matrix_t matrix, int dim);

void t8dg_square3D_matrix_scale(t8dg_square_3D_matrix_t matrix, double scaling_factor, int dim);

void t8dg_square3D_matrix_apply_rotation_reflection_matrix_to_matrix(t8dg_square_3D_matrix_t matrix_invers, t8dg_square_3D_matrix_t matrix, int idx_rotation_reflection, int dim);



#endif /* SRC_T8DG_SQUARE_3D_MATRIX_H_ */
