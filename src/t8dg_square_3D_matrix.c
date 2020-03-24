/*
 * t8dg_square_3D_matrix.c
 *
 *  Created on: Mar 21, 2020
 *      Author: lukas
 */

#include "t8dg.h"
#include "t8dg_square_3D_matrix.h"

void
t8dg_square3D_matrix_invert_sub_matrix (t8dg_square_3D_matrix_t matrix_invers, t8dg_square_3D_matrix_t matrix, int dim)
{
  T8_ASSERT (dim > 0 && dim <= DIM3);
  double              det;
  t8dg_square3D_matrix_determinant_sub_matrix (&det, matrix, dim);
  if (dim == 1) {
    matrix_invers[0][0] = 1. / det;
  }
  else if (dim == 2) {
    matrix_invers[0][0] = 1. / det * matrix[1][1];
    matrix_invers[1][1] = 1. / det * matrix[0][0];
    matrix_invers[0][1] = -1. / det * matrix[0][1];
    matrix_invers[1][0] = -1. / det * matrix[1][0];
  }
  else if (dim == 3) {
    SC_ABORT ("not yet implemented");
  }
}

void
t8dg_square3D_matrix_determinant_sub_matrix (double *det, t8dg_square_3D_matrix_t matrix, int dim)
{
  T8_ASSERT (dim > 0 && dim <= DIM3);
  if (dim == 1) {
    *det = matrix[0][0];
  }
  else if (dim == 2) {
    *det = matrix[0][0] * matrix[1][1] - matrix[0][1] * matrix[1][0];
  }
  else if (dim == 3) {
    SC_ABORT ("not yet implemented");
  }
}

void
t8dg_square3D_matrix_copy (t8dg_square_3D_matrix_t matrix_result, t8dg_square_3D_matrix_t matrix, int dim)
{
  int                 ixdim, iydim;
  for (ixdim = 0; ixdim < dim; ixdim++) {
    for (iydim = 0; iydim < dim; iydim++) {
      matrix_result[ixdim][iydim] = matrix[ixdim][iydim];
    }
  }
}

void
t8dg_square3D_matrix_scale (t8dg_square_3D_matrix_t matrix, double scaling_factor, int dim)
{
  int                 ixdim, iydim;
  for (ixdim = 0; ixdim < dim; ixdim++) {
    for (iydim = 0; iydim < dim; iydim++) {
      matrix[ixdim][iydim] *= scaling_factor;
    }
  }
}

void
t8dg_square3D_matrix_apply_rotation_reflection_matrix_to_matrix (t8dg_square_3D_matrix_t matrix_result,
                                                                 t8dg_square_3D_matrix_t matrix, int idx_rotation_reflection, int dim)
{
  if (idx_rotation_reflection == 0) {
    t8dg_square3D_matrix_copy (matrix_result, matrix, dim);
  }
}
