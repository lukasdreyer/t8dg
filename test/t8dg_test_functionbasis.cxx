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

#include <gtest/gtest.h>
#include "../src/t8dg.h"
#include "../src/t8dg_functionbasis.h"
#include "../src/t8dg_vertexset.h"
#include "../src/t8dg_dmatrix.h"
#include "../src/t8dg_sc_array.h"

TEST (functionbasis, creation)
{
  t8dg_vertexset_t   *lgl_vertexset = t8dg_vertexset_new_1D_LGL (2);
  t8dg_functionbasis_t *functionbasis = t8dg_functionbasis_new_1D_lagrange (lgl_vertexset, 0);
  EXPECT_EQ (t8dg_functionbasis_get_num_dof (functionbasis), 2);
  t8dg_functionbasis_destroy (&functionbasis);
  t8dg_vertexset_destroy (&lgl_vertexset);
}

TEST (functionbasis, face_functionbasis_1D)
{
  int                 num_lgl;
  double              vertex[3];
  t8dg_vertexset_t   *lgl_vertexset;
  t8dg_functionbasis_t *functionbasis;
  t8dg_functionbasis_t *left_face_functionbasis;
  t8dg_functionbasis_t *right_face_functionbasis;
  for (num_lgl = 2; num_lgl <= MAX_LGL_NUMBER; num_lgl++) {
    lgl_vertexset = t8dg_vertexset_new_1D_LGL (num_lgl);
    functionbasis = t8dg_functionbasis_new_1D_lagrange (lgl_vertexset, 1);
    left_face_functionbasis = t8dg_functionbasis_get_face_functionbasis (functionbasis, 0);
    right_face_functionbasis = t8dg_functionbasis_get_face_functionbasis (functionbasis, 1);
    EXPECT_EQ (t8dg_functionbasis_get_num_dof (functionbasis), num_lgl);
    EXPECT_EQ (t8dg_functionbasis_get_num_dof (left_face_functionbasis), 1);
    EXPECT_EQ (t8dg_functionbasis_get_num_dof (right_face_functionbasis), 1);
    t8dg_functionbasis_get_lagrange_vertex (left_face_functionbasis, 0, vertex);
    EXPECT_EQ (vertex[0], 0);
    EXPECT_EQ (vertex[1], 0);
    EXPECT_EQ (vertex[2], 0);
    t8dg_functionbasis_get_lagrange_vertex (right_face_functionbasis, 0, vertex);
    EXPECT_EQ (vertex[0], 1);
    EXPECT_EQ (vertex[1], 0);
    EXPECT_EQ (vertex[2], 0);
    t8dg_functionbasis_destroy (&functionbasis);
    t8dg_vertexset_destroy (&lgl_vertexset);
  }
}

TEST (functionbasis, face_functionbasis_2D)
{
  int                 num_lgl = 2, iface, idof;
  double              vertex[3];
  t8dg_vertexset_t   *lgl_vertexset;
  t8dg_functionbasis_t *functionbasis;
  t8dg_functionbasis_t *face_functionbasis;
  lgl_vertexset = t8dg_vertexset_new_1D_LGL (num_lgl);
  functionbasis = t8dg_functionbasis_new_hypercube_lagrange (2, lgl_vertexset, 1);
  EXPECT_EQ (t8dg_functionbasis_get_num_dof (functionbasis), pow (num_lgl, 2));
  for (iface = 0; iface < 4; iface++) {
    face_functionbasis = t8dg_functionbasis_get_face_functionbasis (functionbasis, iface);
    EXPECT_EQ (t8dg_functionbasis_get_num_dof (face_functionbasis), num_lgl);
    for (idof = 0; idof < 2; idof++) {
      t8dg_functionbasis_get_lagrange_vertex (face_functionbasis, idof, vertex);
      EXPECT_EQ (vertex[0], (iface == 1) || ((iface == 2 || iface == 3) && (idof == 1)));
      EXPECT_EQ (vertex[1], (iface == 3) || ((iface == 0 || iface == 1) && (idof == 1)));
      EXPECT_EQ (vertex[2], 0);
    }
  }
  t8dg_functionbasis_destroy (&functionbasis);
  t8dg_vertexset_destroy (&lgl_vertexset);
}

TEST (functionbasis, face_functionbasis_3D)
{
  int                 num_lgl = 2, iface, idof;
  double              vertex[3];
  t8dg_vertexset_t   *lgl_vertexset;
  t8dg_functionbasis_t *functionbasis;
  t8dg_functionbasis_t *face_functionbasis;
  lgl_vertexset = t8dg_vertexset_new_1D_LGL (num_lgl);
  functionbasis = t8dg_functionbasis_new_hypercube_lagrange (3, lgl_vertexset, 1);
  EXPECT_EQ (t8dg_functionbasis_get_num_dof (functionbasis), pow (num_lgl, 3));
  for (iface = 0; iface < 6; iface++) {
    face_functionbasis = t8dg_functionbasis_get_face_functionbasis (functionbasis, iface);
    EXPECT_EQ (t8dg_functionbasis_get_num_dof (face_functionbasis), pow (num_lgl, 2));
    for (idof = 0; idof < 4; idof++) {
      t8dg_functionbasis_get_lagrange_vertex (face_functionbasis, idof, vertex);
      EXPECT_EQ (vertex[0], (iface == 1) || ((iface >= 2 && iface <= 5) && (idof == 1 || idof == 3)));
      EXPECT_EQ (vertex[1], (iface == 3) || ((iface == 0 || iface == 1) && (idof == 1 || idof == 3))
                 || ((iface == 4 || iface == 5) && (idof == 2 || idof == 3)));
      EXPECT_EQ (vertex[2], (iface == 5) || ((iface >= 0 && iface <= 3) && (idof == 2 || idof == 3)));
    }
  }
  t8dg_functionbasis_destroy (&functionbasis);
  t8dg_vertexset_destroy (&lgl_vertexset);
}

TEST (functionbasis, interpolationMatrix_LGL2)
{
  t8dg_vertexset_t   *lgl_vertexset = t8dg_vertexset_new_1D_LGL (2);
  t8dg_functionbasis_t *functionbasis = t8dg_functionbasis_new_1D_lagrange (lgl_vertexset, 0);

  t8dg_dmatrix_t     *interpolation_matrix_left = t8dg_functionbasis_get_lagrange_child_interpolation_matrix (functionbasis, 0);
  t8dg_dmatrix_t     *interpolation_matrix_right = t8dg_functionbasis_get_lagrange_child_interpolation_matrix (functionbasis, 1);

  EXPECT_EQ (t8dg_dmatrix_at (interpolation_matrix_left, 0, 0), 1);
  EXPECT_EQ (t8dg_dmatrix_at (interpolation_matrix_left, 0, 1), 0);
  EXPECT_EQ (t8dg_dmatrix_at (interpolation_matrix_left, 1, 0), 0.5);
  EXPECT_EQ (t8dg_dmatrix_at (interpolation_matrix_left, 1, 1), 0.5);

  EXPECT_EQ (t8dg_dmatrix_at (interpolation_matrix_right, 0, 0), 0.5);
  EXPECT_EQ (t8dg_dmatrix_at (interpolation_matrix_right, 0, 1), 0.5);
  EXPECT_EQ (t8dg_dmatrix_at (interpolation_matrix_right, 1, 0), 0);
  EXPECT_EQ (t8dg_dmatrix_at (interpolation_matrix_right, 1, 1), 1);

  t8dg_functionbasis_destroy (&functionbasis);
  t8dg_vertexset_destroy (&lgl_vertexset);
}

TEST (functionbasis, derivativeMatrix_LGL3)
{
  t8dg_vertexset_t   *lgl_vertexset = t8dg_vertexset_new_1D_LGL (3);
  t8dg_functionbasis_t *functionbasis = t8dg_functionbasis_new_1D_lagrange (lgl_vertexset, 0);
  t8dg_dmatrix_t     *derivativeMatrix = t8dg_functionbasis_get_lagrange_derivative_matrix (functionbasis, 0);
  EXPECT_EQ (t8dg_dmatrix_at (derivativeMatrix, 0, 0), -3);
  EXPECT_EQ (t8dg_dmatrix_at (derivativeMatrix, 0, 1), 4);
  EXPECT_EQ (t8dg_dmatrix_at (derivativeMatrix, 0, 2), -1);

  EXPECT_EQ (t8dg_dmatrix_at (derivativeMatrix, 1, 0), -1);
  EXPECT_EQ (t8dg_dmatrix_at (derivativeMatrix, 1, 1), 0);
  EXPECT_EQ (t8dg_dmatrix_at (derivativeMatrix, 1, 2), 1);

  EXPECT_EQ (t8dg_dmatrix_at (derivativeMatrix, 2, 0), 1);
  EXPECT_EQ (t8dg_dmatrix_at (derivativeMatrix, 2, 1), -4);
  EXPECT_EQ (t8dg_dmatrix_at (derivativeMatrix, 2, 2), 3);
  t8dg_functionbasis_destroy (&functionbasis);
  t8dg_vertexset_destroy (&lgl_vertexset);
}

static double
coord_identity (const double vec[DIM3], void *icoord)
{
  return vec[*(int *) icoord];
}

TEST (functionbasis, derivativeMatrix_hypercube)
{
  int                 num_lgl;
  int                 dim;
  int                 idof;
  t8dg_vertexset_t   *lgl_vertexset;
  t8dg_functionbasis_t *functionbasis;
  t8dg_element_dof_values_t *derivative_dof_values;
  t8dg_element_dof_values_t *dof_values;
  int                 icoord_gradient = 0;
  int                 icoord_derivative = 0;
  for (dim = 2; dim <= 3; dim++) {
    for (num_lgl = 2; num_lgl <= MAX_LGL_NUMBER; num_lgl++) {
      lgl_vertexset = t8dg_vertexset_new_1D_LGL (num_lgl);
      functionbasis = t8dg_functionbasis_new_hypercube_lagrange (dim, lgl_vertexset, 0);
      EXPECT_EQ (t8dg_functionbasis_get_num_dof (functionbasis), pow (num_lgl, dim));
      derivative_dof_values = sc_array_new_count (sizeof (double), t8dg_functionbasis_get_num_dof (functionbasis));
      dof_values = t8dg_element_dof_values_duplicate (derivative_dof_values);

      for (icoord_gradient = 0; icoord_gradient < dim; icoord_gradient++) {
        for (icoord_derivative = 0; icoord_derivative < dim; icoord_derivative++) {
          t8dg_functionbasis_interpolate_scalar_fn (functionbasis, coord_identity, &icoord_gradient, dof_values);
          t8dg_functionbasis_apply_derivative_matrix (functionbasis, icoord_derivative, dof_values, derivative_dof_values);
          for (idof = 0; idof < t8dg_functionbasis_get_num_dof (functionbasis); idof++) {
            EXPECT_NEAR (*(double *) sc_array_index_int (derivative_dof_values, idof), (icoord_gradient == icoord_derivative), 1e-10);
          }
        }
      }

      sc_array_destroy (derivative_dof_values);
      sc_array_destroy (dof_values);
      t8dg_functionbasis_destroy (&functionbasis);
      t8dg_vertexset_destroy (&lgl_vertexset);
    }
  }
}

TEST (functionbasis, derivativeMatrix3d_different_num_dof)
{
  int                 dim = 3;
  int                 idof;
  t8dg_vertexset_t   *lgl3_vertexset;
  t8dg_vertexset_t   *lgl5_vertexset;
  t8dg_vertexset_t   *lgl7_vertexset;
  t8dg_functionbasis_t *fb_tensor[3];
  t8dg_functionbasis_t *functionbasis_tensor2D;
  t8dg_functionbasis_t *functionbasis_tensor3D;
  sc_array_t         *derivative_dof_values;
  sc_array_t         *dof_values;
  int                 icoord_gradient = 0;
  int                 icoord_derivative = 0;

  lgl3_vertexset = t8dg_vertexset_new_1D_LGL (3);
  lgl5_vertexset = t8dg_vertexset_new_1D_LGL (5);
  lgl7_vertexset = t8dg_vertexset_new_1D_LGL (7);
  fb_tensor[0] = t8dg_functionbasis_new_1D_lagrange (lgl3_vertexset, 0);
  fb_tensor[1] = t8dg_functionbasis_new_1D_lagrange (lgl5_vertexset, 0);
  fb_tensor[2] = t8dg_functionbasis_new_1D_lagrange (lgl7_vertexset, 0);
  functionbasis_tensor2D = t8dg_functionbasis_new_tensor (fb_tensor[0], fb_tensor[1], 0);
  functionbasis_tensor3D = t8dg_functionbasis_new_tensor (functionbasis_tensor2D, fb_tensor[2], 0);

  derivative_dof_values = sc_array_new_count (sizeof (double), t8dg_functionbasis_get_num_dof (functionbasis_tensor3D));
  dof_values = t8dg_element_dof_values_duplicate (derivative_dof_values);

  for (icoord_gradient = 0; icoord_gradient < dim; icoord_gradient++) {
    for (icoord_derivative = 0; icoord_derivative < dim; icoord_derivative++) {
      printf ("%i %i\n", icoord_gradient, icoord_derivative);
      t8dg_functionbasis_interpolate_scalar_fn (functionbasis_tensor3D, coord_identity, &icoord_gradient, dof_values);
      t8dg_functionbasis_apply_derivative_matrix (functionbasis_tensor3D, icoord_derivative, dof_values, derivative_dof_values);
      for (idof = 0; idof < t8dg_functionbasis_get_num_dof (functionbasis_tensor3D); idof++) {
        EXPECT_NEAR (*(double *) sc_array_index_int (derivative_dof_values, idof), (icoord_gradient == icoord_derivative), 1e-10);
      }
    }
  }

  sc_array_destroy (derivative_dof_values);
  sc_array_destroy (dof_values);
  t8dg_functionbasis_destroy (&functionbasis_tensor3D);
  t8dg_functionbasis_destroy (&functionbasis_tensor2D);
  t8dg_functionbasis_destroy (&fb_tensor[0]);
  t8dg_functionbasis_destroy (&fb_tensor[1]);
  t8dg_functionbasis_destroy (&fb_tensor[2]);
  t8dg_vertexset_destroy (&lgl3_vertexset);
  t8dg_vertexset_destroy (&lgl5_vertexset);
  t8dg_vertexset_destroy (&lgl7_vertexset);
}
