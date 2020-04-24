#include <gtest/gtest.h>
#include "../src/t8dg.h"
#include "../src/t8dg_functionbasis.h"
#include "../src/t8dg_vertexset.h"
#include "../src/t8dg_dmatrix.h"
#include "../src/t8dg_sc_array.h"

TEST (functionbasis, interpolationMatrix_LGL2)
{
  t8dg_vertexset_t   *lgl_vertexset = t8dg_vertexset_new_1D_LGL (2);
  t8dg_vertexset_t   *lgl_left, *lgl_right;
  t8dg_functionbasis_t *functionbasis = t8dg_functionbasis_new_1D_Lagrange (lgl_vertexset);
  lgl_left = t8dg_vertexset_new_childvertexset_1D (lgl_vertexset, 0);
  lgl_right = t8dg_vertexset_new_childvertexset_1D (lgl_vertexset, 1);

  t8dg_dmatrix_t     *interpolation_matrix_left = t8dg_functionbasis_Lagrange_interpolation_matrix (functionbasis, lgl_left);
  t8dg_dmatrix_t     *interpolation_matrix_right = t8dg_functionbasis_Lagrange_interpolation_matrix (functionbasis, lgl_right);

  EXPECT_EQ (t8dg_dmatrix_at (interpolation_matrix_left, 0, 0), 1);
  EXPECT_EQ (t8dg_dmatrix_at (interpolation_matrix_left, 0, 1), 0);
  EXPECT_EQ (t8dg_dmatrix_at (interpolation_matrix_left, 1, 0), 0.5);
  EXPECT_EQ (t8dg_dmatrix_at (interpolation_matrix_left, 1, 1), 0.5);

  EXPECT_EQ (t8dg_dmatrix_at (interpolation_matrix_right, 0, 0), 0.5);
  EXPECT_EQ (t8dg_dmatrix_at (interpolation_matrix_right, 0, 1), 0.5);
  EXPECT_EQ (t8dg_dmatrix_at (interpolation_matrix_right, 1, 0), 0);
  EXPECT_EQ (t8dg_dmatrix_at (interpolation_matrix_right, 1, 1), 1);

  t8dg_dmatrix_destroy (&interpolation_matrix_left);
  t8dg_dmatrix_destroy (&interpolation_matrix_right);
  t8dg_vertexset_destroy (&lgl_left);
  t8dg_vertexset_destroy (&lgl_right);
  t8dg_functionbasis_destroy (&functionbasis);
  t8dg_vertexset_destroy (&lgl_vertexset);
}

TEST (functionbasis, derivativeMatrix_LGL3)
{
  t8dg_vertexset_t   *lgl_vertexset = t8dg_vertexset_new_1D_LGL (3);
  t8dg_functionbasis_t *functionbasis = t8dg_functionbasis_new_1D_Lagrange (lgl_vertexset);
  t8dg_dmatrix_t     *derivativeMatrix = t8dg_functionbasis_Lagrange_derivative_matrix (functionbasis);
  EXPECT_EQ (t8dg_dmatrix_at (derivativeMatrix, 0, 0), -3);
  EXPECT_EQ (t8dg_dmatrix_at (derivativeMatrix, 0, 1), 4);
  EXPECT_EQ (t8dg_dmatrix_at (derivativeMatrix, 0, 2), -1);

  EXPECT_EQ (t8dg_dmatrix_at (derivativeMatrix, 1, 0), -1);
  EXPECT_EQ (t8dg_dmatrix_at (derivativeMatrix, 1, 1), 0);
  EXPECT_EQ (t8dg_dmatrix_at (derivativeMatrix, 1, 2), 1);

  EXPECT_EQ (t8dg_dmatrix_at (derivativeMatrix, 2, 0), 1);
  EXPECT_EQ (t8dg_dmatrix_at (derivativeMatrix, 2, 1), -4);
  EXPECT_EQ (t8dg_dmatrix_at (derivativeMatrix, 2, 2), 3);
  t8dg_dmatrix_destroy (&derivativeMatrix);
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
  sc_array_t         *derivative_dof_values;
  sc_array_t         *dof_values;
  int                 icoord_gradient = 0;
  int                 icoord_derivative = 0;
  for (dim = 2; dim <= 3; dim++) {
    for (num_lgl = 2; num_lgl < MAX_LGL_NUMBER; num_lgl++) {
      lgl_vertexset = t8dg_vertexset_new_1D_LGL (num_lgl);
      functionbasis = t8dg_functionbasis_new_hypercube_lagrange (dim, lgl_vertexset);
      derivative_dof_values = sc_array_new_count (sizeof (double), t8dg_functionbasis_get_num_dof (functionbasis));
      dof_values = t8dg_sc_array_duplicate (derivative_dof_values);

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
  t8dg_functionbasis_t *functionbasis_tensor;
  sc_array_t         *derivative_dof_values;
  sc_array_t         *dof_values;
  int                 icoord_gradient = 0;
  int                 icoord_derivative = 0;

  lgl3_vertexset = t8dg_vertexset_new_1D_LGL (3);
  lgl5_vertexset = t8dg_vertexset_new_1D_LGL (5);
  lgl7_vertexset = t8dg_vertexset_new_1D_LGL (7);
  fb_tensor[0] = t8dg_functionbasis_new_1D_Lagrange (lgl3_vertexset);
  fb_tensor[1] = t8dg_functionbasis_new_1D_Lagrange (lgl5_vertexset);
  fb_tensor[2] = t8dg_functionbasis_new_1D_Lagrange (lgl7_vertexset);
  functionbasis_tensor = t8dg_functionbasis_new_tensor (3, fb_tensor);
  derivative_dof_values = sc_array_new_count (sizeof (double), t8dg_functionbasis_get_num_dof (functionbasis_tensor));
  dof_values = t8dg_sc_array_duplicate (derivative_dof_values);

  for (icoord_gradient = 0; icoord_gradient < dim; icoord_gradient++) {
    for (icoord_derivative = 0; icoord_derivative < dim; icoord_derivative++) {
      t8dg_functionbasis_interpolate_scalar_fn (functionbasis_tensor, coord_identity, &icoord_gradient, dof_values);
      t8dg_functionbasis_apply_derivative_matrix (functionbasis_tensor, icoord_derivative, dof_values, derivative_dof_values);
      for (idof = 0; idof < t8dg_functionbasis_get_num_dof (functionbasis_tensor); idof++) {
        EXPECT_NEAR (*(double *) sc_array_index_int (derivative_dof_values, idof), (icoord_gradient == icoord_derivative), 1e-10);
      }
    }
  }

  sc_array_destroy (derivative_dof_values);
  sc_array_destroy (dof_values);
  t8dg_functionbasis_destroy (&functionbasis_tensor);
  t8dg_functionbasis_destroy (&fb_tensor[0]);
  t8dg_functionbasis_destroy (&fb_tensor[1]);
  t8dg_functionbasis_destroy (&fb_tensor[2]);
  t8dg_vertexset_destroy (&lgl3_vertexset);
  t8dg_vertexset_destroy (&lgl5_vertexset);
  t8dg_vertexset_destroy (&lgl7_vertexset);
}
