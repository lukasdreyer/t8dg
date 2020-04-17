#include <gtest/gtest.h>
#include "../src/t8dg.h"
#include "../src/t8dg_functionbasis.h"
#include "../src/t8dg_vertexset.h"
#include "../src/t8dg_dmatrix.h"

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
