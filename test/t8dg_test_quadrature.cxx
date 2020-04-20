#include <gtest/gtest.h>
#include "../src/t8dg.h"
#include "../src/t8dg_quadrature.h"
#include "../src/t8dg_vertexset.h"
/* This fail just demonstrates that the test framework is working. */

static double
const_one (const double x[3], void *data)
{
  return 1;
}

static double
x_pow (const double x[3], void *data)
{
  int                 power;
  power = *(int *) data;
  return pow (x[0], power);
}

static double
x_times_one_minus_y_pow (const double x[3], void *data)
{
  int                 power;
  power = *(int *) data;
  return pow (x[0] * (1 - x[1]), power);
}

/* If this test fails, something is really wrong. */
TEST (quadrature1D, integrate_const_one)
{
  int                 num_LGL;
  t8dg_quadrature_t  *quadrature;
  t8dg_vertexset_t   *vertexset;
  for (num_LGL = 1; num_LGL <= MAX_LGL_NUMBER; num_LGL++) {
    vertexset = t8dg_vertexset_new_1D_LGL (num_LGL);
    EXPECT_EQ (t8dg_vertexset_get_num_vertices (vertexset), num_LGL);
    quadrature = t8dg_quadrature_new_vertexset (vertexset);
    EXPECT_EQ (t8dg_quadrature_get_num_element_vertices (quadrature), num_LGL);
    EXPECT_NEAR (t8dg_quadrature_integrate_reference_element (quadrature, const_one, NULL), 1, 1e-10);
    t8dg_quadrature_destroy (&quadrature);
    t8dg_vertexset_destroy (&vertexset);
  }
}

TEST (quadrature1D, integrate_max_order)
{
  int                 num_LGL, power;
  t8dg_quadrature_t  *quadrature;
  t8dg_vertexset_t   *vertexset;
  for (num_LGL = 2; num_LGL <= MAX_LGL_NUMBER; num_LGL++) {
    vertexset = t8dg_vertexset_new_1D_LGL (num_LGL);
    quadrature = t8dg_quadrature_new_vertexset (vertexset);
    power = 2 * num_LGL - 3;
    printf ("%i\n", num_LGL);
    EXPECT_NEAR (t8dg_quadrature_integrate_reference_element (quadrature, x_pow, &power), 1.0 / (power + 1), 1e-10);
    t8dg_quadrature_destroy (&quadrature);
    t8dg_vertexset_destroy (&vertexset);
  }
}

TEST (quadrature2D, integrate_const_one)
{
  int                 num_LGL;
  t8dg_quadrature_t  *quadrature_line;
  t8dg_quadrature_t  *quadrature_square;
  t8dg_vertexset_t   *vertexset;
  for (num_LGL = 1; num_LGL <= MAX_LGL_NUMBER; num_LGL++) {
    vertexset = t8dg_vertexset_new_1D_LGL (num_LGL);
    EXPECT_EQ (t8dg_vertexset_get_num_vertices (vertexset), num_LGL);
    quadrature_line = t8dg_quadrature_new_vertexset (vertexset);
    quadrature_square = t8dg_quadrature_new_tensor_square (quadrature_line);
    EXPECT_EQ (t8dg_quadrature_get_num_element_vertices (quadrature_square), num_LGL * num_LGL);
    EXPECT_NEAR (t8dg_quadrature_integrate_reference_element (quadrature_square, const_one, NULL), 1, 1e-10);
    t8dg_quadrature_destroy (&quadrature_square);
    t8dg_quadrature_destroy (&quadrature_line);
    t8dg_vertexset_destroy (&vertexset);
  }
}

TEST (quadrature2D, integrate_max_order)
{
  int                 num_LGL;
  t8dg_quadrature_t  *quadrature_line;
  t8dg_quadrature_t  *quadrature_square;
  t8dg_vertexset_t   *vertexset;
  int                 power;
  double              result;
  for (num_LGL = 2; num_LGL <= MAX_LGL_NUMBER; num_LGL++) {
    vertexset = t8dg_vertexset_new_1D_LGL (num_LGL);
    EXPECT_EQ (t8dg_vertexset_get_num_vertices (vertexset), num_LGL);
    quadrature_line = t8dg_quadrature_new_vertexset (vertexset);
    quadrature_square = t8dg_quadrature_new_tensor_square (quadrature_line);
    EXPECT_EQ (t8dg_quadrature_get_num_element_vertices (quadrature_square), num_LGL * num_LGL);
    power = 2 * num_LGL - 3;
    result = 1. / ((power + 1) * (power + 1));
    EXPECT_NEAR (t8dg_quadrature_integrate_reference_element (quadrature_square, x_times_one_minus_y_pow, &power), result, 1e-10);
    t8dg_quadrature_destroy (&quadrature_square);
    t8dg_quadrature_destroy (&quadrature_line);
    t8dg_vertexset_destroy (&vertexset);
  }
}
