#include <gtest/gtest.h>
#include "../src/t8dg_quadrature.h"
#include "../src/t8dg_vertexset.h"
/* This fail just demonstrates that the test framework is working. */

static double
const_one (const double x[3])
{
  return 1;
}

static double
x_cubed (const double x[3])
{
  return x[0] * x[0] * x[0];
}

static double
x_to_5 (const double x[3])
{
  return x[0] * x[0] * x[0] * x[0] * x[0];
}

/* If this test fails, something is really wrong. */
TEST (quadrature1D, integrate_const_one)
{
  int                 num_LGL;
  t8dg_quadrature_t  *quadrature;
  t8dg_vertexset_t   *vertexset;
  for (num_LGL = 1; num_LGL <= 4; num_LGL++) {
    vertexset = t8dg_vertexset_new_1D_LGL (num_LGL);
    quadrature = t8dg_quadrature_new (vertexset);
    EXPECT_NEAR (t8dg_quadrature_integrate_reference_element (quadrature, const_one), 1, 1e-10);
    t8dg_quadrature_destroy (&quadrature);
    t8dg_vertexset_destroy (&vertexset);
  }
}

TEST (quadrature1D, integrate_max_order)
{
  int                 num_LGL;
  t8dg_quadrature_t  *quadrature;
  t8dg_vertexset_t   *vertexset;
  t8dg_scalar_function_3d_fn integrand_fn;
  for (num_LGL = 3; num_LGL <= 4; num_LGL++) {
    vertexset = t8dg_vertexset_new_1D_LGL (num_LGL);
    quadrature = t8dg_quadrature_new (vertexset);
    switch (num_LGL) {
    case (3):
      integrand_fn = x_cubed;
      break;
    case (4):
      integrand_fn = x_to_5;
      break;
    }
    EXPECT_NEAR (t8dg_quadrature_integrate_reference_element (quadrature, integrand_fn), 1.0 / (2 * num_LGL - 2), 1e-10);
    t8dg_quadrature_destroy (&quadrature);
    t8dg_vertexset_destroy (&vertexset);
  }
}
