#include "t8dg.h"
#include <t8_cmesh.h>
#include <t8_vec.h>
#include <t8_cmesh_vtk.h>
#include "t8dg_common.h"
#include <sc_mpi.h>

#include <../example/common/t8_example_common.h>

t8dg_scalar_function_3d_time_fn
t8dg_common_initial_cond_fn (int initial_cond_arg)
{
  switch (initial_cond_arg) {
  case (0):
    return t8_scalar3d_constant_one;
  case (1):
    return t8dg_scalar1d_hat_function;
  case (2):
    return t8_scalar3d_step_function;
  case (3):
    return t8_scalar3d_sinx;
  case (4):
    return t8dg_scalar3d_norm_function;
  case (5):
    return t8dg_scalar2d_hat_function;
  case (6):
    return t8dg_scalar2d_step_function;
  case (7):
    return t8dg_scalar2d_triangle_step_function;
  default:
    return t8_scalar3d_constant_zero;
  }
}

double
t8dg_scalar1d_hat_function (const double x[3], const double t)
{
  return 0.5 - (fabs (0.5 - x[0]));
}

double
t8dg_scalar2d_hat_function (const double x[3], const double t)
{
  double              center[3] = { 0.5, 0.5, 0 };
  return sqrt (0.5) - t8_vec_dist (x, center);
}

double
t8dg_scalar3d_norm_function (const double x[3], const double t)
{
  return t8_vec_norm (x);
}

double
t8dg_scalar2d_step_function (const double x[3], const double t)
{
  double              center[3] = { 0.5, 0.5, 0 };
  return t8_vec_dist (x, center) < 0.15;
}

double
t8dg_scalar2d_triangle_step_function (const double x[3], const double t)
{
  return x[0] > 0.3 && x[1] > 0.3 && x[0] + x[1] < 0.9;
}
