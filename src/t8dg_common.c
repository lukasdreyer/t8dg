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

#include "t8dg.h"
#include <t8_cmesh.h>
#include <t8_vec.h>
#include <t8_cmesh_vtk.h>
#include "t8dg_common.h"
#include <sc_mpi.h>

t8dg_scalar_function_3d_time_fn
t8dg_common_initial_cond_fn (int initial_cond_arg)
{
  switch (initial_cond_arg) {
  case (0):
    return t8dg_scalar3d_constant_one;
  case (1):
    return t8dg_scalar1d_hat_function;
  case (2):
    return t8dg_scalar1d_step_function;
  case (3):
    return t8dg_scalar3d_cos_product;
  case (4):
    return t8dg_scalar3d_norm_function;
  case (5):
    return t8dg_scalar2d_hat_function;
  case (6):
    return t8dg_scalar2d_step_function;
  case (7):
    return t8dg_scalar2d_triangle_step_function;
  case (8):
    return t8dg_scalar3d_step_function;
  case (9):
    return t8dg_circle_ring_step_function;
  case (10):
    return t8dg_scalar2d_angle;
  case (11):
    return t8dg_cylinder_ring_sin_product_fn;
  case (12):
    return t8dg_cylinder_ring_step_function;
  case (13):
    return t8dg_smooth_indicator1Dfn;
  case (14):
    return t8dg_smooth_indicator2Dfn;
  case (15):
    return t8dg_smooth_indicator3Dfn;
  case (16):
    return t8dg_scalar3d_constant_zero;
  case (17):
    return t8dg_circle_ring_sin_product_fn;
  default:
    return NULL;
  }
}

t8dg_scalar_function_3d_time_fn
t8dg_common_analytic_solution_fn (int initial_cond_arg, double diffusion_coefficient)
{
  if (diffusion_coefficient > 0) {
    switch (initial_cond_arg) {
    case (0):
      return t8dg_scalar3d_constant_one;
    case (3):
      return t8dg_scalar3d_cos_product;
    default:
      return NULL;
    }
  }
  else {
    return t8dg_common_initial_cond_fn (initial_cond_arg);
  }
}

double
t8dg_scalar3d_constant_one (const double x[3], const double t, void *fn_data)
{
  return 1;
}

double
t8dg_scalar2d_angle (const double x[3], const double t, void *fn_data)
{
  return fabs (atan2 (x[1], x[0]));
}

double
t8dg_scalar3d_cos_product (const double x[3], const double t, void *fn_data)
{
  t8dg_scalar3d_cos_product_data_t *cos_data = fn_data;
  int                 dimension = cos_data->dim;
  double              diffusion_coefficient = cos_data->diffusion_coefficient;
  return exp (-diffusion_coefficient * dimension * 4 * M_PI * M_PI * t) * cos (2 * M_PI * x[0]) * cos (2 * M_PI * x[1]) * cos (2 * M_PI *
                                                                                                                               x[2]);
}

double
t8dg_scalar3d_constant_zero (const double x[3], const double t, void *fn_data)
{
  return 0;
}

double
t8dg_scalar1d_hat_function (const double x[3], const double t, void *fn_data)
{
  return 0.5 - (fabs (0.5 - x[0]));
}

double
t8dg_scalar2d_hat_function (const double x[3], const double t, void *fn_data)
{
  double              center[3] = { 0.5, 0.5, 0 };
  return sqrt (0.5) - t8_vec_dist (x, center);
}

double
t8dg_scalar3d_norm_function (const double x[3], const double t, void *fn_data)
{
  return t8_vec_norm (x);
}

double
t8dg_scalar1d_step_function (const double x[3], const double t, void *fn_data)
{
  double              center[3] = { 0.5, 0, 0 };
  return t8_vec_dist (x, center) < 0.2;
}

double
t8dg_scalar2d_step_function (const double x[3], const double t, void *fn_data)
{
  double              center[3] = { 0.5, 0.5, 0 };
  return t8_vec_dist (x, center) < 0.15;
}

double
t8dg_scalar2d_triangle_step_function (const double x[3], const double t, void *fn_data)
{
  return x[0] > 0.3 && x[1] > 0.3 && x[0] + x[1] < 0.9;
}

double
t8dg_scalar3d_step_function (const double x[3], const double t, void *fn_data)
{
  double              center[3] = { 0.5, 0.5, 0.5 };
  return t8_vec_dist (x, center) < 0.15;
}

double
t8dg_circle_ring_step_function (const double x[3], const double t, void *fn_data)
{
  double              center[3] = { 1.5, 0.0, 0.0 };
  double              dist = t8_vec_dist (x, center);
  if (dist < 0.1)
    return 1;
  if (dist > 0.2)
    return 0;
  dist = (dist - 0.1) * 10;     /* transform to [0,1] */
  return (cos (dist * M_PI) + 1) / 2;
}

double
t8dg_cylinder_ring_sin_product_fn (const double x[3], const double t, void *fn_data)
{
  double              angle = atan2 (x[1], x[0]);
  double              radius = sqrt (x[0] * x[0] + x[1] * x[1]);
  double              h = x[2];
  return sin (angle) * sin ((radius - 1.5) * 2 * M_PI) * sin ((h - 0.5) * 2 * M_PI);
}

double
t8dg_circle_ring_sin_product_fn (const double x[3], const double t, void *fn_data)
{
  double              angle = atan2 (x[1], x[0]);
  double              radius = sqrt (x[0] * x[0] + x[1] * x[1]);
  return sin (angle) * sin ((radius - 1.5) * 2 * M_PI);
}

static double
t8dg_smooth_h (const double x)
{
  if (x <= 0)
    return 0;
  return exp (-1 / x);
}

static double
t8dg_smooth_g (const double x)
{
  return t8dg_smooth_h (1 - x) / (t8dg_smooth_h (x) + t8dg_smooth_h (1 - x));
}

double
t8dg_cylinder_ring_step_function (const double x[3], const double t, void *fn_data)
{
  double              center[3] = { 1.5, 0.0, 0.5 };
  double              dist = t8_vec_dist (x, center);
  if (dist < 0.1)
    return 1;
  if (dist > 0.2)
    return 0;
  dist = (dist - 0.1) * 10;     /* transform to [0,1] */
  return t8dg_smooth_g (dist);
}

double
t8dg_cylinder_ring_source_fn (const double x[3], const double t, void *fn_data)
{
  return 30 * t8dg_cylinder_ring_step_function (x, t, fn_data);
}

double
t8dg_cos_indicator1Dfn (const double x[3], const double t, void *fn_data)
{
  double              radius = 0.2;
  double              smoothing_factor = 0.1;
  double              center[3] = { 0.5, 0.0, 0.0 };
  double              dist = t8_vec_dist (x, center);
  if (dist < radius)
    return 1;
  if (dist > (1 + smoothing_factor) * radius)
    return 0;
  dist = (dist - radius) / (radius * smoothing_factor); /* transform to [0,1] */
  return (cos (dist * M_PI) + 1) / 2;
}

double
t8dg_cos_indicator2Dfn (const double x[3], const double t, void *fn_data)
{
  double              radius = 0.2;
  double              smoothing_factor = 0.1;
  double              center[3] = { 0.5, 0.5, 0.0 };
  double              dist = t8_vec_dist (x, center);
  if (dist < radius)
    return 1;
  if (dist > (1 + smoothing_factor) * radius)
    return 0;
  dist = (dist - radius) / (radius * smoothing_factor); /* transform to [0,1] */
  return (cos (dist * M_PI) + 1) / 2;
}

double
t8dg_cos_indicator3Dfn (const double x[3], const double t, void *fn_data)
{
  double              radius = 0.2;
  double              smoothing_factor = 0.1;
  double              center[3] = { 0.5, 0.5, 0.5 };
  double              dist = t8_vec_dist (x, center);
  if (dist < radius)
    return 1;
  if (dist > (1 + smoothing_factor) * radius)
    return 0;
  dist = (dist - radius) / (radius * smoothing_factor); /* transform to [0,1] */
  return (cos (dist * M_PI) + 1) / 2;
}

double
t8dg_smooth_indicator1Dfn (const double x[3], const double t, void *fn_data)
{
  double              radius = 0.2;
  double              smoothing_factor = 0.5;
  double              center[3] = { 0.5, 0.0, 0.0 };
  double              dist = t8_vec_dist (x, center);
  if (dist < radius)
    return 1;
  if (dist > (1 + smoothing_factor) * radius)
    return 0;
  dist = (dist - radius) / (radius * smoothing_factor); /* transform to [0,1] */
  return t8dg_smooth_g (dist);
}

double
t8dg_smooth_indicator2Dfn (const double x[3], const double t, void *fn_data)
{
  double              radius = 0.2;
  double              smoothing_factor = 0.5;
  double              center[3] = { 0.5, 0.5, 0.0 };
  double              dist = t8_vec_dist (x, center);
  if (dist < radius)
    return 1;
  if (dist > (1 + smoothing_factor) * radius)
    return 0;
  dist = (dist - radius) / (radius * smoothing_factor); /* transform to [0,1] */
  return t8dg_smooth_g (dist);
}

double
t8dg_smooth_indicator3Dfn (const double x[3], const double t, void *fn_data)
{
  double              radius = 0.2;
  double              smoothing_factor = 0.5;
  double              center[3] = { 0.5, 0.5, 0.5 };
  double              dist = t8_vec_dist (x, center);
  if (dist < radius)
    return 1;
  if (dist > (1 + smoothing_factor) * radius)
    return 0;
  dist = (dist - radius) / (radius * smoothing_factor); /* transform to [0,1] */
  return t8dg_smooth_g (dist);
}
