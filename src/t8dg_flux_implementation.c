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

#include "t8_vec.h"
#include "t8dg_flux.h"
#include "t8dg_flux_implementation.h"

void
t8dg_linear_flux3D_constant_flux_fn (double x_vec[3], double flux_vec[3], double t, void *flux_data)
{
  /*In this case independent of x_vec and t */
  double             *tangent = ((t8dg_linear_flux3D_constant_flux_data_t *) flux_data)->flow_direction;
  double              flow_velocity = ((t8dg_linear_flux3D_constant_flux_data_t *) flux_data)->flow_velocity;
  t8_vec_axb (tangent, flux_vec, flow_velocity, 0);
}

void
t8dg_rotating_flux_2D_fn (double x_vec[3], double flux_vec[3], double t, void *flux_data)
{
  flux_vec[0] = x_vec[1];
  flux_vec[1] = -x_vec[0];
  flux_vec[2] = 0;
}

void
t8dg_spiral_flux_3D_fn (double x_vec[3], double flux_vec[3], double t, void *flux_data)
{
  flux_vec[0] = 2 * M_PI * x_vec[1];
  flux_vec[1] = -2 * M_PI * x_vec[0];
  flux_vec[2] = 1;
}

void
t8dg_linear_flux_1D_from_3D_fn (double x_vec[3], double *flux_value, double t, void *flux_data)
{
  double              flux_vec[3];
  t8dg_linear_flux_1D_from_3D_data_t *flux_1D_data = (t8dg_linear_flux_1D_from_3D_data_t *) flux_data;
  flux_1D_data->flux_3D_fn (x_vec, flux_vec, t, flux_1D_data->flux_3D_data);
  *flux_value = flux_vec[flux_1D_data->component];
}

double
t8dg_linear_numerical_flux3D_lax_friedrich_fn (const double u_minus, const double u_plus, const double flow_vector[3],
                                               const double normal_vector[3], const void *numerical_flux_data)
{
  double              flow_velocity_bound;
  double              average, jump;

  flow_velocity_bound = *(double *) numerical_flux_data;

  average = (u_minus + u_plus) / 2;
  jump = u_minus - u_plus;

  return t8_vec_dot (flow_vector, normal_vector) * average + flow_velocity_bound / 2 * jump;
}

double
t8dg_numerical_flux1D_central (const double u_minus, const double u_plus, const double normal_component, const void *numerical_flux_data)
{
  double              average;
  average = (u_minus + u_plus) / 2;
  return average * normal_component;
}

double
t8dg_numerical_flux1D_left (const double u_minus, const double u_plus, const double normal_component, const void *numerical_flux_data)
{
  double              u_value = u_minus;
  int                 reverse_factor = 1;
  if (numerical_flux_data != NULL && *(int *) numerical_flux_data == 1) {
    reverse_factor = 1;
    u_value = u_plus;
  }
  return reverse_factor * u_value * normal_component;
}

double
t8dg_numerical_flux1D_right (const double u_minus, const double u_plus, const double normal_component, const void *numerical_flux_data)
{
  double              u_value = u_plus;
  int                 reverse_factor = 1;
  if (numerical_flux_data != NULL && *(int *) numerical_flux_data == 1) {
    reverse_factor = 1;
    u_value = u_minus;
  }
  return reverse_factor * u_value * normal_component;
}
