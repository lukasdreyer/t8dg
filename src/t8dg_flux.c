/*
 * t8dg_numerical_flux.c
 *
 *  Created on: Mar 17, 2020
 *      Author: lukas
 */

#include "t8dg.h"
#include "t8dg_flux.h"

#include <t8_vec.h>

#if 0
struct t8dg_linear_flux3D
{
  t8dg_linear_flux3D_fn flux_velocity_fn;                         /**< Divergence free vector field */
  t8dg_numerical_linear_flux3D_fn numerical_flux_fn;
  void               *flux_data;
};

void
t8dg_flux_destroy (t8dg_flux_t ** pflux)
{
  t8dg_linear_flux3D_t *flux = *pflux;
  if (flux->flux_data != NULL) {
    T8DG_FREE (flux->flux_data);
    flux->flux_data = NULL;
  }
  flux->flux_velocity_fn = NULL;
  flux->numerical_flux_fn = NULL;
  T8DG_FREE (flux);
  *pflux = NULL;
}

void
t8dg_flux_calulate_flux (const t8dg_flux_t * flux, const double x_vec[3], double flux_vec[3], const double t)
{
  flux->flux_velocity_fn (x_vec, flux_vec, t, flux->flux_data);
}

double
t8dg_flux_calculate_numerical_flux_value (const t8dg_flux_t * flux, const double u_minus, const double u_plus, const double flow_vector[3],
                                          const double normal_vector[3])
{
  double              numerical_flux[3];
  flux->numerical_flux_fn (u_minus, u_plus, flow_vector, normal_vector, numerical_flux);
  return t8_vec_dot (numerical_flux, normal_vector);
}
#endif
