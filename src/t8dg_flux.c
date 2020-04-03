/*
 * t8dg_numerical_flux.c
 *
 *  Created on: Mar 17, 2020
 *      Author: lukas
 */

#include "t8dg.h"
#include "t8dg_flux.h"
#include <t8_vec.h>

struct t8dg_linear_flux
{
  t8dg_linear_flux_velocity_time_fn flux_velocity_fn;        /**< Divergence free vector field */
  void               *flux_data;
};

typedef struct t8dg_linear_flux_1D_linear_geometry_data
{
  double              tangential_vector[3];
  double              flow_velocity;
} t8dg_linear_flux_1D_linear_geometry_data_t;

static void
t8dg_linear_flux_1D_linear_geometry_flux_fn (const double x_vec[3], double flux_vec[3], const double t, const void *flux_data)
{
  /*In this case independent of x_vec and t */
  double             *tangent = ((t8dg_linear_flux_1D_linear_geometry_data_t *) flux_data)->tangential_vector;
  double              flow_velocity = ((t8dg_linear_flux_1D_linear_geometry_data_t *) flux_data)->flow_velocity;
  t8_vec_axb (tangent, flux_vec, flow_velocity, 0);
}

t8dg_linear_flux_t *
t8dg_linear_flux_new_1D_linear_geometry (const double tangential_vector[3], const double flow_velocity)
{
  t8dg_linear_flux_t *flux = T8DG_ALLOC (t8dg_linear_flux_t, 1);
  t8dg_linear_flux_1D_linear_geometry_data_t *data = T8_ALLOC (t8dg_linear_flux_1D_linear_geometry_data_t, 1);
  flux->flux_velocity_fn = t8dg_linear_flux_1D_linear_geometry_flux_fn;
  data->flow_velocity = flow_velocity;
  t8_vec_axb (tangential_vector, data->tangential_vector, 1. / t8_vec_norm (tangential_vector), 0);
  flux->flux_data = data;
  return flux;
}

void
t8dg_linear_flux_destroy (t8dg_linear_flux_t ** pflux)
{
  t8dg_linear_flux_t *flux = *pflux;
  if (flux->flux_data != NULL) {
    T8DG_FREE (flux->flux_data);
    flux->flux_data = NULL;
  }
  flux->flux_velocity_fn = NULL;
  T8DG_FREE (flux);
  *pflux = NULL;
}

double
t8dg_linear_numerical_flux_upwind_1D (const double u_minus, const double u_plus, const double flow_vector[3], const double normal_vector[3])
{
  /*TODO: ASSERT linear dependence of flow_vector and normal_vector */
  if (t8_vec_dot (flow_vector, normal_vector) > 0)
    return u_minus * t8_vec_dot (flow_vector, normal_vector);
  else
    return u_plus * t8_vec_dot (flow_vector, normal_vector);
}

void
t8dg_linear_flux_calulate_flux (t8dg_linear_flux_t * linear_flux, double x_vec[3], double flux_vec[3], double t)
{
  linear_flux->flux_velocity_fn (x_vec, flux_vec, t, linear_flux->flux_data);
}
