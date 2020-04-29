/*
 * t8dg_numerical_flux.c
 *
 *  Created on: Mar 17, 2020
 *      Author: lukas
 */

#include "t8dg.h"
#include "t8dg_flux.h"

#include <t8_vec.h>

typedef void        (*t8dg_numerical_flux_fn) (const double u_minus, const double u_plus,
                                               const double flow_vector[3], const double normal_vector[3], double numerical_flux[3]);

typedef void        (*t8dg_flux_fn) (const double x_vec[3], double flux_velocity[3], const double t, const void *flux_data);

struct t8dg_flux
{
  t8dg_flux_fn        flux_velocity_fn;                  /**< Divergence free vector field */
  t8dg_numerical_flux_fn numerical_flux_fn;
  void               *flux_data;
};

typedef struct t8dg_linear_flux_constant_flux_data
{
  double              flow_direction[3];
  double              flow_velocity;
} t8dg_linear_flux_constant_flux_data_t;

static void
t8dg_linear_flux_constant_flux_fn (const double x_vec[3], double flux_vec[3], const double t, const void *flux_data)
{
  /*In this case independent of x_vec and t */
  double             *tangent = ((t8dg_linear_flux_constant_flux_data_t *) flux_data)->flow_direction;
  double              flow_velocity = ((t8dg_linear_flux_constant_flux_data_t *) flux_data)->flow_velocity;
  t8_vec_axb (tangent, flux_vec, flow_velocity, 0);
}

static void
t8dg_linear_numerical_flux_upwind (const double u_minus, const double u_plus, const double flow_vector[3], const double normal_vector[3],
                                   double numerical_flux[3])
{
  double              flow_velocity;
  flow_velocity = t8_vec_norm (flow_vector);

  t8_vec_axb (flow_vector, numerical_flux, (u_minus + u_plus) / 2, 0);
  t8_vec_axpy (normal_vector, numerical_flux, flow_velocity * (u_minus - u_plus) / 2);
}

t8dg_flux_t        *
t8dg_flux_new_linear_constant_flux (const double flow_direction[3], const double flow_velocity)
{
  t8dg_flux_t        *flux = T8DG_ALLOC (t8dg_flux_t, 1);
  t8dg_linear_flux_constant_flux_data_t *data = T8_ALLOC (t8dg_linear_flux_constant_flux_data_t, 1);
  flux->flux_velocity_fn = t8dg_linear_flux_constant_flux_fn;
  data->flow_velocity = flow_velocity;
  t8_vec_axb (flow_direction, data->flow_direction, 1. / t8_vec_norm (flow_direction), 0);
  flux->flux_data = data;
  flux->numerical_flux_fn = t8dg_linear_numerical_flux_upwind;
  return flux;
}

void
t8dg_flux_destroy (t8dg_flux_t ** pflux)
{
  t8dg_flux_t        *flux = *pflux;
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
