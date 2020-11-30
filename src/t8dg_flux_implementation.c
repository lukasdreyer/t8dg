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
t8dg_linear_numerical_flux1D_central (const double u_minus, const double u_plus,
                                      const double flow_constant, const double normal_component, const void *numerical_flux_data)
{
  double              average;
  average = (u_minus + u_plus) / 2;
  return flow_constant * average * normal_component;
}
