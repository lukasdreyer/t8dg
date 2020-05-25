/*
 * t8dg_numerical_flux.h
 *
 *  Created on: Mar 17, 2020
 *      Author: lukas
 */
/** @file t8dg_numerical_flux.h */

#ifndef SRC_T8DG_NUMERICAL_FLUX_H_
#define SRC_T8DG_NUMERICAL_FLUX_H_

#include <t8.h>
#include <t8_forest.h>

#include "t8dg.h"
#include "t8dg_local_precomputed_values.h"
#include "t8dg_quadrature.h"
#include "t8dg_coarse_geometry.h"

T8DG_EXTERN_C_BEGIN ();

typedef struct t8dg_linear_flux t8dg_linear_flux_t;

/** linear numerical flux function for the 1D case  */
typedef double      (*t8dg_linear_numerical_flux_fn) (const double u_minus, const double u_plus,
                                                      const double flow_vector[3], const double normal_vector[3]);
/** time dependent flow velocity vector field*/
typedef void        (*t8dg_linear_flux_velocity_time_fn) (const double x_vec[3], double flux_velocity[3], const double t,
                                                          const void *flux_data);

t8dg_linear_flux_t *t8dg_linear_flux_new_1D_linear_geometry (const double tangential_vector[3], const double flow_velocity);

void                t8dg_linear_flux_destroy (t8dg_linear_flux_t ** pflux);

void                t8dg_linear_flux_calulate_flux (t8dg_linear_flux_t * linear_flux, double x_vec[3], double flux_vec[3], double t);

/** For linear dependent flow- and normal vector, values of u at both adjacent elements,
 * calculates the 1D upwind flux.
 * The normal vector points from element_minus to element_plus.
 * All values should be evaluated at the same point, F_e(x_q)
 *
 * \param [in] u_minus            		value of u at quadpoint on element_minus
 * \param [in] u_plus            		value of u at quadpoint on element_plus
 * \param [in] flow_vector            		flow vector
 * \param [in] normal_vector            	normal vector
 *
 * \return                      		Flux value
 */

double              t8dg_linear_numerical_flux_upwind_1D (const double u_minus, const double u_plus, const double flow_vector[3],
                                                          const double normal_vector[3]);

void                t8dg_flux_element_multiply_flux_value (const t8dg_linear_flux_t * linear_flux, sc_array_t * element_quad_values,
                                                           double current_time, t8dg_local_precomputed_values_t * local_values,
                                                           t8_forest_t forest, t8_locidx_t itree, t8_locidx_t ielement,
                                                           t8dg_quadrature_t * quadrature, t8dg_coarse_geometry_t * coarse_geometry);

T8DG_EXTERN_C_END ();

#endif /* SRC_T8DG_NUMERICAL_FLUX_H_ */
