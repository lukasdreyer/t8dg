/*
 * t8dg_numerical_flux.h
 *
 *  Created on: Mar 17, 2020
 *      Author: lukas
 */

#ifndef SRC_T8DG_NUMERICAL_FLUX_H_
#define SRC_T8DG_NUMERICAL_FLUX_H_

#include <sc_containers.h>
#include "t8dg.h"

typedef struct t8dg_mortar t8dg_mortar_t;
typedef double      (*t8dg_linear_numerical_flux_1D_fn) (const double u_minus, const double u_plus,
                                                         const double flow_vector[3], const double normal_vector[3]);
typedef double      (*t8dg_linear_flux_velocity_3D_time_fn) (const double flux_velocity[3], const double x_vec[3], double t);

struct t8dg_mortar
{
  int                 number_face_quadrature_points;
  t8dg_locidx_t       elem_idx_minus, elem_idx_plus;

  /*one value for each quadrature point */
  sc_array_t         *u_minus;
  sc_array_t         *u_plus;

  sc_array_t         *fluxes;

};

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
double              t8dg_upwind_flux_1D (const double u_minus, const double u_plus, const double flow_vector[3],
                                         const double normal_vector[3]);

#endif /* SRC_T8DG_NUMERICAL_FLUX_H_ */
