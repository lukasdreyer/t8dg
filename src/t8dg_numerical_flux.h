/*
 * t8dg_numerical_flux.h
 *
 *  Created on: Mar 17, 2020
 *      Author: lukas
 */
/** @file t8dg_numerical_flux.h */

#ifndef SRC_T8DG_NUMERICAL_FLUX_H_
#define SRC_T8DG_NUMERICAL_FLUX_H_

#include <sc_containers.h>
#include "t8dg.h"
#include <t8.h>
#include <t8_forest.h>

/** linear numerical flux function for the 1D case  */
typedef double      (*t8dg_linear_numerical_flux_1D_fn) (const double u_minus, const double u_plus,
                                                         const double flow_vector[3], const double normal_vector[3]);
/** time dependent flow velocity vector field*/
typedef double      (*t8dg_linear_flux_velocity_3D_time_fn) (const double flux_velocity[3], const double x_vec[3], double t);

/** struct used to save the calculated numerical fluxes at quadrature points*/
typedef struct t8dg_mortar
{
  int                 number_face_quadrature_points;    /**< The number of face quadrature points*/
  t8_locidx_t         elem_idata_minus;                   /**< Local index of the element corresponding to u_minus */
  t8_locidx_t         elem_idata_plus;                    /**< Local index of the element corresponding to u_plus */

  int                 iface_plus, iface_minus;

  /*one value for each quadrature point */
  sc_array_t         *u_minus;                          /**< value of u on elem_minus at face quadrature points */
  sc_array_t         *u_plus;                           /**< value of u on elem_plus at face quadrature points */

  sc_array_t         *fluxes;                           /**< value of (cu)*.n at face quadrature points */
  int                 valid;                            /**< indicates wether the fluxes are already newly calculated this timestep*/

} t8dg_mortar_t;                /*maybe change to opaque handle */

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

sc_array_t         *t8dg_mortar_get_flux (t8dg_mortar_t * mortar);

t8dg_mortar_t      *t8dg_mortar_new (t8_forest_t forest, t8_locidx_t itree, t8_locidx_t ielement, int iface);

void                t8dg_mortar_destroy (t8dg_mortar_t ** pmortar);

void                t8dg_mortar_get_idata_iface (t8dg_mortar_t * mortar, t8_locidx_t * pidata, int *piface, int side);

#endif /* SRC_T8DG_NUMERICAL_FLUX_H_ */
