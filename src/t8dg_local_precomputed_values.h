/*
 * t8dg_local_precomputed_values.h
 *
 *  Created on: Apr 5, 2020
 *      Author: lukas
 */

#ifndef SRC_T8DG_LOCAL_PRECOMPUTED_VALUES_H_
#define SRC_T8DG_LOCAL_PRECOMPUTED_VALUES_H_

#include <t8.h>
#include <t8_forest.h>
#include <t8_element.h>

#include "t8dg_quadrature.h"
#include "t8dg_geometry.h"
#include "t8dg_flux.h"
#include "t8dg_global_precomputed_values.h"

typedef struct t8dg_local_precomputed_values t8dg_local_precomputed_values_t;

/*TODO: normal vector on functionbasis!!*/

void                t8dg_local_precomputed_values_copy_element_values
  (t8dg_local_precomputed_values_t * incoming_values, t8_locidx_t incoming_idata,
   t8dg_local_precomputed_values_t * outgoing_values, t8_locidx_t outgoing_idata);

void                t8dg_local_precomputed_values_set_element (t8dg_local_precomputed_values_t * values,
                                                               const t8dg_geometry_transformation_data_t * geometry_data,
                                                               const t8dg_global_precomputed_values_t * global_values);

t8dg_local_precomputed_values_t *t8dg_local_precomputed_values_new (t8_forest_t forest, t8dg_global_precomputed_values_t * global_values);

void                t8dg_local_precomputed_values_destroy (t8dg_local_precomputed_values_t ** pvalues);

double             *t8dg_local_precomputed_values_get_transformed_gradient_tangential_vector (const t8dg_local_precomputed_values_t *
                                                                                              values, const t8_locidx_t idata,
                                                                                              const int iquad, const int idim);

/*TODO: Change to idof*/
double             *t8dg_local_precomputed_values_get_face_normal_vector (const t8dg_local_precomputed_values_t * values,
                                                                          const t8_locidx_t idata, const int iface, const int iquad);

void                t8dg_local_precomputed_values_element_multiply_trafo_quad_weight
  (const t8dg_local_precomputed_values_t * local_values, t8_locidx_t idata, sc_array_t * element_quad, sc_array_t * result_element_quad);

void                t8dg_local_precomputed_values_element_divide_trafo_quad_weight
  (const t8dg_local_precomputed_values_t * local_values, t8_locidx_t idata, sc_array_t * element_quad, sc_array_t * result_element_quad);

void                t8dg_local_precomputed_values_face_multiply_trafo_quad_weight
  (const t8dg_local_precomputed_values_t * local_values,
   t8_locidx_t idata, const int iface, sc_array_t * face_quad, sc_array_t * result_face_quad);

void                t8dg_local_precomputed_values_face_divide_trafo_quad_weight
  (const t8dg_local_precomputed_values_t * local_values,
   t8_locidx_t idata, const int iface, sc_array_t * face_quad, sc_array_t * result_face_quad);

void                t8dg_local_precomputed_values_element_multiply_flux_value
  (const t8dg_local_precomputed_values_t * local_values, const t8dg_flux_t * flux,
   const t8dg_geometry_transformation_data_t * geometry_data,
   t8dg_quadrature_t * quadrature, double current_time, int idim, sc_array_t * element_quad_values, sc_array_t * element_flux_quad_values);

void                t8dg_local_precomputed_values_ghost_exchange (t8_forest_t forest, t8dg_local_precomputed_values_t * local_values);

void                t8dg_local_precomputed_values_partition (t8_forest_t forest_old, t8_forest_t forest_partition,
                                                             t8dg_local_precomputed_values_t * local_values_old,
                                                             t8dg_local_precomputed_values_t * local_values_partition);

#endif /* SRC_T8DG_LOCAL_PRECOMPUTED_VALUES_H_ */
