/*
 * t8dg_precomputed_values.h
 *
 *  Created on: Apr 17, 2020
 *      Author: lukas
 */

#ifndef SRC_T8DG_PRECOMPUTED_VALUES_H_
#define SRC_T8DG_PRECOMPUTED_VALUES_H_

#include "t8dg_local_precomputed_values.h"
#include "t8dg_global_precomputed_values.h"

void                t8dg_precomputed_values_apply_element_inverse_mass_matrix
  (const t8dg_global_precomputed_values_t * global_values,
   const t8dg_local_precomputed_values_t * local_values,
   const t8dg_locidx_t idata, const sc_array_t * dof_values, sc_array_t * result_dof_values);

void                t8dg_precomputed_values_apply_element_mass_matrix (const t8dg_global_precomputed_values_t * global_values,
                                                                       const t8dg_local_precomputed_values_t * local_values,
                                                                       const t8dg_locidx_t idata, const sc_array_t * dof_values,
                                                                       sc_array_t * result_dof_values);

void                t8dg_precomputed_values_transform_child_dof_to_parent_dof (const t8dg_global_precomputed_values_t * global_values,
                                                                               const sc_array_t * const child_dof[MAX_SUBELEMENTS],
                                                                               sc_array_t * parent_dof, const int num_children,
                                                                               t8dg_local_precomputed_values_t * local_values_old,
                                                                               t8dg_local_precomputed_values_t * local_values_new,
                                                                               t8_locidx_t idata_first_child, t8_locidx_t idata_parent);

double              t8dg_precomputed_values_element_norm_infty (sc_array_t * element_dof_values);

double
 
 
 
 
 
 
 t8dg_precomputed_values_element_norm_l2_squared (sc_array_t * element_dof_values, t8dg_global_precomputed_values_t * global_values,
                                                  t8dg_local_precomputed_values_t * local_values, t8_locidx_t idata);

#endif /* SRC_T8DG_PRECOMPUTED_VALUES_H_ */
