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
#include "t8dg_geometry.h"
#include "t8dg_mortar.h"

typedef struct t8dg_precomputed_values_fn_evaluation_data
{
  t8dg_geometry_transformation_data_t *geometry_data;
  t8dg_scalar_function_3d_time_fn function;
  double              time;
} t8dg_precomputed_values_fn_evaluation_data_t;

void                t8dg_precomputed_values_apply_element_boundary_integral (t8dg_global_precomputed_values_t * global_values,
                                                                             t8dg_local_precomputed_values_t * local_values,
                                                                             t8dg_mortar_array_t * mortar_array, int idata,
                                                                             sc_array_t * element_result_dof);

void                t8dg_precomputed_values_apply_face_inverse_mass_matrix (const t8dg_global_precomputed_values_t * global_values,
                                                                            const t8dg_local_precomputed_values_t * local_values,
                                                                            const t8dg_locidx_t idata, const int iface,
                                                                            sc_array_t * dof_values, sc_array_t * result_dof_values);

void                t8dg_precomputed_values_apply_face_mass_matrix (const t8dg_global_precomputed_values_t * global_values,
                                                                    const t8dg_local_precomputed_values_t * local_values,
                                                                    const t8dg_locidx_t idata, const int iface, sc_array_t * dof_values,
                                                                    sc_array_t * result_dof_values);

void                t8dg_precomputed_values_apply_element_inverse_mass_matrix
  (const t8dg_global_precomputed_values_t * global_values,
   const t8dg_local_precomputed_values_t * local_values,
   const t8dg_locidx_t idata, sc_array_t * dof_values, sc_array_t * result_dof_values);

void                t8dg_precomputed_values_apply_element_mass_matrix (const t8dg_global_precomputed_values_t * global_values,
                                                                       const t8dg_local_precomputed_values_t * local_values,
                                                                       const t8dg_locidx_t idata, sc_array_t * dof_values,
                                                                       sc_array_t * result_dof_values);

void                t8dg_precomputed_values_transform_child_dof_to_parent_dof (const t8dg_global_precomputed_values_t * global_values,
                                                                               sc_array_t * child_dof[MAX_SUBELEMENTS],
                                                                               sc_array_t * parent_dof, const int num_children,
                                                                               const t8dg_local_precomputed_values_t * local_values_old,
                                                                               const t8dg_local_precomputed_values_t * local_values_new,
                                                                               t8_locidx_t idata_first_child, t8_locidx_t idata_parent);

void                t8dg_precomputed_values_transform_face_child_dof_to_parent_dof
  (const t8dg_global_precomputed_values_t * global_values,
   sc_array_t * child_face_dof[MAX_SUBFACES],
   sc_array_t * parent_face_dof, const int num_face_children,
   const t8dg_local_precomputed_values_t * local_values, t8_locidx_t idata_child[MAX_SUBFACES], t8_locidx_t idata_parent, int iface_parent);

double              t8dg_precomputed_values_element_norm_infty (sc_array_t * element_dof_values);

double              t8dg_precomputed_values_element_norm_l2_squared
  (sc_array_t * element_dof_values, t8dg_global_precomputed_values_t * global_values,
   t8dg_local_precomputed_values_t * local_values, t8_locidx_t idata);

double
             t8dg_precomputed_values_transform_reference_vertex_and_evaluate (const double reference_vertex[3], void *scalar_fn_data);

#endif /* SRC_T8DG_PRECOMPUTED_VALUES_H_ */
