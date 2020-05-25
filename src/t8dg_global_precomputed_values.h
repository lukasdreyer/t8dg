/*
 * t8dg_global_precomputed_values.h
 *
 *  Created on: Apr 2, 2020
 *      Author: lukas
 */
#ifndef SRC_T8DG_GLOBAL_PRECOMPUTED_VALUES_H_
#define SRC_T8DG_GLOBAL_PRECOMPUTED_VALUES_H_

#include "t8dg_quadrature.h"
#include "t8dg_functionbasis.h"

T8DG_EXTERN_C_BEGIN ();

typedef struct t8dg_global_precomputed_values t8dg_global_precomputed_values_t;

t8dg_global_precomputed_values_t *t8dg_global_precomputed_values_new_1D_LGL (const int number_of_LGL_vertices);

void                t8dg_global_precomputed_values_destroy (t8dg_global_precomputed_values_t ** pvalues);

void                t8dg_global_precomputed_values_transform_element_dof_to_face_quad (const t8dg_global_precomputed_values_t * values,
                                                                                       const int iface,
                                                                                       const sc_array_t * element_dof_array,
                                                                                       sc_array_t * face_quad_array);

void                t8dg_global_precomputed_values_transform_face_quad_to_element_dof (t8dg_global_precomputed_values_t * values,
                                                                                       const int iface,
                                                                                       const sc_array_t * face_quad_array,
                                                                                       sc_array_t * element_dof_array);

void                t8dg_global_precomputed_values_transform_element_dof_to_element_quad (const t8dg_global_precomputed_values_t * values,
                                                                                          const sc_array_t * element_dof_array,
                                                                                          sc_array_t * element_quad_array);

void                t8dg_global_precomputed_values_transform_element_quad_to_element_dof (const t8dg_global_precomputed_values_t * values,
                                                                                          const sc_array_t * element_quad_array,
                                                                                          sc_array_t * element_dof_array);

void                t8dg_global_precomputed_values_element_apply_derivative_matrix_transpose
  (const t8dg_global_precomputed_values_t * global_values, sc_array_t * derivative_dof_values, sc_array_t * dof_values);

void                t8dg_global_precomputed_values_transform_element_dof_to_child_dof (const t8dg_global_precomputed_values_t *
                                                                                       global_values, const sc_array_t * element_dof,
                                                                                       sc_array_t * child_dof, const int ichild);

void                t8dg_global_precomputed_values_transform_child_dof_to_parent_dof (const t8dg_global_precomputed_values_t *
                                                                                      global_values,
                                                                                      const sc_array_t * const child_dof[MAX_SUBELEMENTS],
                                                                                      sc_array_t * parent_dof, const int num_children);

int                 t8dg_global_precomputed_values_get_num_dof (const t8dg_global_precomputed_values_t * values);

int                 t8dg_global_precomputed_values_get_num_faces (const t8dg_global_precomputed_values_t * values);

t8dg_quad_idx_t     t8dg_global_precomputed_values_get_num_elem_quad (const t8dg_global_precomputed_values_t * values);

t8dg_functionbasis_t *t8dg_global_precomputed_values_get_functionbasis (const t8dg_global_precomputed_values_t * values);

t8dg_quadrature_t  *t8dg_global_precomputed_values_get_quadrature (const t8dg_global_precomputed_values_t * values);

int                 t8dg_global_precomputed_values_get_dim (const t8dg_global_precomputed_values_t * values);

T8DG_EXTERN_C_END ();

#endif /* SRC_T8DG_GLOBAL_PRECOMPUTED_VALUES_H_ */
