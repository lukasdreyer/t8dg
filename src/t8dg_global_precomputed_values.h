/*
 * t8dg_global_precomputed_values.h
 *
 *  Created on: Apr 2, 2020
 *      Author: lukas
 */

#ifndef SRC_T8DG_GLOBAL_PRECOMPUTED_VALUES_H_
#define SRC_T8DG_GLOBAL_PRECOMPUTED_VALUES_H_

typedef struct t8dg_global_precomputed_values t8dg_global_precomputed_values_t;

void
               t8dg_global_precomputed_values_new_1D_LGL (t8dg_global_precomputed_values_t ** pvalues, const int number_of_LGL_vertices);

void                t8dg_global_precomputed_values_destroy (const t8dg_global_precomputed_values_t ** pvalues);

void
 
 
 
 
 
 
 
 t8dg_global_precomputed_values_transform_element_dof_to_face_quad (const t8dg_global_precomputed_values_t * values,
                                                                    const int iface,
                                                                    const sc_array_t * element_dof_array, sc_array_t * face_quad_array);

void
 
 
 
 
 
 
 
 t8dg_global_precomputed_values_transform_face_quad_to_element_dof (t8dg_global_precomputed_values_t * values,
                                                                    const int iface,
                                                                    const sc_array_t * face_quad_array, sc_array_t * element_dof_array);

#endif /* SRC_T8DG_GLOBAL_PRECOMPUTED_VALUES_H_ */
