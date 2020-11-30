/*
 * t8dg_local_values.h
 *
 *  Created on: Apr 5, 2020
 *      Author: lukas
 */

#ifndef SRC_T8DG_LOCAL_VALUES_H_
#define SRC_T8DG_LOCAL_VALUES_H_

#include <t8.h>
#include <t8_forest.h>
#include <t8_element.h>

#include "t8dg_quadrature.h"
#include "t8dg_geometry.h"
#include "t8dg_flux.h"
#include "t8dg_global_values.h"
#include "t8dg_dof.h"
#include "t8dg_quad.h"

T8DG_EXTERN_C_BEGIN ();

typedef struct t8dg_local_values t8dg_local_values_t;

t8dg_local_values_t *t8dg_local_values_new (t8_forest_t forest, t8dg_global_values_t * global_values[T8_ECLASS_COUNT],
                                            t8dg_coarse_geometry_t * coarse_geometry);

void                t8dg_local_values_destroy (t8dg_local_values_t ** pvalues);

/*Setter*/
void
 
 
 
 
 
 
 
 t8dg_local_values_copy_element_values (t8dg_local_values_t * src_values, t8_locidx_t src_idata,
                                        t8dg_local_values_t * dest_values, t8_locidx_t dest_idata);

/*TODO!!!!
void                t8dg_local_values_set_element (t8dg_local_values_t * values,
                                                               const t8dg_geometry_transformation_data_t * geometry_data,
                                                               const t8dg_global_values_t * global_values);
*/

/* Operations on quad_values */

void                t8dg_local_values_element_multiply_trafo_quad_weight
  (const t8dg_local_values_t * local_values, t8_locidx_t itree, t8_locidx_t ielement, t8dg_element_quad_values_t * src_element_quad,
   t8dg_element_quad_values_t * dest_element_quad);

void                t8dg_local_values_element_divide_trafo_quad_weight
  (const t8dg_local_values_t * local_values, t8_locidx_t itree, t8_locidx_t ielement, t8dg_element_quad_values_t * src_element_quad,
   t8dg_element_quad_values_t * dest_element_quad);

void                t8dg_local_values_face_multiply_trafo_quad_weight
  (const t8dg_local_values_t * local_values,
   t8_locidx_t itree, t8_locidx_t ielement, const int iface, t8dg_face_quad_values_t * src_face_quad,
   t8dg_face_quad_values_t * dest_face_quad);

void                t8dg_local_values_face_divide_trafo_quad_weight
  (const t8dg_local_values_t * local_values,
   t8_locidx_t itree, t8_locidx_t ielement, const int iface, t8dg_face_quad_values_t * src_face_quad,
   t8dg_face_quad_values_t * dest_face_quad);

/*TODO
void                t8dg_local_values_element_multiply_flux_value
  (const t8dg_local_values_t * local_values, const t8dg_flux_t * flux,
   const t8dg_geometry_transformation_data_t * geometry_data,
   t8dg_quadrature_t * quadrature, double current_time, int idim, sc_array_t * element_quad_values, sc_array_t * element_flux_quad_values);

*/

/*TODO, wer ist zustandig */
void                t8dg_values_transform_face_child_dof_to_parent_dof
  (const t8dg_global_values_t * global_values,
   sc_array_t * child_face_dof[MAX_SUBFACES],
   sc_array_t * parent_face_dof, const int num_face_children,
   const t8dg_local_values_t * local_values, t8_locidx_t idata_child[MAX_SUBFACES], t8_locidx_t idata_parent, int iface_parent,
   t8dg_local_values_t * local_values_partition);

void                t8dg_local_values_transform_child_dof_to_parent_dof (const t8dg_global_values_t * global_values,
                                                                         sc_array_t * child_dof[MAX_SUBELEMENTS],
                                                                         sc_array_t * parent_dof, const int num_children,
                                                                         const t8dg_local_values_t * local_values_old,
                                                                         const t8dg_local_values_t * local_values_new,
                                                                         t8_locidx_t idata_first_child, t8_locidx_t idata_parent);

/*new:*/

void                t8dg_local_values_apply_element_mass_matrix (t8dg_local_values_t * values, t8_locidx_t itree, t8_locidx_t ielement,
                                                                 t8dg_element_dof_values_t * src_element_dof,
                                                                 t8dg_element_dof_values_t * dest_element_dof);

void                t8dg_local_values_apply_element_inverse_mass_matrix (t8dg_local_values_t * values, t8_locidx_t itree,
                                                                         t8_locidx_t ielement, t8dg_element_dof_values_t * src_element_dof,
                                                                         t8dg_element_dof_values_t * dest_element_dof);

void                t8dg_local_values_apply_face_mass_matrix (t8dg_local_values_t * values, t8_locidx_t itree, t8_locidx_t ielement,
                                                              int iface, t8dg_face_dof_values_t * src_face_dof,
                                                              t8dg_face_dof_values_t * dest_face_dof);

void                t8dg_local_values_apply_face_inverse_mass_matrix (t8dg_local_values_t * values, t8_locidx_t itree, t8_locidx_t ielement,
                                                                      int iface, t8dg_face_dof_values_t * src_face_dof,
                                                                      t8dg_face_dof_values_t * dest_face_dof);

/* data exchange and partition */

void                t8dg_local_values_ghost_exchange (t8dg_local_values_t * local_values);

void                t8dg_local_values_partition (t8dg_local_values_t * local_values_old, t8dg_local_values_t * local_values_partition);

/*Getter*/

double             *t8dg_local_values_get_transformed_gradient_tangential_vector (const t8dg_local_values_t *
                                                                                  values, const t8_locidx_t idata,
                                                                                  const int iquad, const int idim);

double             *t8dg_local_values_get_face_normal_vector (const t8dg_local_values_t * values,
                                                              const t8_locidx_t idata, const int iface, const int iquad);

t8dg_global_values_t *t8dg_local_values_get_global_values_from_itree (t8dg_local_values_t * local_values, t8_locidx_t itree);

int                 t8dg_local_values_get_max_num_faces (const t8dg_local_values_t * local_values);

T8DG_EXTERN_C_END ();

#endif /* SRC_T8DG_LOCAL_VALUES_H_ */
