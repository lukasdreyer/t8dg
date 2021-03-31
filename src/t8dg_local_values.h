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
void                t8dg_local_values_copy_element_values (t8dg_local_values_t * src_values, t8_locidx_t src_idata,
                                                           t8dg_local_values_t * dest_values, t8_locidx_t dest_idata);

void                t8dg_local_values_set_all_elements (t8dg_local_values_t * local_values);

void                t8dg_local_values_set_all_local_elements (t8dg_local_values_t * local_values);

void                t8dg_local_values_set_element (t8dg_local_values_t * local_values, t8_locidx_t itree, t8_locidx_t ielement);

void                t8dg_local_values_set_all_ghost_elements (t8dg_local_values_t * local_values);

void                t8dg_local_values_set_ghost_element (t8dg_local_values_t * local_values, t8_locidx_t ighosttree, t8_locidx_t ielement);

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

void
 
 
 
 
 
 
 
 t8dg_local_values_face_divide_trafo_quad_weight_idata_eclass (const t8dg_local_values_t * local_values,
                                                               t8_locidx_t idata, t8_eclass_t element_eclass, const int iface,
                                                               t8dg_face_quad_values_t * src_face_quad,
                                                               t8dg_face_quad_values_t * dest_face_quad);

void
 
 
 
 
 
 
 
 t8dg_local_values_face_multiply_trafo_quad_weight_idata_eclass (const t8dg_local_values_t * local_values,
                                                                 t8_locidx_t idata, t8_eclass_t element_eclass, const int iface,
                                                                 t8dg_face_quad_values_t * src_face_quad,
                                                                 t8dg_face_quad_values_t * dest_face_quad);

void                t8dg_local_values_element_multiply_directional_transformed_gradient_tangential_vector_component
  (const t8dg_local_values_t * local_values,
   t8_locidx_t itree, t8_locidx_t ielement,
   int derivative_direction, int icomp,
   t8dg_element_quad_values_t * src_element_quad_values, t8dg_element_quad_values_t * dest_element_quad_values);

void                t8dg_local_values_transform_child_dof_to_parent_dof (t8dg_local_values_t * local_values_old,
                                                                         t8dg_local_values_t * local_values_new,
                                                                         t8dg_element_dof_values_t * child_dof[MAX_SUBELEMENTS],
                                                                         t8dg_element_dof_values_t * parent_dof, int num_children,
                                                                         t8_locidx_t itree, t8_locidx_t ielem_first_child,
                                                                         t8_locidx_t ielem_parent);

void                t8dg_local_values_transform_orient_face_child_dof_to_parent_dof_hanging_nodes
  (t8dg_local_values_t * local_values,
   t8dg_face_dof_values_t * child_face_dof[MAX_SUBFACES],
   t8dg_face_dof_values_t * parent_face_dof,
   const int num_face_children,
   t8_locidx_t idata_child_neighbour[MAX_SUBFACES],
   t8_eclass_t element_eclass_children[MAX_SUBFACES],
   int iface_child_neighbour[MAX_SUBFACES], t8_locidx_t idata_parent, t8_eclass_t element_eclass_parent, int iface_parent, int orientation);

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

void                t8dg_local_values_apply_face_inverse_mass_matrix_idata_eclass
  (t8dg_local_values_t * values, t8_locidx_t idata, t8_eclass_t element_eclass, int iface,
   t8dg_face_dof_values_t * src_face_dof, t8dg_face_dof_values_t * dest_face_dof);

void                t8dg_local_values_apply_face_mass_matrix_idata_eclass
  (t8dg_local_values_t * values, t8_locidx_t idata, t8_eclass_t element_eclass, int iface,
   t8dg_face_dof_values_t * src_face_dof, t8dg_face_dof_values_t * dest_face_dof);

void                t8dg_local_values_apply_element_component_stiffness_matrix_dof
  (t8dg_local_values_t * local_values, t8_locidx_t itree,
   t8_locidx_t ielement, int icomp, t8dg_element_dof_values_t * src_element_dof, t8dg_element_dof_values_t * dest_element_dof);

/* data exchange and partition */

void                t8dg_local_values_ghost_exchange (t8dg_local_values_t * local_values);

void                t8dg_local_values_partition (t8dg_local_values_t * local_values_old, t8dg_local_values_t * local_values_partition);

/*Getter*/
int                 t8dg_local_values_get_max_num_faces (const t8dg_local_values_t * local_values);
t8dg_global_values_t *t8dg_local_values_get_global_values (t8dg_local_values_t * local_values, t8_locidx_t itree, t8_locidx_t ielement);
void
 
 
 
 
 
 
 
 t8dg_local_values_get_face_normal_vector (const t8dg_local_values_t * values,
                                           const t8_locidx_t itree, const t8_locidx_t ielement, const int iface, const int idof,
                                           double vector[3]);

void
 
 
 
 
 
 
 
 t8dg_local_values_get_face_normal_vector_idata_eclass (const t8dg_local_values_t * values,
                                                        const t8_locidx_t idata, const t8_eclass_t element_eclass, const int iface,
                                                        const int idof, double vector[3]);

t8dg_coarse_geometry_t *t8dg_local_values_get_coarse_geometry (t8dg_local_values_t * local_values);

double
 
 
 
 
 
 
 t8dg_local_values_element_norm_l2_squared (t8dg_local_values_t * local_values, t8dg_element_dof_values_t * element_dof_values,
                                            t8_locidx_t itree, t8_locidx_t ielement);

void                t8dg_local_values_element_error_ana_l2_squared (t8dg_local_values_t * local_values, t8dg_dof_values_t * dof_values,
                                                                    t8dg_dof_values_t * analytical_sol_dof, t8_locidx_t itree,
                                                                    t8_locidx_t ielement, double time, double *error_squared_summand,
                                                                    double *ana_norm_squared_summand);

T8DG_EXTERN_C_END ();

#endif /* SRC_T8DG_LOCAL_VALUES_H_ */
