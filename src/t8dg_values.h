/*
 * t8dg_values.h
 *
 *  Created on: Apr 17, 2020
 *      Author: lukas
 */

#ifndef SRC_T8DG_VALUES_H_
#define SRC_T8DG_VALUES_H_

#include "t8dg_local_values.h"
#include "t8dg_global_values.h"
#include "t8dg_geometry.h"
#include "t8dg_mortar.h"
#include "t8dg_dof.h"
T8DG_EXTERN_C_BEGIN ();

typedef struct t8dg_values t8dg_values_t;

t8dg_values_t      *t8dg_values_new_LGL_hypercube (int dim, int num_LGL_vertices, t8dg_coarse_geometry_t * geometry, t8_forest_t forest);

void                t8dg_values_destroy (t8dg_values_t ** p_values);

void                t8dg_values_ghost_exchange (t8dg_values_t * values);

void                t8dg_values_partition (t8dg_values_t * values, t8_forest_t forest_partition);

/*Always call together with replace_iterate inbetween*/
void                t8dg_values_allocate_adapt (t8dg_values_t * values, t8_forest_t forest_adapt);
void                t8dg_values_cleanup_adapt (t8dg_values_t * values);
/*call these only between allocate_adapt and cleanup_adapt*/

void                t8dg_values_set_element_adapt (t8dg_values_t * values, t8_locidx_t itree, t8_locidx_t ielement_new);

void                t8dg_values_copy_element_values (t8dg_values_t * values, t8_locidx_t idata_old, t8_locidx_t idata_new);

void                t8dg_values_transform_parent_dof_to_child_dof (t8dg_values_t * values, t8dg_dof_values_t * dof_values,
                                                                   t8dg_dof_values_t * dof_values_adapt, t8_locidx_t itree,
                                                                   t8_locidx_t ielem_parent_old, t8_locidx_t ielem_child_new, int ichild);

void                t8dg_values_transform_child_dof_to_parent_dof (t8dg_values_t * values, t8dg_dof_values_t * dof_values,
                                                                   t8dg_dof_values_t * dof_values_adapt, t8_locidx_t itree,
                                                                   t8_locidx_t num_children, t8_locidx_t ielem_first_child_old,
                                                                   t8_locidx_t ielem_parent_new);

/*Matrix applications*/

void                t8dg_values_apply_mass_matrix (t8dg_values_t * values, t8dg_dof_values_t * dof_src, t8dg_dof_values_t * dof_result);

void                t8dg_values_apply_inverse_mass_matrix (t8dg_values_t * values, t8dg_dof_values_t * dof_src,
                                                           t8dg_dof_values_t * dof_result);

void                t8dg_values_apply_component_stiffness_matrix_dof (t8dg_values_t * values, int icomp, t8dg_dof_values_t * flux_dof,
                                                                      t8dg_dof_values_t * dest_dof);

void                t8dg_values_apply_stiffness_matrix_linear_flux_fn3D
  (t8dg_values_t * values, t8dg_linear_flux3D_fn flux_fn, t8dg_flux_data_base *flux_data, double time, t8dg_dof_values_t * src_dof,
   t8dg_dof_values_t * dest_dof);

void                t8dg_values_apply_boundary_integrals
  (t8dg_values_t * values, t8dg_dof_values_t * src_dof, t8dg_dof_values_t * dest_dof,
   t8dg_linear_flux3D_fn linear_flux, t8dg_flux_data_base *flux_data, t8dg_numerical_linear_flux3D_fn numerical_flux,
   void *numerical_flux_data, double time);

void                t8dg_values_apply_component_boundary_integrals (t8dg_values_t * values, t8dg_dof_values_t * src_dof,
                                                                    t8dg_dof_values_t * dest_dof, int icomp,
                                                                    t8dg_numerical_flux1D_fn numerical_flux, void *numerical_flux_data,
                                                                    double time);

/* Getter */
t8dg_global_values_t *t8dg_values_get_global_values (t8dg_values_t * values, t8_locidx_t itree, t8_locidx_t ielement);
t8dg_global_values_t *t8dg_values_get_global_values_adapt (t8dg_values_t * values, t8_locidx_t itree, t8_locidx_t ielement);

t8dg_global_values_t **t8dg_values_get_global_values_array (t8dg_values_t * values);

t8_forest_t         t8dg_values_get_forest (t8dg_values_t * values);

double              t8dg_values_get_ghost_exchange_time (t8dg_values_t * values);

int                 t8dg_values_get_dim (t8dg_values_t * values);

double              t8dg_values_element_norm_l2_squared (t8dg_values_t * values, t8dg_element_dof_values_t * element_dof, t8_locidx_t itree,
                                                         t8_locidx_t ielement);
double              t8dg_values_element_area (t8dg_values_t * values, t8_locidx_t itree, t8_locidx_t ielement);

double              t8dg_values_norm_l2 (t8dg_values_t * values, t8dg_dof_values_t * dof_values, sc_MPI_Comm comm);

double              t8dg_values_norm_l2_rel (t8dg_values_t * values, t8dg_dof_values_t * dof_values,
                                             t8dg_scalar_function_3d_time_fn analytical_sol_fn, double time, void *analytical_sol_data,
                                             sc_MPI_Comm comm);
double              t8dg_values_norm_l_infty_rel (t8dg_values_t * values, t8dg_dof_values_t * dof_values,
                                                  t8dg_scalar_function_3d_time_fn analytical_sol_fn, double time, void *analytical_sol_data,
                                                  sc_MPI_Comm comm);

void                t8dg_values_interpolate_scalar_function_3d_time
  (t8dg_values_t * values, t8dg_scalar_function_3d_time_fn function, double time, void *function_data, t8dg_dof_values_t * dof_values);

#if T8_WITH_PETSC

void                t8dg_values_block_precon_apply_boundary_integrals (t8dg_values_t * values, t8dg_dof_values_t * src_dof,
                                                                       t8dg_dof_values_t * dest_dof, t8dg_linear_flux3D_fn linear_flux,
                                                                       t8dg_flux_data_base *flux_data, t8dg_numerical_linear_flux3D_fn numerical_flux,
                                                                       void *numerical_flux_data, double time, int selector);

void                t8dg_values_block_precon_apply_component_boundary_integrals (t8dg_values_t * values, t8dg_dof_values_t * src_dof,
                                                                                 t8dg_dof_values_t * dest_dof, int icomp,
                                                                                 t8dg_numerical_flux1D_fn numerical_flux,
                                                                                 void *numerical_flux_data, double time, int selector);

t8_locidx_t         t8dg_values_count_num_local_dofs (t8dg_values_t * values);

void                t8dg_values_mg_lvl_prepare_next_interpolation_step (t8dg_values_t * values, t8dg_mortar_array_t * mortar_array_lvl);

void
 
 
 
 
 
 
 
 t8dg_values_mg_lvl_set_interpolation_step (t8dg_values_t * values, t8_forest_t forest, t8_forest_t forest_adapt,
                                            t8dg_local_values_t * local_values, t8dg_local_values_t * local_values_adapt);

void
 
 
 
 
 
 
 
 t8dg_values_mg_lvl_allocate_properties (t8dg_values_t * values, int num_mg_lvls, t8_forest_t * forests,
                                         t8dg_local_values_t ** local_values_lvl, t8dg_mortar_array_t ** mortar_array_lvl);

t8dg_local_values_t **t8dg_values_get_local_values (t8dg_values_t * values);

t8dg_mortar_array_t **t8dg_values_get_mortar_array (t8dg_values_t * values);

t8dg_coarse_geometry_t *t8dg_values_get_coarse_geometry (t8dg_values_t * values);

#endif

T8DG_EXTERN_C_END ();

#endif /* SRC_T8DG_VALUES_H_ */
