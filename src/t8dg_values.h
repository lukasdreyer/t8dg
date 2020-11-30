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

void                t8dg_values_ghost_exchange (t8dg_values_t * values);

/*Always call together with replace_iterate inbetween*/

void                t8dg_values_allocate_adapt (t8dg_values_t * values, t8_forest_t forest_adapt);

void                t8dg_values_cleanup_adapt (t8dg_values_t * values, t8_forest_t forest_adapt);

void                t8dg_values_adapt (t8dg_values_t * values, t8_forest_t forest_adapt);

void                t8dg_values_partition (t8dg_values_t * values, t8_forest_t forest_partition);

void                t8dg_values_apply_mass_matrix (t8dg_values_t * values, t8dg_dof_values_t * dof_src, t8dg_dof_values_t * dof_result);

void                t8dg_values_apply_inverse_mass_matrix (t8dg_values_t * values, t8dg_dof_values_t * dof_src,
                                                           t8dg_dof_values_t * dof_result);

void
 
 
 
 
 
 
 
 t8dg_values_apply_stiffness_matrix_linear_flux_fn3D (t8dg_values_t * values, t8dg_linear_flux3D_fn flux_fn, void *flux_data,
                                                      t8dg_dof_values_t * src_dof, t8dg_dof_values_t * dest_dof);

void
 
 
 
 
 
 
 
 t8dg_values_apply_boundary_integrals (t8dg_values_t * values, t8dg_dof_values_t * src_dof, t8dg_dof_values_t * dest_dof,
                                       t8dg_linear_flux3D_fn linear_flux, void *flux_data, t8dg_numerical_linear_flux3D_fn numerical_flux,
                                       void *numerical_flux_data, double time);

t8dg_global_values_t *t8dg_values_get_global_values_from_itree (t8dg_values_t * values, t8_locidx_t itree);

t8dg_global_values_t **t8dg_values_get_global_values_array (t8dg_values_t * values);

void                t8dg_values_copy_element_values (t8dg_values_t * values, t8_locidx_t idata_old, t8_locidx_t idata_new);

//TODO: move implementation to local_values
double              t8dg_values_norm_infty_rel (t8dg_values_t * values, t8dg_dof_values_t * numerical_dof,
                                                t8dg_dof_values_t * analytical_sol);

double              t8dg_values_norm_l2_squared_rel (t8dg_values_t * values, t8dg_dof_values_t * numerical_dof,
                                                     t8dg_dof_values_t * analytical_sol);

/*split into itree ielement? */
double              t8dg_values_element_norm_l2_squared (t8dg_values_t * values, t8dg_element_dof_values_t * element_dof,
                                                         t8_locidx_t idata);

//Used 
double              t8dg_values_transform_reference_vertex_and_evaluate (const double reference_vertex[3], void *scalar_fn_data);

T8DG_EXTERN_C_END ();

#endif /* SRC_T8DG_VALUES_H_ */
