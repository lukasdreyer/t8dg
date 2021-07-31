/** @file t8dg_preconditioner.h */
/* This header hold the structs needed by the matrix preconditioner routines which provide the contexts of the several methods */

#ifndef SRC_T8DG_PRECONDITIONER_H_
#define SRC_T8DG_PRECONDITIONER_H_
#include <sc_containers.h>
#include "t8dg.h"
#include "t8dg_dof.h"
#include "t8dg_values.h"
#include "t8dg_advect_diff_problem.h"

#if T8_WITH_PETSC
#include <petscksp.h>
#endif

T8DG_EXTERN_C_BEGIN ();

#if T8_WITH_PETSC

typedef int         (*t8dg_corase_lvl_adapt_func_t) (t8_forest_t forest, t8_forest_t forest_from, t8_locidx_t itree, t8_locidx_t ielement,
                                                     t8_eclass_scheme_c * ts, int num_elements, t8_element_t * elements[]);

typedef struct
{
  t8dg_linear_advection_diffusion_problem_t *problem;
  t8_forest_t         forest_coarsened;
  t8dg_dof_values_t **dof_values;
  t8dg_values_t      *dg_values;
  t8dg_adapt_data_t  *adapt_data;
  size_t              num_local_dofs;
  PetscInt           *global_indexing;
  t8dg_mortar_array_t *tmp_mortar_coarse_lvl;
  t8_forest_t         problem_forest;
  t8dg_dof_values_t **dof_values_adapt;
} t8dg_mg_coarse_lvl_t;

typedef struct
{
  t8dg_mg_coarse_lvl_t *coarse_lvl;
  void               *user_data;
  t8dg_dof_values_t  *problem_dofs;
  t8dg_dof_values_t  *problem_dofs_derivation;
  t8dg_time_matrix_application time_derivative_func;
  double              current_point_in_time;
} t8dg_coarse_matrix_ctx_t;

typedef struct
{
  t8dg_linear_advection_diffusion_problem_t *problem;
  t8dg_mg_coarse_lvl_t *coarse_lvl;
  size_t              num_local_dofs_coarse_grid;
  size_t              num_local_dofs_fine_grid;
  PetscInt           *fine_lvl_global_indexing;
} t8dg_mg_interpolating_ctx_t;

void                t8dg_precon_write_dof_to_vec (t8dg_dof_values_t * dofs, Vec * p_vec, int *indexing);

void                t8dg_precon_write_vec_to_dof (Vec * p_vec, t8dg_dof_values_t * dofs);

void
 
 
 
 
 
 
 
 t8dg_mg_set_up_two_lvl_precon (t8dg_linear_advection_diffusion_problem_t * problem, t8dg_dof_values_t ** pdof_array,
                                t8dg_time_matrix_application time_derivative, t8dg_mg_coarse_lvl_t * coarse_lvl_mesh,
                                t8dg_corase_lvl_adapt_func_t coarsening_func, t8dg_mg_interpolating_ctx_t * res_prol_ctx,
                                t8dg_coarse_matrix_ctx_t * cmat_ctx, PetscInt * fine_forest_indexing);

void                t8dg_mg_swap_coarse_fine_data (t8dg_mg_coarse_lvl_t * coarse_lvl);

void
 
 
 
 
 
 
 
 t8dg_mg_construct_coarse_lvl (t8dg_linear_advection_diffusion_problem_t * problem, t8dg_dof_values_t ** pdof_array,
                               t8dg_mg_coarse_lvl_t * coarse_lvl_mesh, t8dg_corase_lvl_adapt_func_t coarsening_func);

void
               t8dg_mg_create_coarse_lvl_matrix (Mat * A_coarse, t8dg_mg_coarse_lvl_t * coarse_lvl, t8dg_coarse_matrix_ctx_t * cmat_ctx);

void                t8dg_mg_create_restriction_matrix (Mat * Restriction, t8dg_mg_interpolating_ctx_t * res_prol_ctx);

void                t8dg_mg_create_prolongation_matrix (Mat * Prolongation, t8dg_mg_interpolating_ctx_t * res_prol_ctx);

#endif

#if 0
/* neuer versuch, neues glueck */
/* Struct that provides a context for matrix-free preconditioning */
#if T8_WITH_PETSC
typedef struct
{
  t8dg_linear_advection_diffusion_problem_t *problem;
  t8_forest_t         forest_coarsened;
  t8_forest_t         tmp_dg_values_forest;
  t8_forest_t         tmp_forest;
  t8_forest_t         problem_forest;
  t8dg_values_t      *dg_values;
  t8dg_local_values_t tmp_local_values;
  t8dg_dof_values_t  *dof_values_fine_mesh;
  t8dg_dof_values_t  *dof_prolongated_coarse_grid_correction;
  t8dg_dof_values_t  *dof_current_coarse_grid_correction;
  t8dg_dof_values_t  *dofs_initial_problem;
  t8dg_mortar_array_t *tmp_mortar_array;
  t8dg_mortar_array_t *mortar_array_coarsened;
  //t8dg_dof_values_t *adapt_data_dof_values;
  t8dg_adapt_data_t  *adapt_data;
  PetscInt           *global_indexing_coarse_grid;
  double              current_time_point;
  size_t              num_local_dofs_coarse_grid;
  t8dg_time_matrix_application time_derivative_func;
  //t8dg_timestepping_data_t *time_data;
  //void               *user_data;

} t8dg_preconditioner_mg_ctx_t;
#endif

#if T8_WITH_PETSC
/* Creates a two level multigrid preconditioner with corresponding 'preconditioning context' */
PetscErrorCode      MGShellPCCreate (t8dg_preconditioner_mg_ctx_t ** pcshell);

/* Sets up the multigrid preconditioner - calculates a coarse mesh and other components */
PetscErrorCode      MGShellPCSetUp (PC pc, t8dg_advect_diff_problem_t * problem, t8dg_dof_values_t ** problem_dofs,
                                    t8dg_dof_values_t ** problem_dofs_adapt, t8dg_adapt_data_t ** adapt_data,
                                    t8dg_dg_values_t ** initial_dg_values, t8_forest_t * initial_problem_forest,
                                    t8dg_time_matrix_application time_derivative_func, double current_time_point);

/* Applies the two level multigrid preconditioner to a vector */
PetscErrorCode
MGShellPCApply (PC pc, Vec in, Vec out):
/* Frees the resources used by the multigrid preconditioning method */
PetscErrorCode
MGShellPCDestroy (PC pc);
#if 0
     void                t8dg_multigrid_coarse_level_solver (t8dg_linear_advection_diffusion_problem_t * problem, PC pc);
#endif
#endif

     void                t8dg_mulitgrid_preconditioner (t8dg_linear_advection_diffusion_problem_t * problem);

#endif
T8DG_EXTERN_C_END ();

#endif /* SRC_T8DG_PRECONDITIONER_H_ */
