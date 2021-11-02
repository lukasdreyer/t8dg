#include <t8_forest.h>
#include <t8_forest/t8_forest_iterate.h>
#include "t8dg_adapt.h"
#include "t8dg_preconditioner.h"

#if T8_WITH_PETSC

/* The maximum number of coarse levels to use within multigrid preconditioning is defined in 't8dg_preconditioner.h' */

/* The accuracy of convergence within the GMRES routines used during Block-Jacobi (/Block Gauss-Seidel) preconditioning */
#define T8DG_PRECON_BJ_GMRES_ACC 1.0e-7
/* The accuracy of convergence within the GMRES routines used during multigrid preconditioning */
#define T8DG_PRECON_MG_GMRES_ACC 1.0e-7
/* Whether the coarse level solver should be preconditioned with Block-Jacobi or not (true or false) */
#define T8DG_PRECON_MG_COARSE_GRID_PRECONDITIONER 0
/* Wheter a V-Cycle or W-Cycle should be used in the multigrid preconditioner (Options: PC_MG_CYCLE_V, PC_MG_CYCLE_W) */
#define T8DG_PRECON_MG_CYCLE_TYPE PC_MG_CYCLE_V
/* Number of Pre- and Post-Smoothing iterations on each multigrid level (except the finest) */
#define T8DG_PRECON_MG_NUM_SMOOTH_ITERATIONS 3
/* Whether the coarse level forest created during multigrid preconditioning should be saved within a vtk-file or not (true or false) */
#define T8DG_PRECON_MG_WRITE_OUT_COARSE_LVL 0
/* If the application time of the preconditioner should be measured (true or false) */
#define T8DG_PRECON_MEASURE_PC_TIME 1

/* Block-Preconditioner routines for matrix-free operations */
/* currently their name imply it is only block-jacobi which is implemented, but block-gauss-seidel is also possible, it uses the exact same structure and members; renaming these routines still has to be done */
extern PetscErrorCode JacobiShellPCCreate (t8dg_block_preconditioner_ctx_t **, t8dg_linear_advection_diffusion_problem_t *,
                                           t8dg_dof_values_t **, PetscInt *, int);
extern PetscErrorCode JacobiShellPCSetUp (PC);
extern PetscErrorCode JacobiShellPCApply (PC, Vec, Vec);
extern PetscErrorCode JacobiShellPCDestroy (PC);
extern PetscErrorCode MatMult_MF_Jacobi_Preconditioner (Mat, Vec, Vec);

/* Matrix-free routines used by a multigrid preconditioner which can handle an arbitrary number of coarse levels */
extern PetscErrorCode MatMult_MF_MG_LVL_Restriction (Mat, Vec, Vec);
extern PetscErrorCode MatMult_MF_MG_LVL_Prolongation (Mat, Vec, Vec);
extern PetscErrorCode MatMult_MF_MG_LVL_Coarse_Solver (Mat, Vec, Vec);
extern PetscErrorCode MatMult_MF_MG_LVL_Smoothing (Mat, Vec, Vec);

/* Matrix free p-Multigrid routines (not yet implmeneted) */
extern PetscErrorCode MatMult_MF_P_MG_LVL_Coarse_Solver (Mat, Vec, Vec);
extern PetscErrorCode MatMult_MF_P_MG_LVL_Restriction (Mat, Vec, Vec);
extern PetscErrorCode MatMult_MF_P_MG_LVL_Prolongation (Mat, Vec, Vec);
extern PetscErrorCode MatMult_MF_P_MG_LVL_Smoothing (Mat, Vec, Vec);

/* This function initializes a selected preconditioner */
void
t8dg_precon_initialize_preconditioner (PC * pc, int selector, t8dg_precon_general_preconditioner_t * general_precon, void *problem,
                                       t8dg_dof_values_t ** problem_dofs, t8dg_time_derivation_matrix_application_t time_derivative,
                                       Mat * A, PetscInt * vec_global_index)
{
#if T8DG_PRECON_MEASURE_PC_TIME
  general_precon->preconditioner_setup_time = -sc_MPI_Wtime ();
#endif

  /* Fill the general preconditioner with the selected preconditioning routine */
  switch (selector) {
  case 0:
    /* No preconditioning is selected */
    t8dg_precon_init_without_preconditioning (pc);
    break;
  case 1:
    /* Block-Jacobi-preconditioner is selected */
    t8dg_precon_init_jacobi (problem, problem_dofs, A, pc, &general_precon->jacobi_preconditioner_ctx, vec_global_index, selector);
    break;
  case 2:
    /* Currently, no preconditioner is assigned to '2' */
    t8dg_precon_init_without_preconditioning (pc);
    break;
  case 3:
    /* Block-Gauss-Seidel-preconditioner is selected */
    t8dg_precon_init_jacobi (problem, problem_dofs, A, pc, &general_precon->jacobi_preconditioner_ctx, vec_global_index, selector);
    break;
  case 4:
    /* MUltiple Level MG */
    general_precon->multiple_mg_lvls_ctx = T8DG_ALLOC_ZERO (t8dg_mg_levels_ctx_t, 1);
    t8dg_precon_init_multiple_level_mg (problem, problem_dofs, time_derivative, A, pc, vec_global_index,
                                        t8dg_timestepping_data_get_multigrid_levels (t8dg_advect_diff_problem_get_time_data
                                                                                     ((t8dg_linear_advection_diffusion_problem_t *)
                                                                                      problem)), &(general_precon->multiple_mg_lvls_ctx));
    break;
  default:
    /* Default equals no preconditioning */
    t8dg_precon_init_without_preconditioning (pc);
    break;
  }
#if T8DG_PRECON_MEASURE_PC_TIME
  general_precon->preconditioner_setup_time += sc_MPI_Wtime ();
#endif
}

/* This function desctroys a preconditioner and frees allocated memory which was used during the preconditioning routine */
void
t8dg_precon_destroy_preconditioner (PC * pc, int selector, t8dg_precon_general_preconditioner_t * general_precon)
{
  PetscErrorCode      ierr;

#if T8DG_PRECON_MEASURE_PC_TIME
  double              destroy_time = -sc_MPI_Wtime ();
#endif

  switch (selector) {
  case 0:
    /* No preconditioning was selected, therefore, nothing needs to be cleaned up */
    break;
  case 1:
    /* Block-Jacobi-preconditioner is selected  */
#if T8DG_PRECON_MEASURE_PC_TIME
    t8dg_global_essentialf ("The Preconditioner application took approx.:\n\t%f\n",
                            general_precon->jacobi_preconditioner_ctx->block_precon_application_time);
#endif
    t8dg_precon_destroy_block_preconditioner (general_precon->jacobi_preconditioner_ctx);
    ierr = JacobiShellPCDestroy (*pc);
    break;
  case 2:
    /* Currently, no preconditioner is assigned to '2' */
    break;
  case 3:
    /* Block-Gauss-Seidel-preconditioner is selected  */
#if T8DG_PRECON_MEASURE_PC_TIME
    t8dg_global_essentialf ("The Preconditioner application took approx.:\n\t%f\n",
                            general_precon->jacobi_preconditioner_ctx->block_precon_application_time);
#endif
    t8dg_precon_destroy_block_preconditioner (general_precon->jacobi_preconditioner_ctx);
    ierr = JacobiShellPCDestroy (*pc);
    break;
  case 4:
    /* Multiple Level MG */
#if T8DG_PRECON_MEASURE_PC_TIME
    t8dg_global_essentialf
      ("The Preconditioner application (wo smooth. on fine. lvl) took approx.:\n\tSmoothing (coarse lvls): %fs\n\tInterpolation: %fs\n\tCoarse LVL solve: %fs\n",
       general_precon->multiple_mg_lvls_ctx->mg_general_data->preconditioner_application_time_smoothing,
       general_precon->multiple_mg_lvls_ctx->mg_general_data->preconditioner_application_time_interpolation,
       general_precon->multiple_mg_lvls_ctx->mg_general_data->preconditioner_application_time_coarse_lvl);
#endif
    t8dg_precon_destroy_mg_levels (&(general_precon->multiple_mg_lvls_ctx));
    T8DG_FREE (general_precon->multiple_mg_lvls_ctx);
    break;
  default:
    /* No preconditioning was selected, therefore, nothing needs to be cleaned up */
    break;
  }

#if T8DG_PRECON_MEASURE_PC_TIME
  destroy_time += sc_MPI_Wtime ();
  //t8dg_global_essentialf ("Time to destroy the preconditioner took:\n\t%fs\n", destroy_time);
#endif
}

/* Initializes the without-preconditioning option */
PetscErrorCode
t8dg_precon_init_without_preconditioning (PC * pc)
{
  PetscErrorCode      ierr;
  ierr = PCSetType (*pc, PCNONE);
  CHKERRQ (ierr);

  return 0;
}

/**********************************************************************************************/
/************************ Beginning of Jacobi-Preconditioner-Routines *************************/
/**********************************************************************************************/
/****************************************** Take 2 ********************************************/

/* Initializes a Block-Jacobi preconditioner */
PetscErrorCode
t8dg_precon_init_jacobi (void *problem, t8dg_dof_values_t ** problem_dofs, Mat * A, PC * pc,
                         t8dg_block_preconditioner_ctx_t ** precon_jacobi_ctx, PetscInt * vec_global_index, int selector)
{
  PetscErrorCode      ierr;

  /* Create a PCSHELL which resembles a Block-Preconditioner */
  ierr = PCSetType (*pc, PCSHELL);
  CHKERRQ (ierr);

  /* Create an application context */
  JacobiShellPCCreate (precon_jacobi_ctx, (t8dg_linear_advection_diffusion_problem_t *) problem, problem_dofs, vec_global_index, selector);

  /* Assign the block preconditioner context to the preconditioner */
  ierr = PCShellSetContext (*pc, *precon_jacobi_ctx);
  CHKERRQ (ierr);

  /* Set the function that applies the preconditioner */
  ierr = PCShellSetApply (*pc, JacobiShellPCApply);
  CHKERRQ (ierr);

  /* Set the Destroy Function */
  ierr = PCShellSetDestroy (*pc, JacobiShellPCDestroy);
  CHKERRQ (ierr);

  /* Assign an appropriate name */
  if (selector == 1) {
    ierr = PCShellSetName (*pc, "Block-Jacobi-Preconditioner");
    CHKERRQ (ierr);
  }
  else if (selector == 3) {
    ierr = PCShellSetName (*pc, "Block-Gauss-Seidel-Preconditioner");
    CHKERRQ (ierr);
  }
  /* Set up the preconditioner */
  JacobiShellPCSetUp (*pc);

  return 0;
}

/* Creates a Jacobi preconditioner with corresponding 'preconditioning context' */
PetscErrorCode
JacobiShellPCCreate (t8dg_block_preconditioner_ctx_t ** precon_ctx, t8dg_linear_advection_diffusion_problem_t * problem,
                     t8dg_dof_values_t ** pdof_array, PetscInt * indexing, int selector)
{
  PetscErrorCode      ierr;
  t8dg_block_preconditioner_ctx_t *ctx;
  int                 iter;

  /* Allocate a new preconditioning context */
  ierr = PetscNew (&ctx);
  CHKERRQ (ierr);

  /* Assign values to the preconditiong context */
  (ctx->jacobi_ctx).current_point_in_time = t8dg_timestepping_data_get_current_time (t8dg_advect_diff_problem_get_time_data (problem));
  (ctx->jacobi_ctx).problem = (void *) problem;;
  (ctx->jacobi_ctx).problem_dofs = pdof_array;
  (ctx->jacobi_ctx).problem_dofs_derivation = t8dg_dof_values_duplicate (*pdof_array);
  /* This function is declared within 't8dg_advect_diff_problem_cxx' */
  (ctx->jacobi_ctx).jacobi_preconditioner_advect_diff_application = t8dg_advect_diff_problem_block_precon_time_derivative_variant;
  (ctx->jacobi_ctx).timestep = t8dg_timestepping_data_get_time_step (t8dg_advect_diff_problem_get_time_data (problem));
  (ctx->jacobi_ctx).current_implicit_coefficient = 1.0; /* in the case of the implicit euler method */

  /* process-local wise */
  (ctx->jacobi_ctx).num_local_dofs =
    (size_t) t8dg_dof_get_num_local_elements (*pdof_array) * t8dg_dof_get_max_num_element_dof (*pdof_array);

  /* Whether block-jacobi (== 1) or block-gauss-seidel (==3) has been choosen */
  (ctx->jacobi_ctx).selection = selector;

  /* Create a local indexing scheme */
  ierr = PetscMalloc1 ((ctx->jacobi_ctx).num_local_dofs, (&((ctx->jacobi_ctx).local_indexing)));
  CHKERRQ (ierr);
  for (iter = 0; iter < (ctx->jacobi_ctx).num_local_dofs; ++iter) {
    ((ctx->jacobi_ctx).local_indexing)[iter] = iter;
  }

  /* Global indexng scheme for the initial implicit system resulting from an implicit timestepping method */
  ctx->global_indexing = indexing;

  /* Assign the new preconditioning context */
  *precon_ctx = ctx;

  return 0;
}

/* Sets up a Block preconditioner - calculates the preconditioning matrix */
PetscErrorCode
JacobiShellPCSetUp (PC pc)
{
  PetscErrorCode      ierr;
  t8dg_block_preconditioner_ctx_t *ctx;

  ierr = PCShellGetContext (pc, (void **) &ctx);
  CHKERRQ (ierr);

  /* Declare the linear system which has to be solved in order to mimic the application of the inverse preconditioning-matrix */
  /* Create a local vector which will hold the result of the application of the inverse preconditioning matrix */
  ierr = VecCreateSeq (PETSC_COMM_SELF, (ctx->jacobi_ctx).num_local_dofs, &ctx->u);
  CHKERRQ (ierr);

  /* Create a local vector which holds the not yet preconditioned vector; it is going to be the right hand side of the linear system */
  ierr = VecCreateSeq (PETSC_COMM_SELF, (ctx->jacobi_ctx).num_local_dofs, &ctx->f);
  CHKERRQ (ierr);

  /* Create alocal matrix which applies the preconditioner to a vector */
  ierr =
    MatCreateShell (PETSC_COMM_SELF, (ctx->jacobi_ctx).num_local_dofs, (ctx->jacobi_ctx).num_local_dofs, (ctx->jacobi_ctx).num_local_dofs,
                    (ctx->jacobi_ctx).num_local_dofs, &ctx->jacobi_ctx, &ctx->M_jacobi);
  CHKERRQ (ierr);
  /* Set the matrix application of a multiplication in which this matrix takes place in */
  ierr = MatShellSetOperation (ctx->M_jacobi, MATOP_MULT, (void (*)(void)) MatMult_MF_Jacobi_Preconditioner);
  CHKERRQ (ierr);

  /* Create a krylov subspace method as solver regarding the linear system */
  ierr = KSPCreate (PETSC_COMM_SELF, &ctx->jacobi_preconditioner);
  CHKERRQ (ierr);
  ierr = KSPSetOperators (ctx->jacobi_preconditioner, ctx->M_jacobi, ctx->M_jacobi);
  CHKERRQ (ierr);
  ierr = KSPGetPC (ctx->jacobi_preconditioner, &ctx->pc_jacobi_preconditioner);
  CHKERRQ (ierr);
  /* Choose no preconditioer for this GMRES */
  t8dg_precon_init_without_preconditioning (&ctx->pc_jacobi_preconditioner);
  ierr = KSPSetType (ctx->jacobi_preconditioner, KSPGMRES);
  CHKERRQ (ierr);
  /* Set convergence tolerances of the preconditioning solver */
  ierr = KSPSetTolerances (ctx->jacobi_preconditioner, T8DG_PRECON_BJ_GMRES_ACC, T8DG_PRECON_BJ_GMRES_ACC, PETSC_DEFAULT, PETSC_DEFAULT);
  CHKERRQ (ierr);
  /* Set the solver up ready to use */
  ierr = KSPSetUp (ctx->jacobi_preconditioner);
  CHKERRQ (ierr);

  return 0;
}

/* Applies the Jacobi preconditioner to a vector */
PetscErrorCode
JacobiShellPCApply (PC pc, Vec in, Vec out)
{

  PetscErrorCode      ierr;
  t8dg_block_preconditioner_ctx_t *ctx;

  /* Get the preconditioning context */
  ierr = PCShellGetContext (pc, (void **) &ctx);
  CHKERRQ (ierr);

#if T8DG_PRECON_MEASURE_PC_TIME
  ctx->block_precon_application_time -= sc_MPI_Wtime ();
#endif

  /* Take the 'in' vector as the right hand side, so that the GMRES application mimics the application of the inverse of the preconditioning matrix */

  /* It is allowed to copy from a parallel to a sequential vec */
  /* Receive the 'in' vector as the right hand side */
  ierr = VecCopy (in, ctx->f);
  CHKERRQ (ierr);

  /* Solve the system in order to obtain the the application of the inverse preconditioning matrix; the result is stored in u */
  ierr = KSPSolve (ctx->jacobi_preconditioner, ctx->f, ctx->u);
  CHKERRQ (ierr);
#if 0
  /* View whether or not the GMRES did converge */
  ierr = KSPConvergedRateView (ctx->jacobi_preconditioner, PETSC_VIEWER_STDOUT_WORLD);
  CHKERRQ (ierr);
#endif
  /* Copy the result to the 'out' vector */
  ierr = VecCopy (ctx->u, out);
  CHKERRQ (ierr);

#if T8DG_PRECON_MEASURE_PC_TIME
  ctx->block_precon_application_time += sc_MPI_Wtime ();
#endif

  return 0;
}

/* Frees the resources used by the preconditioning methods */
PetscErrorCode
JacobiShellPCDestroy (PC pc)
{
  PetscErrorCode      ierr;
  t8dg_block_preconditioner_ctx_t *ctx;

  /* Get the preconditioning context */
  ierr = PCShellGetContext (pc, (void **) &ctx);
  CHKERRQ (ierr);

  /* Free the preconditioning context */
  ierr = PetscFree (ctx);
  CHKERRQ (ierr);

  return 0;
}

/* Destroys/frees the allocated memory of the block preconditioner */
PetscErrorCode
t8dg_precon_destroy_block_preconditioner (t8dg_block_preconditioner_ctx_t * ctx)
{
  PetscErrorCode      ierr;

  /* Free dof values */
  t8dg_dof_values_destroy (&((ctx->jacobi_ctx).problem_dofs_derivation));

  /* Free the allocated space of PETSc objects */
  ierr = PetscFree ((ctx->jacobi_ctx).local_indexing);
  CHKERRQ (ierr);
  ierr = MatDestroy (&ctx->M_jacobi);
  CHKERRQ (ierr);
  ierr = VecDestroy (&ctx->f);
  CHKERRQ (ierr);
  ierr = VecDestroy (&ctx->u);
  CHKERRQ (ierr);
  ierr = KSPDestroy (&ctx->jacobi_preconditioner);
  CHKERRQ (ierr);

  return 0;
}

/* Matrix-free application which resembles the application of the function describing the block jacobi preconditioning process */
PetscErrorCode
MatMult_MF_Jacobi_Preconditioner (Mat A, Vec in, Vec out)
{

  PetscErrorCode      ierr;
  t8dg_precon_block_matrix_ctx_t *jacobi_ctx;

  ierr = MatShellGetContext (A, &jacobi_ctx);
  CHKERRQ (ierr);

  /* Copy the petsc vec entries into the dof_values of the problem in order to compute their derivation */
  t8dg_precon_write_vec_to_dof (&in, *(jacobi_ctx->problem_dofs), jacobi_ctx->num_local_dofs);

  /* Calculate the time derivation of the degrees of freedom in the jacobi-preconditioning variant */
  (jacobi_ctx->jacobi_preconditioner_advect_diff_application) (*(jacobi_ctx->problem_dofs), jacobi_ctx->problem_dofs_derivation,
                                                               jacobi_ctx->current_point_in_time, jacobi_ctx->problem,
                                                               jacobi_ctx->selection);

  /* calculate the sum of the initially passed degrees of freedom ('in' Vector) and their (scaled) derivation */
  t8dg_dof_values_axpy (-(jacobi_ctx->timestep * jacobi_ctx->current_implicit_coefficient), jacobi_ctx->problem_dofs_derivation,
                        *(jacobi_ctx->problem_dofs));

  /* Copy the dof_values back to a petsc_vector */
  t8dg_precon_write_dof_to_vec (*(jacobi_ctx->problem_dofs), &out, jacobi_ctx->local_indexing, jacobi_ctx->num_local_dofs);

  return 0;
}

/**********************************************************************************************/
/*************************** End of Jacobi-Preconditioner-Routines ****************************/
/**********************************************************************************************/

/**********************************************************************************************/
/****************************** Beginning of Multigrid-Routines *******************************/
/**********************************************************************************************/

/******************/
/******************/
/* p MG */
/******************/
/******************/
/* Not yet implemented */
#if 0
void
t8dg_precon_init_p_mg (void *problem, t8dg_dof_values_t ** problem_dofs, t8dg_time_derivation_matrix_application_t time_derivative,
                       Mat * A_fine, PC * mg_pc, PetscInt * vec_global_index, int coarse_lvl_order, t8dg_p_mg_lvl_ctx_t ** mg_ctx)
{
  t8dg_p_mg_set_up_coarse_lvl ((t8dg_linear_advection_diffusion_problem_t *) problem, problem_dofs, coarse_lvl_order, *mg_ctx);

  t8dg_p_mg_set_up_restriction_operator (*mg_ctx);

  t8dg_p_mg_set_up_prolongation_operator (*mg_ctx);

  t8dg_p_mg_set_up_coarse_level (*mg_ctx);

  t8dg_p_mg_set_up_smoothing_operator (*mg_ctx, A_fine);

}

void
t8dg_p_mg_set_up_smoothing_operator (t8dg_p_mg_lvl_ctx_t * mg_ctx, Mat * A_fine)
{
  mg_ctx->Smoothing_Mat = *A_fine;
}

PetscErrorCode
t8dg_p_mg_set_up_coarse_level (t8dg_p_mg_lvl_ctx_t * mg_ctx)
{
  PetscErrorCode      ierr;
  /* Create the shell-matrix representing the Restriction operator */
  ierr =
    MatCreateShell (PETSC_COMM_WORLD, mg_ctx->p_mg_mat_ctx->num_coarse_dofs, mg_ctx->p_mg_mat_ctx->num_coarse_dofs, PETSC_DETERMINE,
                    PETSC_DETERMINE, mg_ctx->p_mg_mat_ctx, &(mg_ctx->Coarse_Mat));
  CHKERRQ (ierr);

  /* Define the (multiplicative) MatVec-Operation which resembles the application of the restriction matrix */
  ierr = MatShellSetOperation (mg_ctx->Coarse_Mat, MATOP_MULT, (void (*)(void)) MatMult_MF_P_MG_LVL_Coarse_Solver);
  CHKERRQ (ierr);
}

PetscErrorCode
t8dg_p_mg_set_up_prolongation_operator (t8dg_p_mg_lvl_ctx_t * mg_ctx)
{
  PetscErrorCode      ierr;
  /* Create the shell-matrix representing the Restriction operator */
  ierr =
    MatCreateShell (PETSC_COMM_WORLD, mg_ctx->p_mg_mat_ctx->num_fine_dofs, mg_ctx->p_mg_mat_ctx->num_coarse_dofs, PETSC_DETERMINE,
                    PETSC_DETERMINE, mg_ctx->p_mg_mat_ctx, &(mg_ctx->Prolongation));
  CHKERRQ (ierr);

  /* Define the (multiplicative) MatVec-Operation which resembles the application of the restriction matrix */
  ierr = MatShellSetOperation (mg_ctx->Prolongation, MATOP_MULT, (void (*)(void)) MatMult_MF_P_MG_LVL_Prolongation);
  CHKERRQ (ierr);
}

/* Creates the matrix-free restriction routines, interpolating from a fine level onto a coarse level */
PetscErrorCode
t8dg_p_mg_set_up_restriction_operator (t8dg_p_mg_lvl_ctx_t * mg_ctx)
{
  PetscErrorCode      ierr;
  /* Create the shell-matrix representing the Restriction operator */
  ierr =
    MatCreateShell (PETSC_COMM_WORLD, mg_ctx->p_mg_mat_ctx->num_coarse_dofs, mg_ctx->p_mg_mat_ctx->num_fine_dofs, PETSC_DETERMINE,
                    PETSC_DETERMINE, mg_ctx->p_mg_mat_ctx, &(mg_ctx->Restriction));
  CHKERRQ (ierr);

  /* Define the (multiplicative) MatVec-Operation which resembles the application of the restriction matrix */
  ierr = MatShellSetOperation (mg_ctx->Restriction, MATOP_MULT, (void (*)(void)) MatMult_MF_P_MG_LVL_Restriction);
  CHKERRQ (ierr);
}

PetscErrorCode
MatMult_MF_P_MG_LVL_Coarse_Solver (Mat C_Mat, Vec in, Vec out)
{
  PetscErrorCode      ierr;
  t8dg_p_mg_mat_ctx_t *c_appctx;

  /* Get the matrix-free application context */
  ierr = MatShellGetContext (C_Mat, &c_appctx);
  CHKERRQ (ierr);

  return 0;
}

PetscErrorCode
MatMult_MF_P_MG_LVL_Prolongation (Mat Prolongation, Vec in, Vec out)
{
  PetscErrorCode      ierr;
  t8dg_p_mg_mat_ctx_t *prol_ctx;

  /* Get the matrix-free application context */
  ierr = MatShellGetContext (Prolongation, &prol_ctx);
  CHKERRQ (ierr);

  return 0;
}

PetscErrorCode
MatMult_MF_P_MG_LVL_Restriction (Mat Restriction, Vec in, Vec out)
{
  PetscErrorCode      ierr;
  t8dg_p_mg_mat_ctx_t *res_ctx;

  /* Get the matrix-free application context */
  ierr = MatShellGetContext (Restriction, &res_ctx);
  CHKERRQ (ierr);

  return 0;
}

void
t8dg_p_mg_set_up_coarse_lvl (t8dg_linear_advection_diffusion_problem_t * problem, t8dg_dof_values_t ** problem_dofs, int coarse_lvl_order,
                             t8dg_p_mg_lvl_ctx_t * mg_ctx)
{
  /* Allocate a matrix context for restriction, prolongation and coarse level operators */
  mg_ctx->p_mg_mat_ctx = T8DG_ALLOC_ZERO (t8dg_p_mg_mat_ctx_t, 1);

  /* Initialize some members of the p_multigrid context */
  mg_ctx->p_mg_mat_ctx->mg_coarse_order = coarse_lvl_order;
  mg_ctx->p_mg_mat_ctx->problem_dofs = problem_dofs;
  mg_ctx->p_mg_mat_ctx->initial_dg_values = t8dg_advect_diff_problem_get_dg_values (problem);

  /* This constructs new dg_values, a new functionbasis, new globa_values and new local_values; given the new order of LGL points */
  mg_ctx->p_mg_mat_ctx->coarse_dg_values =
    t8dg_values_new_LGL_hypercube (t8dg_values_get_dim (mg_ctx->p_mg_mat_ctx->initial_dg_values), mg_ctx->p_mg_mat_ctx->mg_coarse_order,
                                   t8dg_values_get_coarse_geometry (mg_ctx->p_mg_mat_ctx->initial_dg_values),
                                   t8dg_values_get_forest (mg_ctx->p_mg_mat_ctx->initial_dg_values));
  /* Also, construct corresponding dof_values */
  mg_ctx->p_mg_mat_ctx->coarse_lvl_dofs =
    t8dg_dof_values_new (t8dg_values_get_forest (mg_ctx->p_mg_mat_ctx->initial_dg_values),
                         t8dg_values_get_global_values_array (mg_ctx->p_mg_mat_ctx->coarse_dg_values));

  /* Get the number of degrees of freedom on the coarse level as well as on the fine level */
  mg_ctx->p_mg_mat_ctx->num_coarse_dofs =
    (size_t) t8dg_dof_get_num_local_elements (mg_ctx->p_mg_mat_ctx->coarse_lvl_dofs) *
    t8dg_dof_get_max_num_element_dof (mg_ctx->p_mg_mat_ctx->coarse_lvl_dofs);
  mg_ctx->p_mg_mat_ctx->num_fine_dofs =
    (size_t) t8dg_dof_get_num_local_elements (*(mg_ctx->p_mg_mat_ctx->problem_dofs)) *
    t8dg_dof_get_max_num_element_dof (*(mg_ctx->p_mg_mat_ctx->problem_dofs));
}

PetscErrorCode
t8dg_precon_destroy_p_mg (t8dg_p_mg_lvl_ctx_t ** mg_ctx)
{
  PetscErrorCode      ierr;

  ierr = MatDestroy (&((*mg_ctx)->Restriction));
  CHKERRQ (ierr);
  ierr = MatDestroy (&((*mg_ctx)->Prolongation));
  CHKERRQ (ierr);
  ierr = MatDestroy (&((*mg_ctx)->Smoothing_Mat));
  CHKERRQ (ierr);
  ierr = MatDestroy (&((*mg_ctx)->Coarse_Mat));
  CHKERRQ (ierr);
  t8dg_dof_values_destroy (&((*mg_ctx)->p_mg_mat_ctx->coarse_lvl_dofs));
  t8dg_values_destroy (&((*mg_ctx)->p_mg_mat_ctx->coarse_dg_values));

  T8DG_FREE ((*mg_ctx)->p_mg_mat_ctx);

  return 0;
}
#endif
/******************/
/******************/
/* End of p MG */
/******************/
/******************/

/******************/
/******************/
/* BEGIN Multiple Lvls MG */
/*****************/
/*****************/
/* Initialize a multigrid preconditioner with various coarse levels */
PetscErrorCode
t8dg_precon_init_multiple_level_mg (void *problem, t8dg_dof_values_t ** problem_dofs,
                                    t8dg_time_derivation_matrix_application_t time_derivative, Mat * A_fine, PC * mg_pc,
                                    PetscInt * vec_global_index, int num_mg_levels, t8dg_mg_levels_ctx_t ** mg_ctx)
{
  /* the permitted max. number of coarse levels is defined in 't8dg_preconditioner.h' */
  T8DG_ASSERT (T8DG_PRECON_MAX_MG_LVLS > num_mg_levels);
  if (T8DG_PRECON_MAX_MG_LVLS <= num_mg_levels) {
    num_mg_levels = T8DG_PRECON_MAX_MG_LVLS - 1;
    t8dg_global_essentialf ("Number of Multigrid Levels got corrected to %d\n", num_mg_levels);
  }
  PetscErrorCode      ierr;

  /* Construct all multigrid levels */
  t8dg_mg_set_up_multiple_lvl_precon ((t8dg_linear_advection_diffusion_problem_t *) problem, problem_dofs, time_derivative,
                                      vec_global_index, num_mg_levels, *mg_ctx);

  /* Construct Restriction Operators */
  t8dg_mg_set_up_restriction_operators (*mg_ctx);

  /* Construct Prolongation Operators */
  t8dg_mg_set_up_prolongation_operators (*mg_ctx);

  /* Construct the coarse level solver */
  t8dg_mg_set_up_coarse_lvl_matrix (*mg_ctx);

  /* Construct Smoothing Operators */
  t8dg_mg_set_up_smoothing_operators (*mg_ctx, A_fine);

  /* Set Coarse Solver, Smoother, Matrices and Operators */
  /* Set the Multigrid Preconditioner */
  ierr = PCSetType (*mg_pc, PCMG);
  CHKERRQ (ierr);

  /* Set the nnumber of multigrid levels; the initial/fine level is included in this number */
  ierr = PCMGSetLevels (*mg_pc, num_mg_levels, NULL);
  CHKERRQ (ierr);

  /* Choose a V-Cycle routine for the Preconditioner */
  ierr = PCMGSetType (*mg_pc, PC_MG_MULTIPLICATIVE);
  CHKERRQ (ierr);

  /* Set the cycle rythm; alternatively, 'PC_MG_CYCLE_W' can be used in order to perfom a W-Cycle MG preconditioner */
  ierr = PCMGSetCycleType (*mg_pc, T8DG_PRECON_MG_CYCLE_TYPE);
  CHKERRQ (ierr);

  /* Assign the coarse level solver */
  ierr = PCMGGetCoarseSolve (*mg_pc, &(*mg_ctx)->coarse_solver);
  CHKERRQ (ierr);

  /* Set the operators for the coarse level solver and assign the corresponding matrix-free application */
  ierr = KSPSetOperators ((*mg_ctx)->coarse_solver, (*mg_ctx)->A_coarse_Mat, (*mg_ctx)->A_coarse_Mat);
  CHKERRQ (ierr);

  /* Choose default values for convergence tolerances on the coarse level */
  ierr = KSPSetTolerances ((*mg_ctx)->coarse_solver, T8DG_PRECON_MG_GMRES_ACC, T8DG_PRECON_MG_GMRES_ACC, PETSC_DEFAULT, PETSC_DEFAULT);
  CHKERRQ (ierr);

  /* Choose a GMRES routine as the coarse level solver */
  ierr = KSPSetType ((*mg_ctx)->coarse_solver, KSPGMRES);
  CHKERRQ (ierr);

  /* Extract the preconditioner on the coarse level */
  ierr = KSPGetPC ((*mg_ctx)->coarse_solver, &(*mg_ctx)->coarse_pc);
  CHKERRQ (ierr);

#if T8DG_PRECON_MG_COARSE_GRID_PRECONDITIONER
  /* Use the Block-Jacobi preconditioner at the coarsest level in order to solve the coarse grid correction */
  /* Currently, this results in no efficient use-cases -> rather use 'without_preconditioning */
  t8dg_precon_init_jacobi (problem, &(((*mg_ctx)->interpolation_ctx->dofs_lvl)[(*mg_ctx)->interpolation_ctx->num_mg_levels - 1]),
                           &((*mg_ctx)->A_coarse_Mat), &((*mg_ctx)->coarse_pc), &((*mg_ctx)->precon_jacobi_ctx),
                           ((*mg_ctx)->interpolation_ctx->indexing_scheme)[(*mg_ctx)->interpolation_ctx->num_mg_levels - 1], 1);
#else
  /* Set no preconditioning on the coarse level */
  t8dg_precon_init_without_preconditioning (&(*mg_ctx)->coarse_pc);
#endif

  /* Assign the smoothers on each multigrid level */
  t8dg_mg_initialize_smoothers (mg_pc, mg_ctx);

  /* Amount of pre- and post-smoothing steps on each multigrid level (except the coarsest) (seems like N+1 smoothing steps are performed) */
  ierr = PCMGSetNumberSmooth (*mg_pc, T8DG_PRECON_MG_NUM_SMOOTH_ITERATIONS);
  CHKERRQ (ierr);

  /* Assign the restriction operators between the different levels */
  t8dg_mg_initialize_restrictions (mg_pc, mg_ctx);

  /* Assign the prolongation operators between different levels */
  t8dg_mg_initialize_prolongations (mg_pc, mg_ctx);

  /* Residuals and three further vectors on each level needed by the multigrid preconditioner are automatically handled by PETSc */

  /* Set the application time to zero */
  (*mg_ctx)->mg_general_data->preconditioner_application_time_coarse_lvl = 0.0;
  (*mg_ctx)->mg_general_data->preconditioner_application_time_smoothing = 0.0;
  (*mg_ctx)->mg_general_data->preconditioner_application_time_interpolation = 0.0;

  return 0;
}

/* Assign matrix-free applications to the Prolongation matrices */
PetscErrorCode
t8dg_mg_initialize_prolongations (PC * mg_pc, t8dg_mg_levels_ctx_t ** mg_ctx)
{
  PetscErrorCode      ierr;
  int                 lvl_iter;

  for (lvl_iter = (*mg_ctx)->interpolation_ctx->num_mg_levels - 1; lvl_iter > 0; --lvl_iter) {
    /* Set the Prolongation operator from the coarse level onto the fine level */
    ierr =
      PCMGSetInterpolation (*mg_pc, lvl_iter, (*mg_ctx)->Prolongation_Mats[(*mg_ctx)->interpolation_ctx->num_mg_levels - lvl_iter - 1]);
    CHKERRQ (ierr);
  }
  return 0;
}

/* Assign matrix-free applications to the Restriction matrices */
PetscErrorCode
t8dg_mg_initialize_restrictions (PC * mg_pc, t8dg_mg_levels_ctx_t ** mg_ctx)
{
  PetscErrorCode      ierr;
  int                 lvl_iter;

  for (lvl_iter = (*mg_ctx)->interpolation_ctx->num_mg_levels - 1; lvl_iter > 0; --lvl_iter) {
    ierr = PCMGSetRestriction (*mg_pc, lvl_iter, (*mg_ctx)->Restriction_Mats[(*mg_ctx)->interpolation_ctx->num_mg_levels - lvl_iter - 1]);
    CHKERRQ (ierr);
  }
  return 0;
}

/* Set up smoothing operators on each multigrid level (except the coarsest) */
PetscErrorCode
t8dg_mg_set_up_smoothing_operators (t8dg_mg_levels_ctx_t * mg_lvls, Mat * A_fine)
{
  PetscErrorCode      ierr;
  int                 lvl_iter;
  /* Initialize smoothing matrix on the finest grid; this corresponds to the 'normal' system matrix of the problem */
  mg_lvls->Smoothing_Mats[0] = *A_fine;

  /* Assign smoothing matrices for all multigrid levels; except the coarsest */
  for (lvl_iter = 1; lvl_iter < mg_lvls->interpolation_ctx->num_mg_levels - 1; ++lvl_iter) {
    /* Smoothing matrices need the same information as the interpolation matrices, therefore, the same context will be assigned */
    ierr =
      MatCreateShell (PETSC_COMM_WORLD, (mg_lvls->interpolation_ctx->num_local_dofs[lvl_iter]),
                      (mg_lvls->interpolation_ctx->num_local_dofs[lvl_iter]), PETSC_DETERMINE, PETSC_DETERMINE, mg_lvls->interpolation_ctx,
                      &(mg_lvls->Smoothing_Mats[lvl_iter]));
    CHKERRQ (ierr);

    /* Define the (multiplicative) MatVec-Operation which resembles the application of the amoothing matrix */
    ierr = MatShellSetOperation ((mg_lvls->Smoothing_Mats[lvl_iter]), MATOP_MULT, (void (*)(void)) MatMult_MF_MG_LVL_Smoothing);
    CHKERRQ (ierr);
  }
  return 0;
}

/* Matrix-free application which resembles the smoothing on a multigrid level */
PetscErrorCode
MatMult_MF_MG_LVL_Smoothing (Mat Smoother, Vec in, Vec out)
{
  PetscErrorCode      ierr;
  t8dg_mg_lvl_interpolation_ctx_t *interpolation_ctx;

  /* Get the matrix-free application context */
  ierr = MatShellGetContext (Smoother, &interpolation_ctx);
  CHKERRQ (ierr);

  /* Measure application time */
#if T8DG_PRECON_MEASURE_PC_TIME
  interpolation_ctx->mg_general_data->preconditioner_application_time_smoothing -= sc_MPI_Wtime ();
#endif

  /* Copy the petsc vec entries into the dof_values of the problem in order to compute their derivation */
  t8dg_precon_write_vec_to_dof (&in, interpolation_ctx->dofs_lvl[interpolation_ctx->current_lvl],
                                interpolation_ctx->num_local_dofs[interpolation_ctx->current_lvl]);

  /* Swap dofs with initial problem_dofs */
  t8dg_dof_values_swap (interpolation_ctx->mg_general_data->pdof_array, &(interpolation_ctx->dofs_lvl[interpolation_ctx->current_lvl]));

  /* Calculate the time derivation of the degrees of freedom */
  (interpolation_ctx->mg_general_data->time_derivative_func) (*(interpolation_ctx->mg_general_data->pdof_array),
                                                              interpolation_ctx->dofs_lvl_derivation[interpolation_ctx->current_lvl - 1],
                                                              interpolation_ctx->mg_general_data->current_point_in_time,
                                                              interpolation_ctx->mg_general_data->problem);

  /* Swap dof values back */
  t8dg_dof_values_swap (interpolation_ctx->mg_general_data->pdof_array, &(interpolation_ctx->dofs_lvl[interpolation_ctx->current_lvl]));

  /* calculate the sum of the initially passed degrees of freedom ('in' Vector) and their (scaled) derivation */
  t8dg_dof_values_axpy (-interpolation_ctx->mg_general_data->timestep * interpolation_ctx->mg_general_data->current_coefficient,
                        interpolation_ctx->dofs_lvl_derivation[interpolation_ctx->current_lvl - 1],
                        interpolation_ctx->dofs_lvl[interpolation_ctx->current_lvl]);

  /* Copy the dof_values back to a petsc_vector */
  t8dg_precon_write_dof_to_vec (interpolation_ctx->dofs_lvl[interpolation_ctx->current_lvl], &out,
                                interpolation_ctx->indexing_scheme[interpolation_ctx->current_lvl],
                                interpolation_ctx->num_local_dofs[interpolation_ctx->current_lvl]);

  /* Measure application time */
#if T8DG_PRECON_MEASURE_PC_TIME
  interpolation_ctx->mg_general_data->preconditioner_application_time_smoothing += sc_MPI_Wtime ();
#endif

  return 0;
}

/* Initialize the smoother on each coarse level (except on the coarsest) */
PetscErrorCode
t8dg_mg_initialize_smoothers (PC * mg_pc, t8dg_mg_levels_ctx_t ** mg_ctx)
{
  PetscErrorCode      ierr;
  int                 lvl_iter;
  int                 num_mg_smooth_levels = (*mg_ctx)->interpolation_ctx->num_mg_levels - 1;

  for (lvl_iter = 1; lvl_iter < (*mg_ctx)->interpolation_ctx->num_mg_levels; ++lvl_iter) {
    /* Assign the smoother of the multigrid preconditioner; 0 corresponds to the coarsest level; its reversed in terms of the rest of the implementation of the preconditioner */
    /* on the coarsest level, there is no smoothing performed */
    ierr = PCMGGetSmoother (*mg_pc, lvl_iter, &(*mg_ctx)->smoothers[num_mg_smooth_levels - lvl_iter]);
    CHKERRQ (ierr);
    ierr =
      KSPSetOperators ((*mg_ctx)->smoothers[num_mg_smooth_levels - lvl_iter], (*mg_ctx)->Smoothing_Mats[num_mg_smooth_levels - lvl_iter],
                       (*mg_ctx)->Smoothing_Mats[num_mg_smooth_levels - lvl_iter]);
    CHKERRQ (ierr);
    /* Choose a GMRES routine as the smooter on the fine level */
    ierr = KSPSetType ((*mg_ctx)->smoothers[num_mg_smooth_levels - lvl_iter], KSPGMRES);
    CHKERRQ (ierr);
    /* Extract the preconditioner of the smoother */
    ierr = KSPGetPC ((*mg_ctx)->smoothers[num_mg_smooth_levels - lvl_iter], &((*mg_ctx)->smoother_pcs[num_mg_smooth_levels - lvl_iter]));
    CHKERRQ (ierr);
    /* Set no preconditioning within the smoothing iterations */
    t8dg_precon_init_without_preconditioning (&(*mg_ctx)->smoother_pcs[num_mg_smooth_levels - lvl_iter]);

  }
  return 0;
}

/* Matrix-free routines which performs the restriction between consecutive levels */
PetscErrorCode
MatMult_MF_MG_LVL_Restriction (Mat Restriction, Vec in, Vec out)
{
  PetscErrorCode      ierr;
  t8dg_mg_lvl_interpolation_ctx_t *interpolation_ctx;

  /* Get the matrix-free application context */
  ierr = MatShellGetContext (Restriction, &interpolation_ctx);
  CHKERRQ (ierr);

  /* Measure application time */
#if T8DG_PRECON_MEASURE_PC_TIME
  interpolation_ctx->mg_general_data->preconditioner_application_time_interpolation -= sc_MPI_Wtime ();
#endif

  /* Assign the pointers in adapt_data with a proper size according to the current multigrid level */
  interpolation_ctx->adapt_data->dof_values = interpolation_ctx->dofs_lvl[interpolation_ctx->current_lvl];
  interpolation_ctx->adapt_data->dof_values_adapt = interpolation_ctx->dofs_lvl[interpolation_ctx->current_lvl + 1];

  /* Prepare dg_values for the restriction step */
  t8dg_values_mg_lvl_set_interpolation_step (interpolation_ctx->dg_values,
                                             interpolation_ctx->coarse_level_forests[interpolation_ctx->current_lvl],
                                             interpolation_ctx->coarse_level_forests[interpolation_ctx->current_lvl + 1],
                                             interpolation_ctx->local_values_lvl[interpolation_ctx->current_lvl],
                                             interpolation_ctx->local_values_lvl[interpolation_ctx->current_lvl + 1]);

  /* Write the entries of the 'in' vector in adapt_data->dof_values, these are getting adapted/restricted */
  t8dg_precon_write_vec_to_dof (&in, interpolation_ctx->adapt_data->dof_values,
                                interpolation_ctx->num_local_dofs[interpolation_ctx->current_lvl]);

  /* Calling iterate_replace executes the adaption/restriction */
  t8_forest_iterate_replace (interpolation_ctx->coarse_level_forests[interpolation_ctx->current_lvl + 1],
                             interpolation_ctx->coarse_level_forests[interpolation_ctx->current_lvl], t8dg_adapt_replace);

  /* Write the degrees of freedom to a petsc vec */
  t8dg_precon_write_dof_to_vec (interpolation_ctx->adapt_data->dof_values_adapt, &out,
                                interpolation_ctx->indexing_scheme[interpolation_ctx->current_lvl + 1],
                                interpolation_ctx->num_local_dofs[interpolation_ctx->current_lvl + 1]);

  /* Prepare the next restriction step; currently: set the mortar array */
  t8dg_values_mg_lvl_prepare_next_interpolation_step (interpolation_ctx->dg_values,
                                                      interpolation_ctx->mortar_array_lvl[interpolation_ctx->current_lvl + 1]);

  /* Update the multigrid level index */
  ++(interpolation_ctx->current_lvl);

  /* Measure application time */
#if T8DG_PRECON_MEASURE_PC_TIME
  interpolation_ctx->mg_general_data->preconditioner_application_time_interpolation += sc_MPI_Wtime ();
#endif

  return 0;

}

/* matrix-free routine which performs the prolongation between two consecutive levels */
PetscErrorCode
MatMult_MF_MG_LVL_Prolongation (Mat Prolongation, Vec in, Vec out)
{
  /* The first prolongation step is called after all restriction steps have been called, therefore the current_level within the context has the value (num_mg_levels -1)  */
  PetscErrorCode      ierr;
  t8dg_mg_lvl_interpolation_ctx_t *interpolation_ctx;

  /* Get the matrix-free application context */
  ierr = MatShellGetContext (Prolongation, &interpolation_ctx);
  CHKERRQ (ierr);

  /* Measure application time */
#if T8DG_PRECON_MEASURE_PC_TIME
  interpolation_ctx->mg_general_data->preconditioner_application_time_interpolation -= sc_MPI_Wtime ();
#endif

  /* Assign pointers with a proper size according to the current multigrid level */
  interpolation_ctx->adapt_data->dof_values = interpolation_ctx->dofs_lvl[interpolation_ctx->current_lvl];
  interpolation_ctx->adapt_data->dof_values_adapt = interpolation_ctx->dofs_lvl[interpolation_ctx->current_lvl - 1];

  /* Prepare dg_values for the prolongation step */
  t8dg_values_mg_lvl_set_interpolation_step (interpolation_ctx->dg_values,
                                             interpolation_ctx->coarse_level_forests[interpolation_ctx->current_lvl],
                                             interpolation_ctx->coarse_level_forests[interpolation_ctx->current_lvl - 1],
                                             interpolation_ctx->local_values_lvl[interpolation_ctx->current_lvl],
                                             interpolation_ctx->local_values_lvl[interpolation_ctx->current_lvl - 1]);

  /* Write the entries of the 'in' vector in adapt_data->dof_values, these are getting adapted/prolongated */
  t8dg_precon_write_vec_to_dof (&in, interpolation_ctx->adapt_data->dof_values,
                                interpolation_ctx->num_local_dofs[interpolation_ctx->current_lvl]);

  /* Calling iterate_replace executes the adaption/prolongation */
  t8_forest_iterate_replace (interpolation_ctx->coarse_level_forests[interpolation_ctx->current_lvl - 1],
                             interpolation_ctx->coarse_level_forests[interpolation_ctx->current_lvl], t8dg_adapt_replace);

  /* Write the degrees of freedom to a petsc vec */
  t8dg_precon_write_dof_to_vec (interpolation_ctx->adapt_data->dof_values_adapt, &out,
                                interpolation_ctx->indexing_scheme[interpolation_ctx->current_lvl - 1],
                                interpolation_ctx->num_local_dofs[interpolation_ctx->current_lvl - 1]);

  /* Prepare the next restriction step; currently: set the mortar array */
  t8dg_values_mg_lvl_prepare_next_interpolation_step (interpolation_ctx->dg_values,
                                                      interpolation_ctx->mortar_array_lvl[interpolation_ctx->current_lvl - 1]);

  /* Update the multigrid level index */
  --(interpolation_ctx->current_lvl);

  /* Measure application time */
#if T8DG_PRECON_MEASURE_PC_TIME
  interpolation_ctx->mg_general_data->preconditioner_application_time_interpolation += sc_MPI_Wtime ();
#endif

  return 0;
}

PetscErrorCode
MatMult_MF_MG_LVL_Coarse_Solver (Mat C_Mat, Vec in, Vec out)
{
  /* The restriction steps left the dg_values initialized just right on the coarseest level, therefore, without preparation, the coarse grid correction can be calculated */

  PetscErrorCode      ierr;
  t8dg_mg_lvl_coarse_matrix_ctx_t *c_appctx;

  /* Get the matrix-free application context */
  ierr = MatShellGetContext (C_Mat, &c_appctx);
  CHKERRQ (ierr);

  /* Measure application time */
#if T8DG_PRECON_MEASURE_PC_TIME
  c_appctx->mg_general_data->preconditioner_application_time_coarse_lvl -= sc_MPI_Wtime ();
#endif
  /* Copy the petsc vec entries into the dof_values of the problem in order to compute their derivation */
  t8dg_precon_write_vec_to_dof (&in, *(c_appctx->coarse_lvl_problem_dofs), c_appctx->num_local_dofs_coarsest_lvl);

  /* Swap dofs with initial problem_dofs */
  t8dg_dof_values_swap (c_appctx->mg_general_data->pdof_array, c_appctx->coarse_lvl_problem_dofs);

  (c_appctx->mg_general_data->time_derivative_func) (*(c_appctx->mg_general_data->pdof_array), c_appctx->coarse_lvl_problem_dofs_derivation,
                                                     c_appctx->mg_general_data->current_point_in_time, c_appctx->mg_general_data->problem);

  /* Swap dof values back */
  t8dg_dof_values_swap (c_appctx->mg_general_data->pdof_array, c_appctx->coarse_lvl_problem_dofs);

  /* calculate the sum of the initially passed degrees of freedom ('in' Vector) and their (scaled) derivation */
  t8dg_dof_values_axpy (-c_appctx->mg_general_data->timestep * c_appctx->mg_general_data->current_coefficient,
                        c_appctx->coarse_lvl_problem_dofs_derivation, *(c_appctx->coarse_lvl_problem_dofs));

  /* Copy the dof_values back to a petsc_vector */
  t8dg_precon_write_dof_to_vec (*(c_appctx->coarse_lvl_problem_dofs), &out, c_appctx->indexing_coarsest_level,
                                c_appctx->num_local_dofs_coarsest_lvl);

  /* Measure application time */
#if T8DG_PRECON_MEASURE_PC_TIME
  c_appctx->mg_general_data->preconditioner_application_time_coarse_lvl += sc_MPI_Wtime ();
#endif

  return 0;
}

/* Create the restriction matrices */
PetscErrorCode
t8dg_mg_set_up_restriction_operators (t8dg_mg_levels_ctx_t * mg_lvls)
{
  PetscErrorCode      ierr;
  int                 num_lvl;
  for (num_lvl = 0; num_lvl < mg_lvls->interpolation_ctx->num_mg_levels - 1; ++num_lvl) {
    /* Create the restriction matrix on each level */
    ierr =
      MatCreateShell (PETSC_COMM_WORLD, (mg_lvls->interpolation_ctx->num_local_dofs[num_lvl + 1]),
                      (mg_lvls->interpolation_ctx->num_local_dofs[num_lvl]), PETSC_DETERMINE, PETSC_DETERMINE, mg_lvls->interpolation_ctx,
                      &(mg_lvls->Restriction_Mats[num_lvl]));
    CHKERRQ (ierr);

    /* Define the (multiplicative) MatVec-Operation which resembles the application of the restriction matrix */
    ierr = MatShellSetOperation ((mg_lvls->Restriction_Mats[num_lvl]), MATOP_MULT, (void (*)(void)) MatMult_MF_MG_LVL_Restriction);
    CHKERRQ (ierr);
  }
  return 0;
}

/* Create the prolongation matrices */
PetscErrorCode
t8dg_mg_set_up_prolongation_operators (t8dg_mg_levels_ctx_t * mg_lvls)
{
  PetscErrorCode      ierr;
  int                 num_lvl;

  for (num_lvl = 0; num_lvl < mg_lvls->interpolation_ctx->num_mg_levels - 1; ++num_lvl) {
    /* Create the restriction matrix on each level */
    ierr =
      MatCreateShell (PETSC_COMM_WORLD, mg_lvls->interpolation_ctx->num_local_dofs[num_lvl],
                      mg_lvls->interpolation_ctx->num_local_dofs[num_lvl + 1], PETSC_DETERMINE, PETSC_DETERMINE, mg_lvls->interpolation_ctx,
                      &(mg_lvls->Prolongation_Mats[num_lvl]));
    CHKERRQ (ierr);

    /* Define the (multiplicative) MatVec-Operation which resembles the application of the restriction matrix */
    ierr = MatShellSetOperation ((mg_lvls->Prolongation_Mats[num_lvl]), MATOP_MULT, (void (*)(void)) MatMult_MF_MG_LVL_Prolongation);
    CHKERRQ (ierr);
  }
  return 0;
}

/* Create coarse level matrix (on the coarsest level) */
PetscErrorCode
t8dg_mg_set_up_coarse_lvl_matrix (t8dg_mg_levels_ctx_t * mg_lvls)
{
  PetscErrorCode      ierr;

  ierr =
    MatCreateShell (PETSC_COMM_WORLD, mg_lvls->interpolation_ctx->num_local_dofs[mg_lvls->interpolation_ctx->num_mg_levels - 1],
                    mg_lvls->interpolation_ctx->num_local_dofs[mg_lvls->interpolation_ctx->num_mg_levels - 1], PETSC_DETERMINE,
                    PETSC_DETERMINE, mg_lvls->coarse_matrix_ctx, &(mg_lvls->A_coarse_Mat));
  CHKERRQ (ierr);

  /* Define the (multiplicative) MatVec-Operation on the coarse level mesh */
  ierr = MatShellSetOperation ((mg_lvls->A_coarse_Mat), MATOP_MULT, (void (*)(void)) MatMult_MF_MG_LVL_Coarse_Solver);
  CHKERRQ (ierr);
  return 0;
}

/* A function that sets up the components needed in order to perform a multigrid V-Cycle on various levels */
void
t8dg_mg_set_up_multiple_lvl_precon (t8dg_linear_advection_diffusion_problem_t * problem, t8dg_dof_values_t ** pdof_array,
                                    t8dg_time_matrix_application time_derivative, PetscInt * initial_forest_indexing, int num_mg_levels,
                                    t8dg_mg_levels_ctx_t * mg_ctx)
{
  int                 num_lvl;
  char                file_name_forest[40];

  /* Allocate space for a general data struct */
  mg_ctx->mg_general_data = T8DG_ALLOC_ZERO (t8dg_mg_general_data_t, 1);
  /* Initialize the general data properties */
  mg_ctx->mg_general_data->problem = (void *) problem;
  mg_ctx->mg_general_data->time_derivative_func = time_derivative;
  mg_ctx->mg_general_data->current_point_in_time =
    t8dg_timestepping_data_get_current_time (t8dg_advect_diff_problem_get_time_data (problem));
  mg_ctx->mg_general_data->timestep = t8dg_timestepping_data_get_time_step (t8dg_advect_diff_problem_get_time_data (problem));
  mg_ctx->mg_general_data->current_coefficient = 1.0;   /* in case of implicit euler */
  mg_ctx->mg_general_data->pdof_array = pdof_array;

  mg_ctx->interpolation_ctx = T8DG_ALLOC_ZERO (t8dg_mg_lvl_interpolation_ctx_t, 1);
  mg_ctx->interpolation_ctx->current_lvl = 0;
  mg_ctx->interpolation_ctx->num_mg_levels = num_mg_levels;
  mg_ctx->interpolation_ctx->adapt_data = t8dg_advect_diff_problem_get_adapt_data (problem);
  mg_ctx->interpolation_ctx->dg_values = t8dg_advect_diff_problem_get_dg_values (problem);
  mg_ctx->interpolation_ctx->adapt_data->dg_values = mg_ctx->interpolation_ctx->dg_values;
  mg_ctx->interpolation_ctx->mg_general_data = mg_ctx->mg_general_data;

  /* Initialize members with initial problem */
  (mg_ctx->interpolation_ctx->coarse_level_forests)[0] = t8dg_advect_diff_problem_get_forest (problem);
  t8_forest_set_user_data (mg_ctx->interpolation_ctx->coarse_level_forests[0], mg_ctx->interpolation_ctx->adapt_data);
  t8dg_adapt_data_set_time (mg_ctx->interpolation_ctx->adapt_data, mg_ctx->mg_general_data->current_point_in_time);
  mg_ctx->interpolation_ctx->indexing_scheme[0] = initial_forest_indexing;
  mg_ctx->interpolation_ctx->num_local_dofs[0] =
    (size_t) t8dg_dof_get_num_local_elements (*pdof_array) * t8dg_dof_get_max_num_element_dof (*pdof_array);
  mg_ctx->interpolation_ctx->dofs_lvl[0] = t8dg_dof_values_duplicate (*pdof_array);
  mg_ctx->interpolation_ctx->local_values_lvl[0] = *t8dg_values_get_local_values (mg_ctx->interpolation_ctx->dg_values);
  mg_ctx->interpolation_ctx->mortar_array_lvl[0] = *t8dg_values_get_mortar_array (mg_ctx->interpolation_ctx->dg_values);

  /* Construct each coarse level */
  for (num_lvl = 1; num_lvl < num_mg_levels; ++num_lvl) {
    t8dg_debugf ("Construction of coarse forest on level %d has been called\n", num_lvl);
    t8dg_mg_construct_coarse_level_forest (mg_ctx, t8dg_adapt_multigrid_coarsen_finest_level, num_lvl);
  }

  /* Write out the created forests */
#if T8DG_PRECON_MG_WRITE_OUT_COARSE_LVL
  /* write out the initial mesh */
  t8_forest_write_vtk ((mg_ctx->interpolation_ctx->coarse_level_forests)[0], "t8dg_mg_fine_lvl_forest");
  /* write out the coarse level meshes */
  for (num_lvl = 1; num_lvl < num_mg_levels; ++num_lvl) {
    sprintf (file_name_forest, "t8dg_mg_coarse_lvl_forest_%d", num_lvl);
    t8_forest_write_vtk ((mg_ctx->interpolation_ctx->coarse_level_forests)[num_lvl], file_name_forest);
  }
#endif

  /* Initialize members which are based on the constructed coarse level forests */
  mg_ctx->coarse_matrix_ctx = T8DG_ALLOC_ZERO (t8dg_mg_lvl_coarse_matrix_ctx_t, 1);
  mg_ctx->coarse_matrix_ctx->num_local_dofs_coarsest_lvl = mg_ctx->interpolation_ctx->num_local_dofs[num_mg_levels - 1];
  mg_ctx->coarse_matrix_ctx->mg_general_data = mg_ctx->mg_general_data;
  mg_ctx->coarse_matrix_ctx->coarse_lvl_problem_dofs = &(mg_ctx->interpolation_ctx->dofs_lvl[num_mg_levels - 1]);
  mg_ctx->coarse_matrix_ctx->coarse_lvl_problem_dofs_derivation = mg_ctx->interpolation_ctx->dofs_lvl_derivation[num_mg_levels - 2];
  mg_ctx->coarse_matrix_ctx->indexing_coarsest_level = mg_ctx->interpolation_ctx->indexing_scheme[num_mg_levels - 1];

  /* Allocate space needed by the restriction and interpolation steps */
  t8dg_values_mg_lvl_allocate_properties (mg_ctx->interpolation_ctx->dg_values, mg_ctx->interpolation_ctx->num_mg_levels,
                                          mg_ctx->interpolation_ctx->coarse_level_forests, mg_ctx->interpolation_ctx->local_values_lvl,
                                          mg_ctx->interpolation_ctx->mortar_array_lvl);
}

/* Constructs a coarse level */
/* TODO: Eventually source and sink terms need to be considered */
PetscErrorCode
t8dg_mg_construct_coarse_level_forest (t8dg_mg_levels_ctx_t * mg_lvls, t8dg_corase_lvl_adapt_func_t coarsening_func, int num_lvl)
{
  PetscErrorCode      ierr;
  t8_locidx_t         global_offset_to_first_local_elem;
  size_t              iter;

  /* Keep the forest of the former level */
  t8_forest_ref (mg_lvls->interpolation_ctx->coarse_level_forests[num_lvl - 1]);

  /* Initialize a coarsened forest conatining the elements of the coarser muligrid mesh */
  t8_forest_init (&(mg_lvls->interpolation_ctx->coarse_level_forests[num_lvl]));

  /* Set the user-data pointer of the new forest */
  t8_forest_set_user_data (mg_lvls->interpolation_ctx->coarse_level_forests[num_lvl], mg_lvls->interpolation_ctx->adapt_data);

  /* Set the adapt function coarsening the forest */
  t8_forest_set_adapt (mg_lvls->interpolation_ctx->coarse_level_forests[num_lvl],
                       mg_lvls->interpolation_ctx->coarse_level_forests[num_lvl - 1], coarsening_func, 0);

  /* eventually Balance the forest */
  if (mg_lvls->interpolation_ctx->adapt_data->maximum_refinement_level - mg_lvls->interpolation_ctx->adapt_data->minimum_refinement_level >
      1) {
    t8_forest_set_balance (mg_lvls->interpolation_ctx->coarse_level_forests[num_lvl], NULL, 1);
  }

  /* Ghost values are needed for solving the system on the coarser mesh */
  t8_forest_set_ghost (mg_lvls->interpolation_ctx->coarse_level_forests[num_lvl], 1, T8_GHOST_FACES);

  /* Commit the pre-set forest, so that it will be adapted */
  t8_forest_commit (mg_lvls->interpolation_ctx->coarse_level_forests[num_lvl]);

  /* Allocate space for the degrees of freedom on this level */
  mg_lvls->interpolation_ctx->dofs_lvl[num_lvl] =
    t8dg_dof_values_new (mg_lvls->interpolation_ctx->coarse_level_forests[num_lvl],
                         t8dg_values_get_global_values_array (mg_lvls->interpolation_ctx->dg_values));
  mg_lvls->interpolation_ctx->dofs_lvl_derivation[num_lvl - 1] = t8dg_dof_values_duplicate (mg_lvls->interpolation_ctx->dofs_lvl[num_lvl]);

  /* Get the number of degrees of freedom on this level */
  (mg_lvls->interpolation_ctx->num_local_dofs)[num_lvl] =
    (size_t) t8dg_dof_get_num_local_elements (mg_lvls->interpolation_ctx->dofs_lvl[num_lvl]) *
    t8dg_dof_get_max_num_element_dof (mg_lvls->interpolation_ctx->dofs_lvl[num_lvl]);

  /* Create an indexing scheme on this level for the PETSc Vecs */
  /* Allocate space for an array which holds the global indices of the degrees of freedom */
  ierr = PetscMalloc1 (mg_lvls->interpolation_ctx->num_local_dofs[num_lvl], &(mg_lvls->interpolation_ctx->indexing_scheme[num_lvl]));
  CHKERRQ (ierr);

  /* Compute the offset of the first local element concerning the forest and the amount of degrees of freedom per element */
  global_offset_to_first_local_elem = t8_forest_get_first_local_element_id (mg_lvls->interpolation_ctx->coarse_level_forests[num_lvl]);
  if (global_offset_to_first_local_elem != 0) {
    global_offset_to_first_local_elem =
      global_offset_to_first_local_elem * t8dg_dof_get_max_num_element_dof (mg_lvls->interpolation_ctx->dofs_lvl[num_lvl]);
  }

  /* Fill the array of global indices */
  for (iter = 0; iter < mg_lvls->interpolation_ctx->num_local_dofs[num_lvl]; ++iter) {
    (mg_lvls->interpolation_ctx->indexing_scheme[num_lvl])[iter] = iter + global_offset_to_first_local_elem;
  }

  t8dg_debugf ("Der neue forest hat %d Elemente\nUnd die Anzahl Freiheitsgrade sind %d \n",
               t8_forest_get_num_element (((mg_lvls)->interpolation_ctx->coarse_level_forests)[num_lvl]),
               ((mg_lvls)->interpolation_ctx->num_local_dofs)[num_lvl]);

  return 0;
}

/* Destroys the multigrid preconditioner and frees the allocated space */
PetscErrorCode
t8dg_precon_destroy_mg_levels (t8dg_mg_levels_ctx_t ** mg_ctx)
{
  int                 num_lvl;
  PetscErrorCode      ierr;

  /* Free the general data */
  T8DG_FREE ((*mg_ctx)->mg_general_data);

  /* Free coarse level preconditioner */
#if T8DG_PRECON_MG_COARSE_GRID_PRECONDITIONER
  t8dg_precon_destroy_block_preconditioner ((*mg_ctx)->precon_jacobi_ctx);
#endif

  /* Free all allocated space */
  t8dg_dof_values_destroy (&((*mg_ctx)->interpolation_ctx->dofs_lvl[0]));
  for (num_lvl = 1; num_lvl < (*mg_ctx)->interpolation_ctx->num_mg_levels; ++num_lvl) {
    t8_forest_unref (&((*mg_ctx)->interpolation_ctx->coarse_level_forests[num_lvl]));
    t8dg_dof_values_destroy (&((*mg_ctx)->interpolation_ctx->dofs_lvl[num_lvl]));
    t8dg_dof_values_destroy (&((*mg_ctx)->interpolation_ctx->dofs_lvl_derivation[num_lvl - 1]));
    ierr = PetscFree ((*mg_ctx)->interpolation_ctx->indexing_scheme[num_lvl]);
    CHKERRQ (ierr);
    ierr = MatDestroy (&((*mg_ctx)->Restriction_Mats[num_lvl - 1]));
    CHKERRQ (ierr);
    ierr = MatDestroy (&((*mg_ctx)->Prolongation_Mats[num_lvl - 1]));
    CHKERRQ (ierr);
    if (num_lvl != 1) {
      ierr = MatDestroy (&((*mg_ctx)->Smoothing_Mats[num_lvl - 1]));
      CHKERRQ (ierr);
    }
    ierr = KSPDestroy (&((*mg_ctx)->smoothers[num_lvl - 1]));
    CHKERRQ (ierr);
    t8dg_local_values_destroy (&((*mg_ctx)->interpolation_ctx->local_values_lvl[num_lvl]));
    t8dg_mortar_array_destroy (&((*mg_ctx)->interpolation_ctx->mortar_array_lvl[num_lvl]));
  }
  ierr = MatDestroy (&((*mg_ctx)->A_coarse_Mat));
  CHKERRQ (ierr);
  ierr = KSPDestroy (&((*mg_ctx)->coarse_solver));
  CHKERRQ (ierr);
  T8DG_FREE ((*mg_ctx)->coarse_matrix_ctx);
  T8DG_FREE ((*mg_ctx)->interpolation_ctx);
  return 0;
}

/******************/
/******************/
/* END Try multiple lvls */
/*****************/
/*****************/

/**********************************************************************************************/
/********************************* End of Multigrid-Routines **********************************/
/**********************************************************************************************/

/* A function that writes a t8dg_dof_values_t to a PETSc Vector */
PetscErrorCode
t8dg_precon_write_dof_to_vec (t8dg_dof_values_t * dofs, Vec * p_vec, PetscInt * indexing, size_t num_local_dofs)
{
  PetscErrorCode      ierr;
  /* Sets the values of a petsc vector to the entries of a t8dg_dof_values_t */
  ierr = VecSetValues (*p_vec, num_local_dofs, indexing, t8dg_dof_values_get_double_pointer (dofs, 0), INSERT_VALUES);
  CHKERRQ (ierr);
  /* Assemble the petsc vector */
  ierr = VecAssemblyBegin (*p_vec);
  CHKERRQ (ierr);
  ierr = VecAssemblyEnd (*p_vec);
  CHKERRQ (ierr);
  return 0;
}

/* A function that writes a PETSc Vector to a t8dg_dof_values_t */
PetscErrorCode
t8dg_precon_write_vec_to_dof (Vec * p_vec, t8dg_dof_values_t * dofs, size_t num_local_dofs)
{
  PetscErrorCode      ierr;
  const PetscScalar  *vec_reader;
  double             *dof_pointer;
  int                 dof_iter;

  /* Retrieve the local part of a petsc vector */
  ierr = VecGetArrayRead (*p_vec, &vec_reader);
  CHKERRQ (ierr);

  /* Get a double pointer to the data inside the sc_array_t of t8dg_dof_values_t */
  dof_pointer = t8dg_dof_values_get_double_pointer (dofs, 0);
  /* Overwrite the dof values with the entries of a petsc vector */
  /* it is assumed (not checked!) that the local petsc vector is of length dof_values->num_local_elements */
  for (dof_iter = 0; dof_iter < num_local_dofs; ++dof_iter) {
    dof_pointer[dof_iter] = (double) vec_reader[dof_iter];
  }
  /* Restore the array entries in the petsc vector */
  ierr = VecRestoreArrayRead (*p_vec, &vec_reader);
  CHKERRQ (ierr);

  return 0;
}

/* Thus function updates the preconditioner within the DIRK methods due to the varying coefficients of the different stages */
void
t8dg_precon_dirk_update_preconditioner (int selector, t8dg_precon_general_preconditioner_t * preconditioner,
                                        double stage_related_a_coefficient, double stage_current_time)
{
  switch (selector) {
  case 0:
    /* No preconditioning was selected, therefore, no updating is needed */
    break;
  case 1:
    /* Block-Jacobi-Preconditioner is selected */
    (preconditioner->jacobi_preconditioner_ctx->jacobi_ctx).current_implicit_coefficient = stage_related_a_coefficient;
    (preconditioner->jacobi_preconditioner_ctx->jacobi_ctx).current_point_in_time = stage_current_time;
    break;
  case 2:
    /* Currently, no preconditioner is assigned to '2' */
    break;
  case 3:
    /* Block-Gauss-Seidel-Preconditioner is selected */
    (preconditioner->jacobi_preconditioner_ctx->jacobi_ctx).current_implicit_coefficient = stage_related_a_coefficient;
    (preconditioner->jacobi_preconditioner_ctx->jacobi_ctx).current_point_in_time = stage_current_time;
    break;
  case 4:
    preconditioner->multiple_mg_lvls_ctx->mg_general_data->current_point_in_time = stage_current_time;
    preconditioner->multiple_mg_lvls_ctx->mg_general_data->current_coefficient = stage_related_a_coefficient;
    break;
  default:
    break;
  }
}

double
t8dg_precon_get_setup_time (t8dg_precon_general_preconditioner_t * preconditioner)
{
  return (preconditioner->preconditioner_setup_time);
}

#endif
