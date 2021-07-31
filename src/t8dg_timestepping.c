/*
 * timestepping.c
 *
 *  Created on: Mar 14, 2020
 *      Author: lukas
 */

#include<sc_containers.h>
#include "t8dg_dof.h"
#include "t8dg.h"
#include "t8dg_timestepping.h"
#include "t8dg_timestepping_dirk_coefficients.h"
#include "t8dg_values.h"
#include "t8dg_advect_diff_problem.h"
#include "t8dg_preconditioner.h"

#define TRY_MG_PRECONDITIONER 1

struct t8dg_timestepping_data
{
  int                 time_order;/**< time order of the Runge-kutta timestepping*/
  double              delta_t;  /**< time step */
  double              t;        /**< current time */
  double              T;        /**< end time */
  double              cfl;      /**< cfl number*/
  int                 step_number;
  int                 use_implicit_timestepping;
};

#if T8_WITH_PETSC
/* Matrix-free routine that mimics the application of the matrix resulting from the Implicit/backwards Euler_Method */
//extern PetscErrorCode MatMult_MF_Impl_Euler (Mat, Vec, Vec);
extern PetscErrorCode MatMult_MF_Impl_Euler (Mat, Vec, Vec);
#endif

#if T8_WITH_PETSC
/* Matrix-free routine that mimics the application of the matrices resultung (during the steps) of the second and third order DIRK methods */
extern PetscErrorCode MatMult_MF_DIRK (Mat, Vec, Vec);
#endif

#if T8_WITH_PETSC
/* Jacobi-Preconditioner for matrix-free operations */
extern PetscErrorCode JacobiShellPCCreate (t8dg_timestepping_precon_jacobi_ctx_t **, double, double);
extern PetscErrorCode JacobiShellPCSetUp (PC, Mat, Vec, t8dg_dof_values_t *, size_t);
extern PetscErrorCode JacobiShellPCApply (PC, Vec x, Vec y);
extern PetscErrorCode JacobiShellPCDestroy (PC);
#endif
/*pre computed butcher tableau values for rk with values only on first diagonal*/
double              rk1_b[1] = { 1 };

double              rk2_a[2] = { 1 };
double              rk2_b[2] = { 0.5, 0.5 };

double              rk3_a[3] = { 1. / 3, 2. / 3 };
double              rk3_b[3] = { 1. / 4, 0, 3. / 4 };

double              rk4_a[4] = { 0.5, 0.5, 1 };
double              rk4_b[4] = { 1. / 6, 1. / 3, 1. / 3, 1. / 6 };

void
t8dg_runge_kutta_fill_coefficients (int time_order, double **prk_a, double **prk_b, double **prk_c)
{
  switch (time_order) {
  case 1:
    *prk_a = NULL;
    *prk_b = rk1_b;
    break;
  case 2:
    *prk_a = rk2_a;
    *prk_b = rk2_b;
    break;
  case 3:
    *prk_a = rk3_a;
    *prk_b = rk3_b;
    break;
  case 4:
    *prk_a = rk4_a;
    *prk_b = rk4_b;
    break;
  default:
    break;
  }
  *prk_c = *prk_a;
}

void
t8dg_timestepping_runge_kutta_step (t8dg_time_matrix_application time_derivative,
                                    t8dg_timestepping_data_t * time_data, t8dg_dof_values_t ** pdof_array, void *user_data)
{
  int                 istep;
  double             *rk_a, *rk_b, *rk_c;
  t8dg_dof_values_t  *dof_beginning;
  t8dg_dof_values_t  *dof_change;
  t8dg_dof_values_t  *dof_new;
  t8dg_dof_values_t  *dof_step;
  double              time_beginning, time_current, time_step;
  int                 time_order = t8dg_timestepping_data_get_time_order (time_data);

  time_beginning = t8dg_timestepping_data_get_current_time (time_data);
  time_current = time_beginning;
  time_step = t8dg_timestepping_data_get_time_step (time_data);

  t8dg_runge_kutta_fill_coefficients (time_order, &rk_a, &rk_b, &rk_c);

  dof_beginning = t8dg_dof_values_clone (*pdof_array);
  dof_new = t8dg_dof_values_duplicate (dof_beginning);
  dof_step = t8dg_dof_values_duplicate (dof_beginning);
  dof_change = t8dg_dof_values_duplicate (dof_beginning);

  time_derivative (*pdof_array, dof_change, time_current, user_data);
  t8dg_dof_values_axpyz (rk_b[0] * time_step, dof_change, dof_beginning, dof_new);

  for (istep = 0; istep < time_order - 1; istep++) {
    /*calculate the y-value for which the derivative needs to be evaluated
     * since a has only values on the first minor diagonal, only the k from the step before and the original y is needed*/

    t8dg_dof_values_axpyz (rk_a[istep] * time_step, dof_change, dof_beginning, dof_step);
    t8dg_dof_values_swap (pdof_array, &dof_step);
    /* calculate the derivative at the step time and y value */

    time_current = time_beginning + rk_c[istep] * time_step;
    t8dg_timestepping_data_set_current_time (time_data, time_current);
    time_derivative (*pdof_array, dof_change, time_current, user_data);
    /*add weighted summand to result */
    t8dg_dof_values_axpy (rk_b[istep + 1] * time_step, dof_change, dof_new);
  }

  t8dg_timestepping_data_set_current_time (time_data, time_beginning + time_step);
  t8dg_dof_values_swap (pdof_array, &dof_new);

  t8dg_dof_values_destroy (&dof_beginning);
  t8dg_dof_values_destroy (&dof_change);
  t8dg_dof_values_destroy (&dof_new);
  t8dg_dof_values_destroy (&dof_step);
}

t8dg_timestepping_data_t *
t8dg_timestepping_data_new_cfl (int time_order, double start_time, double end_time, double cfl, int use_implicit_timestepping)
{
  T8DG_ASSERT (time_order > 0);
  t8dg_timestepping_data_t *time_data = T8DG_ALLOC (t8dg_timestepping_data_t, 1);
  time_data->time_order = time_order;
  time_data->t = start_time;
  time_data->T = end_time;
  time_data->cfl = cfl;
  time_data->step_number = 0;
  time_data->delta_t = -1;
  time_data->use_implicit_timestepping = use_implicit_timestepping;
  return time_data;
}

t8dg_timestepping_data_t *
t8dg_timestepping_data_new_constant_timestep (int time_order, double start_time, double end_time, double delta_t,
                                              int use_implicit_timestepping)
{
  T8DG_ASSERT (time_order > 0);
  t8dg_timestepping_data_t *time_data = T8DG_ALLOC (t8dg_timestepping_data_t, 1);
  time_data->time_order = time_order;
  time_data->t = start_time;
  time_data->T = end_time;
  time_data->cfl = 0;
  time_data->step_number = 0;
  time_data->delta_t = delta_t;
  time_data->use_implicit_timestepping = use_implicit_timestepping;
  return time_data;
}

void
t8dg_timestepping_data_destroy (t8dg_timestepping_data_t ** ptime_data)
{
  t8dg_timestepping_data_t *time_data = *ptime_data;
  time_data->time_order = -1;
  time_data->t = -1;
  time_data->T = -1;
  time_data->cfl = -1;
  time_data->delta_t = -1;
  time_data->use_implicit_timestepping = -1;
  T8_FREE (time_data);
  *ptime_data = NULL;
}

double
t8dg_timestepping_data_get_current_time (t8dg_timestepping_data_t * time_data)
{
  return time_data->t;
}

void
t8dg_timestepping_data_set_current_time (t8dg_timestepping_data_t * time_data, double current_time)
{
  time_data->t = current_time;
}

double
t8dg_timestepping_data_get_end_time (t8dg_timestepping_data_t * time_data)
{
  return time_data->T;
}

double
t8dg_timestepping_data_get_cfl (t8dg_timestepping_data_t * time_data)
{
  return time_data->cfl;
}

int
t8dg_timestepping_data_get_time_order (t8dg_timestepping_data_t * time_data)
{
  return time_data->time_order;
}

double
t8dg_timestepping_data_get_time_step (t8dg_timestepping_data_t * time_data)
{
  return time_data->delta_t;
}

void
t8dg_timestepping_data_set_time_step (t8dg_timestepping_data_t * time_data, double delta_t)
{
  if (delta_t < (time_data->T - time_data->t)) {
    time_data->delta_t = delta_t;
  }
  else {
    time_data->delta_t = time_data->T - time_data->t;
  }
}

int
t8dg_timestepping_data_is_endtime_reached (t8dg_timestepping_data_t * time_data)
{
  return !(time_data->t < time_data->T);
}

int
t8dg_timestepping_data_get_step_number (t8dg_timestepping_data_t * time_data)
{
  return time_data->step_number;
}

void
t8dg_timestepping_data_increase_step_number (t8dg_timestepping_data_t * time_data)
{
  time_data->step_number++;
}

double
t8dg_timestepping_data_get_time_left (t8dg_timestepping_data_t * time_data)
{
  return time_data->T - time_data->t;
}

/* Calculates the approximation of the next time step by using the implicit Euler Method in a matrix-free fashion */
void
t8dg_timestepping_implicit_euler (t8dg_time_matrix_application time_derivative,
                                  t8dg_timestepping_data_t * time_data, t8dg_dof_values_t ** pdof_array, void *user_data)
{
#if T8_WITH_PETSC
  size_t              iter;
  //size_t max_local_index;
  //double             *raw_local_dofs;
  double             *dof_values_ptr;
  int                 dof_iter;
  Mat                 A;
  Vec                 u, f;
  KSP                 ksp;
  PC                  pc;
  PetscScalar        *local_dofs;
  PetscScalar        *new_local_approx;
  PetscErrorCode      ierr;
  t8dg_timestepping_impl_euler_ctx_t appctx;
  t8dg_timestepping_precon_jacobi_ctx_t *pcshell_ctx;
  t8dg_dof_values_t  *dof_problem_new;
  PetscInt            its;
  PetscInt           *vec_global_index;
  t8_gloidx_t         global_offset_to_first_local_elem;

  /* Trial of a two-level multigrid preconditioner */
#if TRY_MG_PRECONDITIONER
  KSP                 coarse_solver;
  KSP                 smoother;
  PC                  coarse_pc;
  PC                  smoother_pc;
  Mat                 Restriction, Prolongation;
  t8dg_coarse_matrix_ctx_t cmat_ctx;
  Mat                 A_coarse;
  t8dg_mg_interpolating_ctx_t res_prol_ctx;
  t8dg_mg_coarse_lvl_t coarse_lvl;
#endif

  /* Store the current timestep */
  appctx.timestep = t8dg_timestepping_data_get_time_step (time_data);

  /* Store the next point in time which is going to be calculated */
  appctx.next_time_point = t8dg_timestepping_data_get_current_time (time_data) + (double) appctx.timestep;

  /* Store the time_derivative_function and the problem related data as well as the degrees of freedom */
  appctx.time_derivative_func = time_derivative;
  appctx.time_data = time_data;
  appctx.user_data = user_data;

  /* Store the degrees of freedom in the matrix-free application context (the current (original) degrees of freedom will be directly overwritten) */
  appctx.future_local_dofs = pdof_array;
  /* Is going to store the derivation of the degrees of freedom -> the evaluation of the time_derivative function */
  appctx.future_local_dofs_derivative = t8dg_dof_values_duplicate (*(appctx.future_local_dofs));

  /* the degrees of freedom (of the next time step) which will be filled with the solution */
  dof_problem_new = t8dg_dof_values_duplicate (*(appctx.future_local_dofs));

  /* Advance a time step */
  t8dg_timestepping_data_set_current_time (time_data, appctx.next_time_point);

  t8dg_debugf ("System wird aufgebaut und geloest\n");

        /********** Setting up and Solving the LS ****************/

  /* Calculate the process-local number of degrees of freedom */
  appctx.num_local_dofs =
    (size_t) (t8dg_dof_get_num_local_elements (*(appctx.future_local_dofs)) *
              t8dg_dof_get_max_num_element_dof (*(appctx.future_local_dofs)));

  /* Allocate space for an array which holds the global indices of the degrees of freedom */
  ierr = PetscMalloc1 (appctx.num_local_dofs, &vec_global_index);

  /* compute the offset of the first local element concerning the forest and the amout of degrees of freedom per element */
  global_offset_to_first_local_elem = t8_forest_get_first_local_element_id (t8dg_dof_values_get_forest (*(appctx.future_local_dofs)));
  if (global_offset_to_first_local_elem != 0) {
    global_offset_to_first_local_elem = global_offset_to_first_local_elem * t8dg_dof_get_max_num_element_dof (*(appctx.future_local_dofs));
  }

  /* Fill the array of global indices */
  for (iter = 0; iter < appctx.num_local_dofs; ++iter) {
    vec_global_index[iter] = iter + global_offset_to_first_local_elem;
  }

  /* Store a pointer to these global indices which will be needed for the assembly of the ouput-Vector during the matrix-free MatVec product */
  appctx.global_indexing = vec_global_index;

  /* Create the right-hand-side of the system resulting from the implicit Euler-Method */
  ierr = VecCreateMPI (PETSC_COMM_WORLD, (PetscInt) appctx.num_local_dofs, PETSC_DETERMINE, &f);
  CHKERRQ (ierr);
  /* Fill the vector with the current degrees of freedom resulting from the former time step */
  ierr = VecSetValues (f, appctx.num_local_dofs, vec_global_index, t8dg_dof_get_double_pointer_to_array (*pdof_array), INSERT_VALUES);
  CHKERRQ (ierr);
  ierr = VecAssemblyBegin (f);
  CHKERRQ (ierr);
  ierr = VecAssemblyEnd (f);
  CHKERRQ (ierr);
  ierr = PetscObjectSetName ((PetscObject) f, "Right Hand Side");
  CHKERRQ (ierr);

  /* Create an Approximation Vector storing the solution of the LS */
  /* Use the size and allocation similar to the right-hand-side (-> f) */
  ierr = VecDuplicate (f, &u);
  CHKERRQ (ierr);
  ierr = VecCopy (f, u);
  CHKERRQ (ierr);
  ierr = PetscObjectSetName ((PetscObject) u, "Approximation");
  CHKERRQ (ierr);

  /* Setting up matrix-free Matrix */
  /* Create a matrix shell with local dimensions equal to the dimension of the Vec containing the process-local degrees of freedom and add an application context needed by the matrix-free MatVec multiplication */
  ierr =
    MatCreateShell (PETSC_COMM_WORLD, (PetscInt) appctx.num_local_dofs, (PetscInt) appctx.num_local_dofs, PETSC_DETERMINE, PETSC_DETERMINE,
                    (void *) &appctx, &A);
  CHKERRQ (ierr);
  ierr = PetscObjectSetName ((PetscObject) A, "Matrix-Free Application-Matrix of Imp Euler");
  CHKERRQ (ierr);

  /* Define the (multiplicative) MatVec-Operation */
  ierr = MatShellSetOperation (A, MATOP_MULT, (void (*)(void)) MatMult_MF_Impl_Euler);
  CHKERRQ (ierr);

  /* Setting up KSP Solver */
  KSPCreate (PETSC_COMM_WORLD, &ksp);
  /* Set KSP: Set Operators; here: matrix that defines LS also serve as the matrix that definies the preconditioner */
  ierr = KSPSetOperators (ksp, A, A);
  CHKERRQ (ierr);
  /* Extract KSP and PC Context from  KSP Context */
  ierr = KSPGetPC (ksp, &pc);
  CHKERRQ (ierr);
  /* Select GMRES as Solver; default: Restarted GMRES(30) */
  /* If GMRES is used as smoother -> KSPFGMRES muste be used */
  ierr = KSPSetType (ksp, KSPFGMRES);
  CHKERRQ (ierr);
  /* Set error tolerances in the (outer) GMRES iteration */
  ierr = KSPSetTolerances (ksp, 1.0e-9, 1.0e-9, PETSC_DEFAULT, PETSC_DEFAULT);
  CHKERRQ (ierr);
  t8dg_debugf ("KSPFGMRES was selected\n");
#if 0
  /* Select a Preconditioner - without Preconditioning */
  ierr = PCSetType (pc, PCNONE);
  CHKERRQ (ierr);
#endif
#if 0
  /* Select Jacobi Preconditioner */
  ierr = PCSetType (pc, PCSHELL);
  CHKERRQ (ierr);
  ierr = JacobiShellPCCreate (&pcshell_ctx, appctx.timestep, 1.0);
  CHKERRQ (ierr);
  ierr = PCShellSetApply (pc, JacobiShellPCApply);
  CHKERRQ (ierr);
  ierr = PCShellSetContext (pc, pcshell_ctx);
  CHKERRQ (ierr);
  ierr = PCShellSetDestroy (pc, JacobiShellPCDestroy);
  CHKERRQ (ierr);
  ierr = PCShellSetName (pc, "Try Jacobi Preconditioner");
  CHKERRQ (ierr);
  t8dg_debugf ("Anzahl totaler elemente: %d\n", (int) t8dg_dof_get_num_total_elements (appctx.future_local_dofs));
  ierr = JacobiShellPCSetUp (pc, A, u, appctx.future_local_dofs, (size_t) t8dg_dof_get_num_total_elements (appctx.future_local_dofs));
  CHKERRQ (ierr);
#endif
#if TRY_MG_PRECONDITIONER
  t8dg_debugf ("multigrid_preconditioning begin\n");
  t8dg_mg_set_up_two_lvl_precon ((t8dg_linear_advection_diffusion_problem_t *) user_data, appctx.future_local_dofs, time_derivative,
                                 &coarse_lvl, t8dg_adapt_multigrid_coarsen_finest_level, &res_prol_ctx, &cmat_ctx, vec_global_index);
  t8dg_debugf ("set up complete\n");
  t8dg_mg_create_coarse_lvl_matrix (&A_coarse, &coarse_lvl, &cmat_ctx);
  t8dg_debugf ("coarse matrix complete\n");
  t8dg_mg_create_restriction_matrix (&Restriction, &res_prol_ctx);
  t8dg_debugf ("restriction matrix complete\n");
  t8dg_mg_create_prolongation_matrix (&Prolongation, &res_prol_ctx);
  t8dg_debugf ("prolongation matrix complete\n");
  /* Set Coarse Solver, Smoother, Matrices and Operators */
  ierr = PCSetType (pc, PCMG);
  CHKERRQ (ierr);
  ierr = PCMGSetLevels (pc, 2, NULL);
  CHKERRQ (ierr);
  ierr = PCMGSetType (pc, PC_MG_MULTIPLICATIVE);
  CHKERRQ (ierr);
  ierr = PCMGSetCycleType (pc, PC_MG_CYCLE_V);
  CHKERRQ (ierr);
  ierr = PCMGGetCoarseSolve (pc, &coarse_solver);
  CHKERRQ (ierr);
  ierr = KSPSetOperators (coarse_solver, A_coarse, A_coarse);
  CHKERRQ (ierr);
  ierr = KSPSetTolerances (coarse_solver, 1.0e-8, 1.0e-8, PETSC_DEFAULT, PETSC_DEFAULT);
  CHKERRQ (ierr);
  ierr = KSPSetType (coarse_solver, KSPGMRES);
  CHKERRQ (ierr);
  ierr = KSPGetPC (coarse_solver, &coarse_pc);
  CHKERRQ (ierr);
  ierr = PCSetType (coarse_pc, PCNONE);
  CHKERRQ (ierr);
  ierr = PCMGGetSmoother (pc, 1, &smoother);
  CHKERRQ (ierr);
  ierr = KSPSetOperators (smoother, A, A);
  CHKERRQ (ierr);
  ierr = KSPSetType (smoother, KSPGMRES);
  CHKERRQ (ierr);
  ierr = KSPGetPC (smoother, &smoother_pc);
  CHKERRQ (ierr);
  ierr = PCSetType (smoother_pc, PCNONE);
  CHKERRQ (ierr);
  /* Amount of pre- and post-smoothing steps */
  ierr = PCMGSetNumberSmooth (pc, 3);
  CHKERRQ (ierr);
  /* Set Restriction operator */
  ierr = PCMGSetRestriction (pc, 1, Restriction);
  CHKERRQ (ierr);
  /* Set Prolongation operator */
  ierr = PCMGSetInterpolation (pc, 1, Prolongation);
  CHKERRQ (ierr);
  /* Set residual calculation - my be default since only Mat's and Vec's are used */
  ierr = PCMGSetResidual (pc, 1, PCMGResidualDefault, A);
  CHKERRQ (ierr);
  /* Provide work space vectors */
  /* No they should be automatically managed by PETSc */
  t8dg_debugf ("all multigrid operators were set\n");
#endif

  /* Set initial guess to non-zero (to the last solution) otherwise it will filled with zeros */
  //ierr = KSPSetInitialGuessNonzero (ksp, PETSC_TRUE);
  //CHKERRQ (ierr);
  /* Use eventually further specified properties via the command line, otherwise it is just going to be set up */
  ierr = KSPSetFromOptions (ksp);
  CHKERRQ (ierr);

  t8dg_debugf ("\nGMRES got called\n");

  /* Solve the Linear System */
  ierr = KSPSolve (ksp, f, u);
  CHKERRQ (ierr);

  t8dg_debugf ("\nGMRES solve completed\n");

  /* View whether or not the GMRES did converge */
  ierr = KSPConvergedRateView (ksp, PETSC_VIEWER_STDOUT_WORLD);
  CHKERRQ (ierr);

        /******************* The LS has been solved *********************/

  /* Retrieve the local part of the 'solution' Vector u for reading purposes */
  ierr = VecGetArrayRead (u, &new_local_approx);
  CHKERRQ (ierr);

  /* Get a double pointer to the data inside the sc_array_t of t8dg_dof_values_t */
  dof_values_ptr = t8dg_dof_get_double_pointer_to_array (dof_problem_new);

  /* Overwrite the dof_problem_new with the newly calculated degrees of freedom */
  for (dof_iter = 0; dof_iter < appctx.num_local_dofs; ++dof_iter) {
    dof_values_ptr[dof_iter] = (double) new_local_approx[dof_iter];
  }
  /* Restore the array entries */
  ierr = VecRestoreArrayRead (u, &new_local_approx);
  CHKERRQ (ierr);

  t8dg_debugf ("\nCalculated solution (dofs) of the next timestep:\n");
  t8dg_dof_values_debug_print (dof_problem_new);

  /* Swap the degrees of freedom associated with the problem with the newly calculated degrees of freedom */
  t8dg_dof_values_swap (pdof_array, &dof_problem_new);

  //t8dg_debugf("Dfference between A*u_{k+1} und initial dofs/rhs (u_{k})\n");
  //ierr = VecAXPY(u, -1.0, f); CHKERRQ(ierr);
  //VecView(u, PETSC_VIEWER_STDOUT_WORLD);
  t8dg_debugf ("here memory corruption \n");
  /* Free used resources */
  t8dg_dof_values_destroy (&(appctx.future_local_dofs_derivative));
  t8dg_dof_values_destroy (&dof_problem_new);
  t8dg_debugf ("vor forest unref\n");
#if TRY_MG_PRECONDITIONER
  /* Free multigrid preconditioner data */
  t8_forest_unref (&coarse_lvl.forest_coarsened);
  t8dg_values_destroy_adapt_data (coarse_lvl.dg_values, &coarse_lvl.tmp_mortar_coarse_lvl);
  t8dg_debugf ("hier\n");
  t8dg_dof_values_destroy ((coarse_lvl.dof_values_adapt));
  t8dg_dof_values_destroy (&(cmat_ctx.problem_dofs_derivation));
  t8dg_dof_values_destroy (&((coarse_lvl.adapt_data)->dof_values));
  t8dg_debugf ("oder hier\n");
  t8dg_dof_values_destroy (&((coarse_lvl.adapt_data)->dof_values_adapt));
  t8dg_debugf ("oder hier?\n");
  ierr = MatDestroy (&A_coarse);
  CHKERRQ (ierr);
  t8dg_debugf ("destroyed A\n");
  ierr = MatDestroy (&Restriction);
  CHKERRQ (ierr);
  t8dg_debugf ("destroyed Restriction\n");
  ierr = MatDestroy (&Prolongation);
  CHKERRQ (ierr);
  t8dg_debugf ("destroyed Prolongation\n");
  ierr = KSPDestroy (&coarse_solver);
  CHKERRQ (ierr);
  t8dg_debugf ("destroyed coarse_solver\n");
  ierr = KSPDestroy (&smoother);
  CHKERRQ (ierr);
  t8dg_debugf ("destroyed smoother\n");
  ierr = PetscFree (coarse_lvl.global_indexing);
  CHKERRQ (ierr);
  t8dg_debugf ("destroyed global indexing\n");
#endif
  /* Free PETSc resources */
  ierr = MatDestroy (&A);
  CHKERRQ (ierr);
  ierr = VecDestroy (&f);
  CHKERRQ (ierr);
  ierr = VecDestroy (&u);
  CHKERRQ (ierr);
  //ierr = KSPDestroy (&ksp);
  //CHKERRQ (ierr);
  ierr = PetscFree (vec_global_index);
  CHKERRQ (ierr);

  t8dg_debugf ("\nImplicit Euler-Method has been completed\n");
#else
  t8dg_global_productionf
    ("t8code/t8dg is currently not configured with PETSc.\n\nIn order to use this function t8dg needs to be configured with PETSc.\n");
#endif
}

#if T8_WITH_PETSC
/* Function that mimics the multiplication of the system matrix resulting from the implicit Euler_Method */
PetscErrorCode
MatMult_MF_Impl_Euler (Mat A, Vec in, Vec out)
{
  t8dg_debugf ("in matmultimpleuler\n");
  t8dg_timestepping_impl_euler_ctx_t *appctx;
  PetscErrorCode      ierr;
  PetscScalar        *current_local_approx;
  int                 dof_iter;
  double             *dof_values_ptr;

  /* Get the matrix-free application context */
  ierr = MatShellGetContext (A, &appctx);
  CHKERRQ (ierr);
#if 0
  /* Retrieve the local part of the 'in' Vector for reading purposes */
  ierr = VecGetArrayRead (in, &current_local_approx);
  CHKERRQ (ierr);

  /* Get a pointer to the original degrees of freedom of the probblem */
  dof_values_ptr = t8dg_dof_get_double_pointer_to_array (appctx->future_local_dofs);

  /* Overwrite them with the 'in' Vector of the Matrix Multiplication */
  /* Because this MatMult is only used within the GMRES, it is okay to overwrite the original problem dofs without resetting them */
  for (dof_iter = 0; dof_iter < appctx->num_local_dofs; ++dof_iter) {
    dof_values_ptr[dof_iter] = current_local_approx[dof_iter];
  }

  /* Restore the array entries of the 'in' Vector */
  ierr = VecRestoreArrayRead (in, &current_local_approx);
  CHKERRQ (ierr);
#endif
  t8dg_precon_write_vec_to_dof (&in, *(appctx->future_local_dofs));

  /* Calculate the time derivative (concerning the next time step) of the degrees of freedom passed to the MatMult (= 'in' Vector) and store their derivation inside future_local_dofs_derivative */
  (appctx->time_derivative_func) (*(appctx->future_local_dofs), appctx->future_local_dofs_derivative, appctx->next_time_point,
                                  appctx->user_data);

  /* calculate the sum of the initially passed degrees of freedom ('in' Vector) and their (scaled) derivation */
  t8dg_dof_values_axpy (-(appctx->timestep), appctx->future_local_dofs_derivative, *(appctx->future_local_dofs));

  t8dg_precon_write_dof_to_vec (*(appctx->future_local_dofs), &out, appctx->global_indexing);
#if 0
  /* Write the result of the matrix-application to the 'in' Vector in the 'out' Vector */
  ierr =
    VecSetValues (out, appctx->num_local_dofs, appctx->global_indexing, t8dg_dof_get_double_pointer_to_array (appctx->future_local_dofs),
                  INSERT_VALUES);
  CHKERRQ (ierr);
  ierr = VecAssemblyBegin (out);
  CHKERRQ (ierr);
  ierr = VecAssemblyEnd (out);
  CHKERRQ (ierr);
#endif
  return 0;
}
#endif

/* Calculates the approximation of the next time step using a DIRK method (either DIRK(2,2) or DIRK(3,3) can be selected) */
void
t8dg_timestepping_dirk (t8dg_time_matrix_application time_derivative,
                        t8dg_timestepping_data_t * time_data, t8dg_dof_values_t ** pdof_array, void *user_data, int num_order_stages)
{
#if T8_WITH_PETSC
  T8DG_ASSERT ((num_order_stages == 2 || num_order_stages == 3));

  /* Declare the general DIRK coefficients */
  double              (*dirk_a_coeff)[num_order_stages][num_order_stages];
  double              (*dirk_b_coeff)[num_order_stages];
  double              (*dirk_c_coeff)[num_order_stages];

  int                 step_iter;
  int                 iter;
  PetscErrorCode      ierr;
  t8dg_timestepping_dirk_ctx_t appctx;
  t8dg_timestepping_precon_jacobi_ctx_t *pcshell_ctx;
  /* This array holds the pointer to solution calculated of the different steps */
  t8dg_dof_values_t  *local_dofs_step[num_order_stages];
  double              current_time = t8dg_timestepping_data_get_current_time (time_data);
  PetscInt           *vec_global_index;
  t8_gloidx_t         global_offset_to_first_local_elem;
  Mat                 A;
  Vec                 u, f;
  KSP                 ksp;
  PC                  pc;
  double             *dof_values_ptr;
  PetscScalar        *new_local_approx;

  /* Store the time_derivative_function and the problem related data */
  appctx.time_derivative_func = time_derivative;
  appctx.time_data = time_data;
  appctx.user_data = user_data;
  appctx.timestep = t8dg_timestepping_data_get_time_step (time_data);

  /* Just for preconditioning */
  appctx.next_time_point = current_time;

  /* Clone the current degrees of freedom */
  appctx.local_derivation_degrees = *pdof_array;
  appctx.current_local_dofs = t8dg_dof_values_clone (appctx.local_derivation_degrees);

  /* Duplicate the other t8dg_dof_values_t */
  appctx.future_local_dofs_derivative = t8dg_dof_values_duplicate (appctx.local_derivation_degrees);
  appctx.future_local_dofs_step = t8dg_dof_values_duplicate (appctx.local_derivation_degrees);

  /* Select whether the order 2 or order 3 DIRK method has been choosen */
  if (num_order_stages == 2) {
    dirk_a_coeff = &dirk22_a_coeff;
    dirk_b_coeff = &dirk22_b_coeff;
    dirk_c_coeff = &dirk22_c_coeff;
  }
  else if (num_order_stages == 3) {
    dirk_a_coeff = &dirk33_a_coeff;
    dirk_b_coeff = &dirk33_b_coeff;
    dirk_c_coeff = &dirk33_c_coeff;
  }

  /* Calculate the process-local number of degrees of freedom */
  appctx.num_local_dofs =
    (size_t) (t8dg_dof_get_num_local_elements (appctx.local_derivation_degrees) *
              t8dg_dof_get_max_num_element_dof (appctx.local_derivation_degrees));

  /* Allocate space for an array which holds the global indices of the degrees of freedom */
  ierr = PetscMalloc1 (appctx.num_local_dofs, &vec_global_index);

  /* compute the offset of the first local element concerning the forest and the amout of degrees of freedom per element */
  global_offset_to_first_local_elem = t8_forest_get_first_local_element_id (t8dg_dof_values_get_forest (appctx.local_derivation_degrees));
  if (global_offset_to_first_local_elem != 0) {
    global_offset_to_first_local_elem =
      global_offset_to_first_local_elem * t8dg_dof_get_max_num_element_dof (appctx.local_derivation_degrees);
  }

  /* Fill the process-local array of global indices */
  for (iter = 0; iter < appctx.num_local_dofs; ++iter) {
    vec_global_index[iter] = iter + global_offset_to_first_local_elem;
  }

  /* Store a pointer to these global indices which will be needed for the assembly of the ouput-Vector during the matrix-free MatVec product */
  appctx.global_indexing = vec_global_index;

  /*********************** Setting up the Vectors and Matrix ***************************/

  /* Create the right-hand-side of the system resulting from the implicit Euler-Method */
  ierr = VecCreateMPI (PETSC_COMM_WORLD, (PetscInt) appctx.num_local_dofs, PETSC_DETERMINE, &f);
  CHKERRQ (ierr);
  /* Fill the vector with the current degrees of freedom resulting from the former time step */
  ierr = VecSetValues (f, appctx.num_local_dofs, vec_global_index, t8dg_dof_get_double_pointer_to_array (*pdof_array), INSERT_VALUES);
  CHKERRQ (ierr);
  ierr = VecAssemblyBegin (f);
  CHKERRQ (ierr);
  ierr = VecAssemblyEnd (f);
  CHKERRQ (ierr);
  ierr = PetscObjectSetName ((PetscObject) f, "Right Hand Side");
  CHKERRQ (ierr);

  /* Create an Approximation Vector storing the solution of the LS */
  /* Use the size and allocation similar to the right-hand-side (-> f) */
  ierr = VecDuplicate (f, &u);
  CHKERRQ (ierr);
  /* Copy the entries of f to u -> u is used as an initial guess; otherwise u will be overwitten with zeros during KSPSolve (if non_zero_guess is false) */
  ierr = VecCopy (f, u);
  CHKERRQ (ierr);
  ierr = VecAssemblyBegin (u);
  CHKERRQ (ierr);
  ierr = VecAssemblyEnd (u);
  CHKERRQ (ierr);
  ierr = PetscObjectSetName ((PetscObject) u, "Approximation");
  CHKERRQ (ierr);

  /* Setting up matrix-free Matrix */
  /* Create a matrix shell with local dimensions equal to the dimension of the Vec containing the process-local degrees of freedom and add an application context needed by the matrix-free MatVec multiplication */
  ierr =
    MatCreateShell (PETSC_COMM_WORLD, (PetscInt) appctx.num_local_dofs, (PetscInt) appctx.num_local_dofs, PETSC_DETERMINE, PETSC_DETERMINE,
                    (void *) &appctx, &A);
  CHKERRQ (ierr);
  ierr = PetscObjectSetName ((PetscObject) A, "Matrix-Free Application-Matrix of DIRK methods");
  CHKERRQ (ierr);

  /* Define the (multiplicative) MatVec-Operation */
  ierr = MatShellSetOperation (A, MATOP_MULT, (void (*)(void)) MatMult_MF_DIRK);
  CHKERRQ (ierr);

  /* Setting up KSP Solver */
  KSPCreate (PETSC_COMM_WORLD, &ksp);
  /* Set KSP: Set Operators; here: matrix that defines LS also serve as the matrix that definies the preconditioner */
  ierr = KSPSetOperators (ksp, A, A);
  CHKERRQ (ierr);
  /* Extract KSP and PC Context from  KSP Context */
  ierr = KSPGetPC (ksp, &pc);
  CHKERRQ (ierr);
  /* Select GMRES as Solver; default: Restarted GMRES(30) */
  ierr = KSPSetType (ksp, KSPGMRES);
  CHKERRQ (ierr);
  /* Select a Preconditioner - without Preconditioning */
  //ierr = PCSetType (pc, PCNONE);
  //CHKERRQ (ierr);
  /* Select Jacobi Preconditioner */
  ierr = PCSetType (pc, PCSHELL);
  CHKERRQ (ierr);
  ierr = JacobiShellPCCreate (&pcshell_ctx, appctx.timestep, (*dirk_b_coeff)[0]);
  CHKERRQ (ierr);
  ierr = PCShellSetApply (pc, JacobiShellPCApply);
  CHKERRQ (ierr);
  ierr = PCShellSetContext (pc, pcshell_ctx);
  CHKERRQ (ierr);
  ierr = PCShellSetDestroy (pc, JacobiShellPCDestroy);
  CHKERRQ (ierr);
  ierr = PCShellSetName (pc, "Try Jacobi Preconditioner");
  CHKERRQ (ierr);
  appctx.next_time_point = current_time + (((*dirk_c_coeff)[step_iter]) * appctx.timestep);
  t8dg_timestepping_data_set_current_time (time_data, appctx.next_time_point);
  ierr =
    JacobiShellPCSetUp (pc, A, u, appctx.local_derivation_degrees,
                        (size_t) t8dg_dof_get_num_total_elements (appctx.local_derivation_degrees));
  CHKERRQ (ierr);
  appctx.next_time_point = current_time;
  t8dg_timestepping_data_set_current_time (time_data, appctx.next_time_point);
  ierr = KSPSetTolerances (ksp, 1.0e-8, 1.0e-8, PETSC_DEFAULT, PETSC_DEFAULT);
  CHKERRQ (ierr);
  /* Set initial guess to non-zero (to the last solution) otherwise it will filled with zeros */
  ierr = KSPSetInitialGuessNonzero (ksp, PETSC_TRUE);
  CHKERRQ (ierr);
  /* Use eventually further specified properties via the command line, otherwise it is just going to be set up */
  ierr = KSPSetFromOptions (ksp);
  CHKERRQ (ierr);

  /************************* Creation of Vectors and Matrix finished *************************/

  t8dg_debugf ("\nDIRK Method of order %d has been called\n", num_order_stages);

  /* In every RKV step a LS has to be solved in order to obtain an intermediate solution which is than further processed */
  for (step_iter = 0; step_iter < num_order_stages; ++step_iter) {

    /* Set the next time point of the upcoming RKV step */
    appctx.next_time_point = current_time + (((*dirk_c_coeff)[step_iter]) * appctx.timestep);
    t8dg_timestepping_data_set_current_time (time_data, appctx.next_time_point);

    /* Set the current coefficient needed in the matrix-free application for the next step */
    appctx.dirk_current_a_coeff = (*dirk_a_coeff)[step_iter][step_iter];

    /* Fill the future_local_dofs_step with zeros, so that the axpy product will not be irritated */
    t8dg_dof_values_set_zero (appctx.future_local_dofs_step);

    /* Accumulate the (weighted) derivations of the degrees of freedom of the former stages */
    for (iter = 0; iter < step_iter; ++iter) {
      t8dg_dof_values_axpy ((appctx.timestep * ((*dirk_a_coeff)[step_iter][iter])), local_dofs_step[iter], appctx.future_local_dofs_step);
    }
    /* Add the (weighted) derivations of the degrees of freedom of the former stages to the right-hand side */
    if (step_iter > 0) {
      ierr =
        VecSetValues (f, appctx.num_local_dofs, vec_global_index, t8dg_dof_get_double_pointer_to_array (appctx.future_local_dofs_step),
                      ADD_VALUES);
      CHKERRQ (ierr);
      ierr = VecAssemblyBegin (f);
      CHKERRQ (ierr);
      ierr = VecAssemblyEnd (f);
    }

    /* store the current b coefficient of the dirk method in the preconditoning context */
    //pcshell_ctx->current_b_coeff = (*dirk_b_coeff)[step_iter];
    /* In a second or third step the preconditioner gets updated */
    /*
       if (step_iter > 0)
       {
       appctx.local_derivation_degrees = t8dg_dof_values_clone(appctx.future_local_dofs_step);
       }
     */
    /* Solve the Linear System */
    ierr = KSPSolve (ksp, f, u);
    CHKERRQ (ierr);

    /* View whether or not the GMRES did converge */
    ierr = KSPConvergedRateView (ksp, PETSC_VIEWER_STDOUT_WORLD);
    CHKERRQ (ierr);

    /* Write the solution to future_local_dofs_steps */
    /* Retrieve the local part of the 'solution' Vector u for reading purposes */
    ierr = VecGetArrayRead (u, &new_local_approx);
    CHKERRQ (ierr);

    /* Get a double pointer to the data inside the sc_array_t of t8dg_dof_values_t */
    dof_values_ptr = t8dg_dof_get_double_pointer_to_array (appctx.local_derivation_degrees);

    /* Overwrite the dof_problem_new with the newly calculated degrees of freedom */
    for (iter = 0; iter < appctx.num_local_dofs; ++iter) {
      dof_values_ptr[iter] = (double) new_local_approx[iter];
    }

    /* Restore the array entries */
    ierr = VecRestoreArrayRead (u, &new_local_approx);
    CHKERRQ (ierr);

    /* Create a new empty t8dg_dof_values_t and swap with the current solution */
    local_dofs_step[step_iter] = t8dg_dof_values_duplicate (appctx.local_derivation_degrees);

    /* Store the derivation of degrees of freedom calculated in the current step in local_dofs_step */
    (appctx.time_derivative_func) (appctx.local_derivation_degrees, local_dofs_step[step_iter], appctx.next_time_point, appctx.user_data);

  }

  /* Add the values of the staged to the former approximation to obtain the degrees of freedom of the next step */
  for (iter = 0; iter < num_order_stages; ++iter) {
    t8dg_dof_values_axpy ((appctx.timestep * (*dirk_b_coeff)[iter]), local_dofs_step[iter], appctx.current_local_dofs);
  }

  /* Swap the former degrees of freedom with the newly calculated degrees of freedom */
  t8dg_dof_values_swap (pdof_array, &(appctx.current_local_dofs));

  t8dg_debugf ("\nDIRK Method completed\n");

  /* Free used resources */
  t8dg_dof_values_destroy (&(appctx.future_local_dofs_derivative));
  t8dg_dof_values_destroy (&(appctx.future_local_dofs_step));
  t8dg_dof_values_destroy (&(appctx.current_local_dofs));
  for (iter = 0; iter < num_order_stages; ++iter) {
    t8dg_dof_values_destroy (&local_dofs_step[iter]);
  }

  /* Free PETSc resources */
  ierr = MatDestroy (&A);
  CHKERRQ (ierr);
  ierr = VecDestroy (&f);
  CHKERRQ (ierr);
  ierr = VecDestroy (&u);
  CHKERRQ (ierr);
  ierr = KSPDestroy (&ksp);
  CHKERRQ (ierr);
  ierr = PetscFree (vec_global_index);
  CHKERRQ (ierr);

#else
  t8dg_global_productionf
    ("t8code/t8dg is currently not configured with PETSc.\n\nIn order to use this function t8dg needs to be configured with PETSc.\n");
#endif
}

#if T8_WITH_PETSC
PetscErrorCode
MatMult_MF_DIRK (Mat A, Vec in, Vec out)
{
  t8dg_timestepping_dirk_ctx_t *appctx;
  PetscErrorCode      ierr;
  PetscScalar        *current_local_approx;
  int                 dof_iter;
  double             *dof_values_ptr;

  /* Get the mtarix-free application context */
  ierr = MatShellGetContext (A, &appctx);
  CHKERRQ (ierr);

  /* Retrieve the local part of the 'in' Vector for reading purposes */
  ierr = VecGetArrayRead (in, &current_local_approx);
  CHKERRQ (ierr);

  /* Get a pointer to the original degrees of freedom of the probblem */
  dof_values_ptr = t8dg_dof_get_double_pointer_to_array (appctx->local_derivation_degrees);

  /* Overwrite them with the 'in' Vector of the Matrix Multiplication */
  for (dof_iter = 0; dof_iter < appctx->num_local_dofs; ++dof_iter) {
    dof_values_ptr[dof_iter] = current_local_approx[dof_iter];
  }

  /* Restore the array entries of the 'in' Vector */
  ierr = VecRestoreArrayRead (in, &current_local_approx);
  CHKERRQ (ierr);

  /* Get the derivation of the degrees of freedom of the current RKV step and store them in appctx->future_local_dofs_derivative */
  (appctx->time_derivative_func) (appctx->local_derivation_degrees, appctx->future_local_dofs_derivative, appctx->next_time_point,
                                  appctx->user_data);

  /* Subtract the derivation of the degrees from the dofs of the 'in' Vector */
  t8dg_dof_values_axpy (-(appctx->timestep * appctx->dirk_current_a_coeff), appctx->future_local_dofs_derivative,
                        appctx->local_derivation_degrees);

  /* Subtract the (weighted) derivations of the degrees of freedom calculated in the former stages which are stored in future_local_dofs_step */
  //t8dg_dof_values_axpy (1.0, appctx->future_local_dofs_step, appctx->local_derivation_degrees);

  /* Write the result of the matrix-application to the 'in' Vector in the 'out' Vector */
  ierr =
    VecSetValues (out, appctx->num_local_dofs, appctx->global_indexing,
                  t8dg_dof_get_double_pointer_to_array (appctx->local_derivation_degrees), INSERT_VALUES);
  CHKERRQ (ierr);

  /* Assemble the output vector */
  ierr = VecAssemblyBegin (out);
  CHKERRQ (ierr);
  ierr = VecAssemblyEnd (out);
  CHKERRQ (ierr);

  return 0;
}
#endif

/* Select either an implicit or an explicit Runge Kutta method */
void
t8dg_timestepping_choose_impl_expl_method (t8dg_time_matrix_application time_derivative,
                                           t8dg_timestepping_data_t * time_data, t8dg_dof_values_t ** pdof_array, void *user_data)
{
  if (time_data->use_implicit_timestepping == 0) {
    /* An explicit time stepping method has been chosen */
    t8dg_debugf ("Explicit RKV of order %d has been called.\n", time_data->time_order);
    t8dg_timestepping_runge_kutta_step (time_derivative, time_data, pdof_array, user_data);
  }
  else if (time_data->use_implicit_timestepping == 1) {
    /* An implicit time stepping method has been chosen */
    switch (time_data->time_order) {
    case 1:
      t8dg_debugf ("The implicit Euler method has been called.\n");
      t8dg_timestepping_implicit_euler (time_derivative, time_data, pdof_array, user_data);
      break;
    case 2:
      t8dg_debugf ("The DIRK(2,2) method has been called.\n");
      t8dg_timestepping_dirk (time_derivative, time_data, pdof_array, user_data, 2);
      break;
    case 3:
      t8dg_debugf ("The DIRK(3,3) method has been called.\n");
      t8dg_timestepping_dirk (time_derivative, time_data, pdof_array, user_data, 3);
      break;
    default:
      t8dg_debugf ("An implicit DIRK method of order %d has been called, which is not implemented.\n", time_data->time_order);
    }
  }
}

/************** Preconditioners **************/
/* Jacobi Preconditioning routines */
#if T8_WITH_PETSC
/* Creates a Jacobi preconditioner with corresponding 'preconditioning context' */
PetscErrorCode
JacobiShellPCCreate (t8dg_timestepping_precon_jacobi_ctx_t ** pcshell, double timestep, double b_coeff)
{
  PetscErrorCode      ierr;
  t8dg_timestepping_precon_jacobi_ctx_t *pcctx;

  /* Create a new jacobi preconditioning context */
  ierr = PetscNew (&pcctx);
  CHKERRQ (ierr);
  /* Set Vec that represents the matrix diagonal */
  pcctx->matrix_diagonal = 0;
  pcctx->timestep = timestep;
  pcctx->current_b_coeff = b_coeff;
  /* return the preconditioning context */
  *pcshell = pcctx;

  return 0;
}

/* Sets up the Jacobi preconditioner - calculates main diagonal of the system matrix */
PetscErrorCode
JacobiShellPCSetUp (PC pc, Mat pmat, Vec u, t8dg_dof_values_t * current_dofs, size_t num_total_elements)
{
  PetscErrorCode      ierr;

  t8dg_timestepping_precon_jacobi_ctx_t *pcshell;
  Vec                 mat_diagonal;
  t8dg_timestepping_impl_euler_ctx_t *appctx;
  t8dg_dof_values_t  *identity_vector;
  double             *dof_pointer_identity_vec;

  ierr = MatShellGetContext (pmat, &appctx);
  CHKERRQ (ierr);

  identity_vector = t8dg_dof_values_duplicate (current_dofs);

  ierr = VecDuplicate (u, &mat_diagonal);
  CHKERRQ (ierr);

  ierr = PCShellGetContext (pc, (void **) &pcshell);
  CHKERRQ (ierr);

  t8dg_linear_advection_diffusion_problem_t *problem = (t8dg_linear_advection_diffusion_problem_t *) appctx->user_data;

#if 1
  t8dg_advect_diff_problem_jacobi_precon (problem, current_dofs, identity_vector, pcshell->timestep, appctx->num_local_dofs,
                                          appctx->num_local_dofs);

  dof_pointer_identity_vec = t8dg_dof_get_double_pointer_to_array (identity_vector);

  ierr =
    VecSetValues (mat_diagonal, appctx->num_local_dofs, appctx->global_indexing, t8dg_dof_get_double_pointer_to_array (identity_vector),
                  INSERT_VALUES);
  CHKERRQ (ierr);
  ierr = VecAssemblyBegin (mat_diagonal);
  CHKERRQ (ierr);
  ierr = VecAssemblyEnd (mat_diagonal);
  CHKERRQ (ierr);
  ierr = VecReciprocal (mat_diagonal);
  CHKERRQ (ierr);
#endif

  t8dg_debugf ("\nSetting up preconditioner finished\n");
  t8dg_debugf ("Preconditoner looks like\n*************\n");
  VecView (mat_diagonal, PETSC_VIEWER_STDOUT_WORLD);
  t8dg_debugf ("**************\n\n\n");

  //VecSet(mat_diagonal, 1.0);
  pcshell->matrix_diagonal = mat_diagonal;

  t8dg_dof_values_destroy (&identity_vector);
#if 0
#if 0
  //t8dg_timestepping_dirk_ctx_t *appctx;
  t8dg_timestepping_impl_euler_ctx_t *appctx;
  t8dg_dof_values_t  *dof_values_current_copy;
  t8dg_dof_values_t  *identity_vector;
  double             *dof_pointer;
  double             *dof_pointer_identity_vec;
#endif
  Vec                 mat_diagonal;
#if 0
  /* Get the matrix-free application context */
  ierr = MatShellGetContext (pmat, &appctx);
  CHKERRQ (ierr);
#endif
  /* Get the preconditoning context */
  ierr = PCShellGetContext (pc, (void **) &pcshell);
  CHKERRQ (ierr);

  /* Preallocate the mat_diagonal vector and an identity vector */
  ierr = VecDuplicate (u, &mat_diagonal);
  CHKERRQ (ierr);
  ierr = VecSet (mat_diagonal, 1.0);
  CHKERRQ (ierr);
#if 1
  dof_values_current_copy = t8dg_dof_values_clone (current_dofs);
  identity_vector = t8dg_dof_values_duplicate (current_dofs);
  dof_pointer = t8dg_dof_get_double_pointer_to_array (current_dofs);
  for (int i = 0; i < num_total_elements; ++i) {
    dof_pointer[i] = 0.0;
  }
  /* Calculate the diagonal of the application matrix */
  /******************* !!!!!! nur im seriellen Fall möglich !!!!!!!!! *********/
  dof_pointer_identity_vec = t8dg_dof_get_double_pointer_to_array (identity_vector);
  for (int i = 0; i < appctx->num_local_dofs; ++i) {
    dof_pointer[i] = 1.0;
    (appctx->time_derivative_func) (current_dofs, identity_vector, appctx->next_time_point, appctx->user_data);
    t8dg_debugf ("\n*****\nAbleitung wert von identity_vec auswertung an stelle %d ist %f\n*****\n", i, dof_pointer_identity_vec[i]);
    ierr =
      VecSetValue (mat_diagonal, i, 1.0 - ((pcshell->timestep * pcshell->current_b_coeff) * (dof_pointer_identity_vec[i])), INSERT_VALUES);
    CHKERRQ (ierr);
    ierr = VecAssemblyBegin (mat_diagonal);
    CHKERRQ (ierr);
    ierr = VecAssemblyEnd (mat_diagonal);
    CHKERRQ (ierr);
    dof_pointer[i] = 0.0;
  }
  t8dg_dof_values_swap (&current_dofs, &dof_values_current_copy);
  t8dg_dof_values_destroy (&current_dofs);
  t8dg_dof_values_destroy (&identity_vector);
  /* Calculate the reciprocal of the diagonal */
  ierr = VecReciprocal (mat_diagonal);
  CHKERRQ (ierr);
  //VecSet(mat_diagonal, 1.0);
  /* Store the preconditoner */
#endif
  pcshell->matrix_diagonal = mat_diagonal;
#if 0
  t8dg_debugf ("\nSetting up preconditioner finished\n");
  t8dg_debugf ("Preconditoner looks like\n*************\n");
  VecView (mat_diagonal, PETSC_VIEWER_STDOUT_WORLD);
  t8dg_debugf ("**************\n\n\n");
#endif
#endif
  return 0;
}

/* Applies the Jacobi preconditioner to a vector */
PetscErrorCode
JacobiShellPCApply (PC pc, Vec in, Vec out)
{
  PetscErrorCode      ierr;
  t8dg_timestepping_precon_jacobi_ctx_t *pcshell;

  /* Get the preconditoning context */
  ierr = PCShellGetContext (pc, (void **) &pcshell);
  CHKERRQ (ierr);

  /* Pointwise multiplication of the main diagonal to the 'out' vector, equals preconditioning with the diagonal matrix extracted from the system matrix */
  ierr = VecPointwiseMult (out, in, pcshell->matrix_diagonal);
  CHKERRQ (ierr);

  return 0;
}

/* Frees the resources used by the preconditioning methods */
PetscErrorCode
JacobiShellPCDestroy (PC pc)
{
  PetscErrorCode      ierr;
  t8dg_timestepping_precon_jacobi_ctx_t *pcshell;

  /* Get the preconditoning context */
  ierr = PCShellGetContext (pc, (void **) &pcshell);
  CHKERRQ (ierr);

  /* Destroy the vector which stores the main diagonal entries */
  //ierr = VecDestroy(pcshell->matrix_diagonal); CHKERRQ(ierr);

  /* Free the prconditioning context */
  ierr = PetscFree (pcshell);
  CHKERRQ (ierr);

  return 0;
}

#endif

/*********************************************/
