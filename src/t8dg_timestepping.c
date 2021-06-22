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
struct t8dg_timestepping_data
{
  int                 time_order;/**< time order of the Runge-kutta timestepping*/
  double              delta_t;  /**< time step */
  double              t;        /**< current time */
  double              T;        /**< end time */
  double              cfl;      /**< cfl number*/
  int                 step_number;
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
t8dg_timestepping_data_new_cfl (int time_order, double start_time, double end_time, double cfl)
{
  T8DG_ASSERT (time_order > 0);
  t8dg_timestepping_data_t *time_data = T8DG_ALLOC (t8dg_timestepping_data_t, 1);
  time_data->time_order = time_order;
  time_data->t = start_time;
  time_data->T = end_time;
  time_data->cfl = cfl;
  time_data->step_number = 0;
  time_data->delta_t = -1;
  return time_data;
}

t8dg_timestepping_data_t *
t8dg_timestepping_data_new_constant_timestep (int time_order, double start_time, double end_time, double delta_t)
{
  T8DG_ASSERT (time_order > 0);
  t8dg_timestepping_data_t *time_data = T8DG_ALLOC (t8dg_timestepping_data_t, 1);
  time_data->time_order = time_order;
  time_data->t = start_time;
  time_data->T = end_time;
  time_data->cfl = 0;
  time_data->step_number = 0;
  time_data->delta_t = delta_t;
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
  t8dg_dof_values_t  *dof_problem_new;
  PetscInt            its;
  PetscInt           *vec_global_index;
  t8_gloidx_t         global_offset_to_first_local_elem;

  //t8dg_debugf("\n\nInitial dofs passed to implicit Euler:\n\n");
  //t8dg_dof_values_debug_print(*pdof_array);

  /* Store the current timestep */
  appctx.timestep = t8dg_timestepping_data_get_time_step (time_data);

  /* Store the next point in time which is going to be calculated */
  appctx.next_time_point = t8dg_timestepping_data_get_current_time (time_data) + (double) appctx.timestep;

  /* Store the time_derivative_function and the problem related data as well as the degrees of freedom */
  appctx.time_derivative_func = time_derivative;
  appctx.time_data = time_data;
  appctx.user_data = user_data;

  /* Store the degrees of freedom in the matrix-free application context (the current (original) degrees of freedom will be directly overwritten) */
  appctx.future_local_dofs = *pdof_array;
  /* Is going to store the derivation of the degrees of freedom -> the evaluation of the time_derivative function */
  appctx.future_local_dofs_derivative = t8dg_dof_values_duplicate (appctx.future_local_dofs);

  /* the degrees of freedom (of the next time step) which will be filled with the solution */
  dof_problem_new = t8dg_dof_values_duplicate (appctx.future_local_dofs);

  /* Advance a time step */
  t8dg_timestepping_data_set_current_time (time_data, appctx.next_time_point);

  t8dg_debugf ("System wird aufgebaut und geloest\n");

        /********** Setting up and Solving the LS ****************/

  /* Calculate the process-local number of degrees of freedom */
  appctx.num_local_dofs =
    (size_t) (t8dg_dof_get_num_local_elements (appctx.future_local_dofs) * t8dg_dof_get_max_num_element_dof (appctx.future_local_dofs));

  /* Allocate space for an array which holds the global indices of the degrees of freedom */
  ierr = PetscMalloc1 (appctx.num_local_dofs, &vec_global_index);

  /* compute the offset of the first local element concerning the forest and the amout of degrees of freedom per element */
  global_offset_to_first_local_elem = t8_forest_get_first_local_element_id (t8dg_dof_values_get_forest (appctx.future_local_dofs));
  if (global_offset_to_first_local_elem != 0) {
    global_offset_to_first_local_elem = global_offset_to_first_local_elem * t8dg_dof_get_max_num_element_dof (appctx.future_local_dofs);
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

  //t8dg_debugf("\n Right Hand Side f:\n");
  //ierr = VecView(f, PETSC_VIEWER_STDOUT_WORLD);

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
  ierr = KSPSetType (ksp, KSPGMRES);
  CHKERRQ (ierr);
  /* Select a Preconditioner - without Preconditioning */
  ierr = PCSetType (pc, PCNONE);
  CHKERRQ (ierr);
  ierr = KSPSetTolerances (ksp, 1.0e-9, 1.0e-9, PETSC_DEFAULT, PETSC_DEFAULT);
  CHKERRQ (ierr);
  /* Set initial guess to non-zero (to the last solution) otherwise it will filled with zeros */
  ierr = KSPSetInitialGuessNonzero (ksp, PETSC_TRUE);
  CHKERRQ (ierr);
  /* Use eventually further specified properties via the command line, otherwise it is just going to be set up */
  ierr = KSPSetFromOptions (ksp);
  CHKERRQ (ierr);

  t8dg_debugf ("\nGMRES got called\n");

  /* Solve the Linear System */
  ierr = KSPSolve (ksp, f, u);
  CHKERRQ (ierr);

  t8dg_debugf ("\nGMRES solve completed\n");

  /* KSPGetSolution will most likely not be required in this case */
  //ierr = KSPGetSolution(ksp, &u); CHKERRQ(ierr);

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

  t8dg_debugf ("\ndCalculated solution (dofs) of the next timestep:\n");
  t8dg_dof_values_debug_print (dof_problem_new);

  /* Swap the degrees of freedom associated with the problem with the newly calculated degrees of freedom */
  t8dg_dof_values_swap (pdof_array, &dof_problem_new);

  //t8dg_debugf("Dfference between A*u_{k+1} und initial dofs/rhs (u_{k})\n");
  //ierr = VecAXPY(u, -1.0, f); CHKERRQ(ierr);
  //VecView(u, PETSC_VIEWER_STDOUT_WORLD);

  /* Free used resources */
  t8dg_dof_values_destroy (&(appctx.future_local_dofs_derivative));
  t8dg_dof_values_destroy (&dof_problem_new);

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
  t8dg_timestepping_impl_euler_ctx_t *appctx;
  PetscErrorCode      ierr;
  PetscScalar        *current_local_approx;
  int                 dof_iter;
  double             *dof_values_ptr;

  /* Get the matrix-free application context */
  ierr = MatShellGetContext (A, &appctx);
  CHKERRQ (ierr);

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

  /* Calculate the time derivative (concerning the next time step) of the degrees of freedom passed to the MatMult (= 'in' Vector) and store their derivation inside future_local_dofs_derivative */
  (appctx->time_derivative_func) (appctx->future_local_dofs, appctx->future_local_dofs_derivative, appctx->next_time_point,
                                  appctx->user_data);

  /* calculate the sum of the initially passed degrees of freedom ('in' Vector) and their (scaled) derivation */
  t8dg_dof_values_axpy (-(appctx->timestep), appctx->future_local_dofs_derivative, appctx->future_local_dofs);

  /* Write the result of the matrix-application to the 'in' Vector in the 'out' Vector */
  ierr =
    VecSetValues (out, appctx->num_local_dofs, appctx->global_indexing, t8dg_dof_get_double_pointer_to_array (appctx->future_local_dofs),
                  INSERT_VALUES);
  CHKERRQ (ierr);
  ierr = VecAssemblyBegin (out);
  CHKERRQ (ierr);
  ierr = VecAssemblyEnd (out);
  CHKERRQ (ierr);

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

  /* Clone the current degrees of freedom */
  appctx.local_derivation_degrees = *pdof_array;
  appctx.current_local_dofs = t8dg_dof_values_clone (appctx.local_derivation_degrees);

  /* Duplicate the other t8dg_dof_values_t */
  appctx.future_local_dofs_derivative = t8dg_dof_values_duplicate (appctx.local_derivation_degrees);
  appctx.future_local_dofs_step = t8dg_dof_values_duplicate (appctx.local_derivation_degrees);
  appctx.future_local_dofs = t8dg_dof_values_duplicate (appctx.local_derivation_degrees);
  t8dg_dof_values_set_zero (appctx.future_local_dofs);

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
    (size_t) (t8dg_dof_get_num_local_elements (appctx.future_local_dofs) * t8dg_dof_get_max_num_element_dof (appctx.future_local_dofs));

  /* Allocate space for an array which holds the global indices of the degrees of freedom */
  ierr = PetscMalloc1 (appctx.num_local_dofs, &vec_global_index);

  /* compute the offset of the first local element concerning the forest and the amout of degrees of freedom per element */
  global_offset_to_first_local_elem = t8_forest_get_first_local_element_id (t8dg_dof_values_get_forest (appctx.future_local_dofs));
  if (global_offset_to_first_local_elem != 0) {
    global_offset_to_first_local_elem = global_offset_to_first_local_elem * t8dg_dof_get_max_num_element_dof (appctx.future_local_dofs);
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
  //ierr = VecSetValues(f, appctx.num_local_dofs, vec_global_index, t8dg_dof_get_double_pointer_to_array(*pdof_array), INSERT_VALUES); CHKERRQ(ierr);
  ierr = VecSet (f, 0.0);
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
  ierr = VecSet (u, 1.0);
  CHKERRQ (ierr);
  ierr = VecAssemblyBegin (u);
  CHKERRQ (ierr);
  ierr = VecAssemblyEnd (u);
  CHKERRQ (ierr);
  /* Copy the entries of f to u - u is used as an initial guess; otherwise u will be overwitten with zeros during KSPSolve */
  //ierr = VecCopy (f, u);
  //CHKERRQ (ierr);
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
  ierr = PCSetType (pc, PCNONE);
  CHKERRQ (ierr);
  ierr = KSPSetTolerances (ksp, 1.0e-12, 1.0e-12, PETSC_DEFAULT, PETSC_DEFAULT);
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
    appctx.next_time_point = current_time + ((*dirk_c_coeff)[step_iter]) * appctx.timestep;
    t8dg_timestepping_data_set_current_time (time_data, appctx.next_time_point);

    /* Set the current coefficients needed in the matrix-free application for the next step */
    appctx.dirk_current_a_coeff = (*dirk_a_coeff)[step_iter][step_iter];
    //appctx.dirk_current_b_coeff = (*dirk_b_coeff)[step_iter];

    /* Assemble the (pre) degrees of freedom which will be derivated in the MatVec product */
    //t8dg_dof_values_set_zero (appctx.future_local_dofs_step);
    t8dg_dof_values_copy (appctx.current_local_dofs, appctx.future_local_dofs_step);
    for (iter = 0; iter < step_iter; ++iter) {
      t8dg_dof_values_axpy (appctx.timestep * ((*dirk_a_coeff)[step_iter][iter]), local_dofs_step[iter], appctx.future_local_dofs_step);
    }

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
    dof_values_ptr = t8dg_dof_get_double_pointer_to_array (appctx.future_local_dofs_step);

    /* Overwrite the dof_problem_new with the newly calculated degrees of freedom */
    for (iter = 0; iter < appctx.num_local_dofs; ++iter) {
      dof_values_ptr[iter] = (double) new_local_approx[iter];
    }

    /* Restore the array entries */
    ierr = VecRestoreArrayRead (u, &new_local_approx);
    CHKERRQ (ierr);

    /* Create a new empty t8dg_dof_values_t and swap with the current solution */
    local_dofs_step[step_iter] = t8dg_dof_values_duplicate (appctx.future_local_dofs_step);

    /* Store the approximation of the degrees of freedom of the current step; they are the k_{j} of the different steps */
    t8dg_dof_values_swap (&(appctx.future_local_dofs_step), &local_dofs_step[step_iter]);

  }

  /* Add the values of the stpes to the former approximation to obtain the degrees of freedom of the next step */
  for (iter = 0; iter < num_order_stages; ++iter) {
    t8dg_dof_values_axpy ((appctx.timestep * (*dirk_b_coeff)[iter]), local_dofs_step[iter], appctx.current_local_dofs);
  }

  /* Swap the former degrees of freedom with the newly calculated degrees of freedom */
  t8dg_dof_values_swap (pdof_array, &(appctx.current_local_dofs));

  t8dg_debugf ("\nDIRK Method completed\n");

  /* Free used resources */
  t8dg_dof_values_destroy (&(appctx.future_local_dofs));
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
  dof_values_ptr = t8dg_dof_get_double_pointer_to_array (appctx->future_local_dofs);

  /* Overwrite them with the 'in' Vector of the Matrix Multiplication */
  for (dof_iter = 0; dof_iter < appctx->num_local_dofs; ++dof_iter) {
    dof_values_ptr[dof_iter] = current_local_approx[dof_iter];
  }

  /* Restore the array entries of the 'in' Vector */
  ierr = VecRestoreArrayRead (in, &current_local_approx);
  CHKERRQ (ierr);

  /* Overwrite the degrees of freedom in the problem to obtain the derivation of the degrees of freedom of the 'in' Vector */
  t8dg_dof_values_axpyz ((appctx->timestep * appctx->dirk_current_a_coeff), appctx->future_local_dofs, appctx->future_local_dofs_step,
                         appctx->local_derivation_degrees);

  /* Get the derivation of the degrees of freedom of the current RKV step and store them in appctx->future_local_dofs_derivative */
  (appctx->time_derivative_func) (appctx->local_derivation_degrees, appctx->future_local_dofs_derivative, appctx->next_time_point,
                                  appctx->user_data);

  /* Subtract the derivation of the degrees from the dofs of the 'in' Vector */
  t8dg_dof_values_axpy (-1.0, appctx->future_local_dofs_derivative, appctx->future_local_dofs);

  /* Write the result of the matrix-application to the 'in' Vector in the 'out' Vector */
  ierr =
    VecSetValues (out, appctx->num_local_dofs, appctx->global_indexing, t8dg_dof_get_double_pointer_to_array (appctx->future_local_dofs),
                  INSERT_VALUES);
  CHKERRQ (ierr);
  ierr = VecAssemblyBegin (out);
  CHKERRQ (ierr);
  ierr = VecAssemblyEnd (out);
  CHKERRQ (ierr);

  return 0;
}
#endif
