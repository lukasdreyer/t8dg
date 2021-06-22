#include "t8dg.h"
#include "t8dg_timestepping.h"
#include "t8dg_dof.h"

#if 0
#include <petscksp.h>
#include <petscdmda.h>
#include <petscdm.h>
#endif

/* Help Message */
static char         help[] = "First matrix-free PETSc Example";

#if T8_WITH_PETSC
typedef struct
{
  PetscScalar         example_coeff;
  PetscInt            problem_size;
  Vec                 example_vec;
  DM                  element_data;
} ExampleOperatorCtx;
#endif

/* Struct which keeps information regarding the matrix-free application of system matrix resulting from the implicit Euler-Method */
#if 0
typedef struct
{
  PetscScalar         timestep;
  double next_time_point;
  size_t num_local_dofs;
  t8dg_time_matrix_application time_derivative_func;
  t8dg_dof_values_t  *future_local_dofs;
  t8dg_dof_values_t  *future_local_dofs_derivative;
  t8dg_timestepping_data_t *time_data;
  void *user_data;
} t8dg_timestepping_impl_euler_ctx_t;
#endif

#if T8_WITH_PETSC
/* Matrix free Example operator: applies a diagonal matrix */
extern PetscErrorCode MatMult_MF_Operator (Mat, Vec, Vec);

#if 1
/* Matrix-free routine that mimics the application of the matrix resulting from the Implicit/backwards Euler_Method */
extern PetscErrorCode MatMult_MF_Impl_Euler (Mat, Vec, Vec);
#endif

#endif

#if 1
void
t8dg_timestepping_implicit_euler (t8dg_time_matrix_application time_derivative,
                                  t8dg_timestepping_data_t * time_data, t8dg_dof_values_t ** pdof_array, void *user_data)
{
#if T8_WITH_PETSC
  size_t iter;
  double             *raw_local_dofs;
  Mat                 A;
  Vec                 u, f;
  KSP ksp;
  PC pc;
  PetscScalar        *local_dofs;
  PetscErrorCode ierr;
  t8dg_timestepping_impl_euler_ctx_t appctx;
  
  /* Store the current timestep */
  appctx.timestep = t8dg_timestepping_data_get_time_step (time_data);

  /* Store the next point in time which is going to be calculated */
  appctx.next_time_point = t8dg_timestepping_data_get_current_time (time_data) + (double) appctx.timestep;

  /* Store the time_derivative_function and the problem related data as well as the degrees of freedom */
  appctx.time_derivative_func = time_derivative;
  appctx.time_data = time_data;
  appctx.user_data = user_data;

  /* Store the degrees of freedom in the matrix-free application context */
  appctx.future_local_dofs = t8dg_dof_values_clone (*pdof_array);
  appctx.future_local_dofs_derivative = t8dg_dof_values_duplicate (appctx.future_local_dofs);


        /********** Setting up and Solving the LS ****************/

  /* Setting up Vectors */

  /* Get the whole number of the process-local degrees of freedom */
  /* in all axpy-routines etc in t8dg, num_total_elem * max_num_dofs is used as the size of dof, but what is num_local_elements then? num_total_elem is eventually with ghost cells */
  /* i dont know if the programs runs correctly if the ghost cells are included */
  appctx.num_local_dofs = (size_t) (t8dg_dof_get_num_total_elements(appctx.future_local_dofs) * t8dg_dof_get_max_num_element_dof(appctx.future_local_dofs));
  /* Allocate the space for all process-local degrees of freedom */
  ierr = PetscMalloc1(appctx.num_local_dofs, &local_dofs); CHKERRQ(ierr);

  /* Retrieve a double pointer to the process-local degrees of freedom stored in a sc_array_t */
  raw_local_dofs = t8dg_dof_get_double_pointer_to_array(appctx.future_local_dofs);

  /* Store the local degrees of freedom in the 'local_dofs' array which will be owned by the global Vec */
  /* in all axpy-routines etc in t8dg, num_total_elem * max_num_dofs is used as the size of dof, but what is num_local_elements then? */
  for (iter = 0; iter < appctx.num_local_dofs; ++iter)
  {
    local_dofs[iter] = (PetscScalar) raw_local_dofs[iter];
  }

  /* Create the global Vec */
  ierr = VecCreateMPIWithArray(PETSC_COMM_WORLD, 1, (PetscInt) appctx.num_local_dofs, PETSC_DETERMINE, local_dofs, &f); CHKERRQ(ierr);
  //ierr = VecSetType (f, VECMPI);
  //CHKERRQ (ierr);
  //ierr = VecCreate (PETSC_COMM_WORLD, &f);
  //CHKERRQ (ierr);
  /* is the size correct or does num_local_elements have to be multiplied by the number of degrees of freedom per element */
  //ierr = VecSetSizes (f, (PetscInt) (dof_beginning->num_local_elements * dof_beginning->max_num_element_dof), (PetscInt) (dof_beginning->num_total_elements * dof_beginning->max_num_element_dof));
  //CHKERRQ (ierr);
  /* Only local values are written into the Vec, there shouldn't be any cross-process value exchanges during the assembly */
  /* ierr = VecSetOption(f, VEC_IGNORE_OFF_PROC_ENTRIES, PETSC_TRUE); CHKERRQ(ierr); */
  //ierr = VecSetType (f, VECMPI);
  //CHKERRQ (ierr);
  ierr = PetscObjectSetName ((PetscObject) f, "Right Hand Side");
  CHKERRQ (ierr);
  /* Duplicate the setup-options for the approximation vector but place_array have to be called to change the array dependency? There shouldnt be any dependencies, vec_duplicate ia rather fpr allocation purposees, vec entries will not be copied */
  ierr = VecDuplicate (f, &u);
  CHKERRQ (ierr);
  ierr = PetscObjectSetName ((PetscObject) u, "Approximation");
  CHKERRQ (ierr);



  /* Setting up matrix-free Matrix */
  //ierr = MatCreate (PETSC_COMM_WORLD, &A);
  //CHKERRQ (ierr);
  //ierr = MatSetSizes (A, PETSC_DECIDE, PETSC_DECIDE, PETSC_DETERMINE, PETSC_DETERMINE);
  //CHKERRQ (ierr);               //Sizes of process local or global components has to be determined/calculated
  //ierr = MatSetType (A, MATSHELL);
  //CHKERRQ (ierr);
  //ierr = MatSetUp (A);
  //CHKERRQ (ierr);
  /* Create a matrix shell with local dimensions equal to the dimension of the Vec containing the process-local degrees of freedom and add an application context needed by the matrix-free MatVec multiplication */
  ierr = MatCreateShell(PETSC_COMM_WORLD, (PetscInt) appctx.num_local_dofs, (PetscInt) appctx.num_local_dofs, PETSC_DETERMINE, PETSC_DETERMINE, (void *) &appctx, &A); CHKERRQ(ierr);
  ierr = PetscObjectSetName ((PetscObject) A, "Matrix-Free Application-Matrix of Imp Euler");
  CHKERRQ (ierr);
  /* Set a Context with furhter information/data */
  //ierr = MatShellSetContext (A, (void *) &appctx);
  //CHKERRQ (ierr);
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
  ierr = KSPSetTolerances (ksp, 1.e-7, PETSC_DEFAULT, PETSC_DEFAULT, PETSC_DEFAULT);
  CHKERRQ (ierr);
  ierr = KSPSetFromOptions (ksp);
  CHKERRQ (ierr);

  /* Solve the Linear System */
  ierr = KSPSolve (ksp, f, u);
  CHKERRQ (ierr);
        /*********************************************************/

  /* Advance a time step */
  t8dg_timestepping_data_set_current_time (time_data, appctx.next_time_point);

  /* Update the new degrees of freedom */
  t8dg_dof_values_swap (pdof_array, &(appctx.future_local_dofs_derivative));

  /* Free used resources */
  t8dg_dof_values_destroy (&(appctx.future_local_dofs_derivative));
  t8dg_dof_values_destroy (&(appctx.future_local_dofs));

  /* Free PETSc resources */
  ierr = MatDestroy (&A);
  CHKERRQ (ierr);
  ierr = VecDestroy (&f);
  CHKERRQ (ierr);
  ierr = VecDestroy (&u);
  CHKERRQ (ierr);
  ierr = KSPDestroy (&ksp);
  CHKERRQ (ierr);
  ierr = PetscFree(local_dofs); CHKERRQ(ierr);

#else
  t8dg_global_productionf
    ("t8code/t8dg is currently not configured with PETSc.\n\nIn order to use this function t8dg needs to be configured with PETSc.\n");
#endif
}
#endif

int
main (int argc, char **args)
{
#if T8_WITH_PETSC
  ExampleOperatorCtx  appctx;
  DM                  grid_data;
  Vec                 x, b, u;
  Mat                 A;
  KSP                 ksp;
  PC                  pc;
  PetscErrorCode      ierr;
  PetscInt            num_iterations;
  PetscReal           norm;

  PetscScalar *ex_array;
  Vec ex_vec;
  /* Initialize PETSc */
  ierr = PetscInitialize (&argc, &args, (char *) 0, help);
  if (ierr)
    return ierr;

  ierr = PetscMalloc1(6, &ex_array); CHKERRQ(ierr);
  for (int i=0; i< 6; ++i)
  {
    ex_array[i] = i;
  }
  ierr = VecCreateMPIWithArray(PETSC_COMM_WORLD, 1, 6, 12, ex_array, &ex_vec); CHKERRQ(ierr);
  ierr = VecView(ex_vec, PETSC_VIEWER_STDOUT_WORLD);
  /* Fill the application context with example values */
  appctx.example_coeff = 2.0;
  appctx.problem_size = 20;

  /* Create distributed memory array */
  DMDACreate1d (PETSC_COMM_WORLD, DM_BOUNDARY_NONE, appctx.problem_size, 1, 1, NULL, &appctx.element_data);
  DMSetFromOptions (appctx.element_data);
  DMSetUp (appctx.element_data);

  /* Sizes of matrix-free Matrix */
  PetscInt            m, n, M = appctx.problem_size, N = appctx.problem_size;

  /* Create Vectors */
  ierr = VecCreate (PETSC_COMM_WORLD, &x);
  CHKERRQ (ierr);
  ierr = PetscObjectSetName ((PetscObject) x, "Approximation");
  CHKERRQ (ierr);
  ierr = VecSetSizes (x, PETSC_DECIDE, appctx.problem_size);
  CHKERRQ (ierr);
  ierr = VecSetFromOptions (x);
  CHKERRQ (ierr);

  ierr = VecDuplicate (x, &b);
  CHKERRQ (ierr);
  ierr = PetscObjectSetName ((PetscObject) b, "RightHandSide");
  CHKERRQ (ierr);

  ierr = VecDuplicate (x, &u);
  CHKERRQ (ierr);
  ierr = PetscObjectSetName ((PetscObject) u, "ExactSolution");
  CHKERRQ (ierr);

  /* Create a matrix-free Matrix */
  ierr = MatCreate (PETSC_COMM_WORLD, &A);
  CHKERRQ (ierr);
  ierr = MatSetSizes (A, PETSC_DECIDE, PETSC_DECIDE, M, N);
  CHKERRQ (ierr);
  ierr = MatSetType (A, MATSHELL);
  CHKERRQ (ierr);
  ierr = MatSetUp (A);
  CHKERRQ (ierr);

  /* Set a Context with furhter information/data */
  ierr = MatShellSetContext (A, (void *) &appctx);
  CHKERRQ (ierr);

  /* Define the (multiplicative) MatVec-Operation */
  ierr = MatShellSetOperation (A, MATOP_MULT, (void (*)(void)) MatMult_MF_Operator);
  CHKERRQ (ierr);

  /* Set exact solution and right hand side; solution is set so 2.5; cooresponding to this solution the rhs is calculated */
  ierr = VecSet (u, 2.5);
  CHKERRQ (ierr);
  ierr = MatMult (A, u, b);
  CHKERRQ (ierr);

  ierr = VecView (b, PETSC_VIEWER_STDOUT_WORLD);

  /* Create KSP */
  KSPCreate (PETSC_COMM_WORLD, &ksp);

  /* Set KSP: Set Operators; here: matrix that defines LS also serve as the matrix that definies the preconditioner */
  ierr = KSPSetOperators (ksp, A, A);
  CHKERRQ (ierr);

  /* Extract KSP and PC Context from  KSP Context */
  ierr = KSPGetPC (ksp, &pc);
  CHKERRQ (ierr);

  /* Select GMRES as Solver; default: Restarted GMRES(30) */
  KSPSetType (ksp, KSPGMRES);

  /* Select a Preconditioner - without Preconditioning */
  ierr = PCSetType (pc, PCNONE);
  CHKERRQ (ierr);

  ierr = KSPSetTolerances (ksp, 1.e-7, PETSC_DEFAULT, PETSC_DEFAULT, PETSC_DEFAULT);
  CHKERRQ (ierr);

  ierr = KSPSetFromOptions (ksp);
  CHKERRQ (ierr);

  /* Solve the Linear System */
  ierr = KSPSolve (ksp, b, x);
  CHKERRQ (ierr);

  /* Extract Number of Iterations needed */
  KSPGetIterationNumber (ksp, &num_iterations);

  ierr = (double) VecAXPY (x, -1.0, u);
  CHKERRQ (ierr);
  ierr = VecNorm (x, NORM_2, &norm);
  CHKERRQ (ierr);

  /* Print out norm and iteration count */
  ierr = PetscPrintf (PETSC_COMM_WORLD, "GMRES: Error Norm: %g, Iterations: %d\n\n", (double) norm, num_iterations);
  CHKERRQ (ierr);

  /* View Solver info */
  ierr = KSPView (ksp, PETSC_VIEWER_STDOUT_WORLD);
  CHKERRQ (ierr);

  /* View Matrix info */
  ierr = MatView (A, PETSC_VIEWER_STDOUT_WORLD);
  CHKERRQ (ierr);

  /* Free resources */
  ierr = VecDestroy (&x);
  CHKERRQ (ierr);
  ierr = VecDestroy (&u);
  CHKERRQ (ierr);
  ierr = VecDestroy (&b);
  CHKERRQ (ierr);
  ierr = MatDestroy (&A);
  CHKERRQ (ierr);
  ierr = KSPDestroy (&ksp);
  CHKERRQ (ierr);

  /* Finalize PETSc */
  ierr = PetscFinalize ();
  return ierr;
#else
  t8dg_global_productionf
    ("t8code/t8dg is currently not configured with PETSc.\n\nIn order to use this function t8dg needs to be configured with PETSc.\n");
  return 0;
#endif
}

/* Function that mimics the explicit multiplication of A*in (of the matrix-free matrix) */
#if T8_WITH_PETSC
PetscErrorCode
MatMult_MF_Operator (Mat A, Vec in, Vec out)
{
  ExampleOperatorCtx *appctx;
  Vec                 local_data;
  Vec                 local_values;
  PetscErrorCode      ierr;
  Vec                 local_indices;

  /* Get the matrix-free application context and the local element data */
  ierr = MatShellGetContext (A, &appctx);
  CHKERRQ (ierr);
  ierr = DMGetLocalVector (appctx->element_data, &local_data);
  CHKERRQ (ierr);
  ierr = DMGlobalToLocalBegin (appctx->element_data, in, INSERT_VALUES, local_data);
  CHKERRQ (ierr);
  ierr = DMGlobalToLocalEnd (appctx->element_data, in, INSERT_VALUES, local_data);
  CHKERRQ (ierr);

  ierr = DMGetLocalVector (appctx->element_data, &local_values);
  CHKERRQ (ierr);
  ierr = VecSet (local_values, 0.0);
  CHKERRQ (ierr);
  /* Applying diagonal matrix with a_(i,i) = example_coeff */
  ierr = VecAXPY (local_values, appctx->example_coeff, local_data);
  CHKERRQ (ierr);

  /* Construct global output Vector */
  VecSet (out, 0.0);
  DMLocalToGlobalBegin (appctx->element_data, local_values, INSERT_VALUES, out);
  DMLocalToGlobalEnd (appctx->element_data, local_values, INSERT_VALUES, out);

  /* Free resources */
  DMRestoreLocalVector (appctx->element_data, &local_data);
  DMRestoreLocalVector (appctx->element_data, &local_values);
  return 0;

}
#endif

#if 1
/* Function that mimics the multiplication of the system matrix resulting from the implicit Euler_Method */
#if T8_WITH_PETSC
PetscErrorCode
MatMult_MF_Impl_Euler (Mat A, Vec in, Vec out)
{
  t8dg_timestepping_impl_euler_ctx_t *appctx;
  PetscErrorCode      ierr;
  PetscScalar         impl_euler_coeff;
  PetscScalar *current_local_approx;
  PetscInt local_range_start, local_range_end, iter;
  int dof_iter;
  double *dof_values_ptr;

  /* Get the matrix-free application context */
  ierr = MatShellGetContext (A, &appctx);
  CHKERRQ (ierr);

  /* Get the size of the timestep */
  impl_euler_coeff = -1.0 * appctx->timestep;

  /* Use the 't8dg_dof_values_t' as a shell for new calculated Values to switch between (Petsc) Vec and t8dg_dof_values_t */

  /* u_{k+1} - (delta_t} * R(u_{k+1}) = u_{k} is the system to be solved, therefore A*u_{k+1} equals the application u_{k+1} - delta_{t}*R(u_{k+1}) */
  /* therefore, it's just the execution of VecWAXPY(w,a,x,y) => w = a*x + y */
  /* the evaluation of time_derivative(...) needs the dofs_values, but this application needs an Vec array */

  /* Retrieve the local part of the 'in' Vector for reading purposes */
  ierr = VecGetArrayRead(in, &current_local_approx); CHKERRQ(ierr);

  /* Get a double pointer to the data inside the sc_array_t of t8dg_dof_values_t */
  dof_values_ptr = t8dg_dof_get_double_pointer_to_array(appctx->future_local_dofs);

  /* Overwrite the t8dg_dof_values_t, so that the current approximation can be evaluted by 'time_derivative' */
  for (dof_iter = 0; dof_iter < appctx->num_local_dofs; ++dof_iter)
  {
    dof_values_ptr[dof_iter] = (double) current_local_approx[dof_iter];
  }
  /* Restore the array entries */
  ierr = VecRestoreArrayRead(in, &current_local_approx); CHKERRQ(ierr);

  /* call time_derivative function */
  (appctx->time_derivative_func) (appctx->future_local_dofs, appctx->future_local_dofs_derivative, appctx->next_time_point, appctx->user_data);

  /* Convert the time_derivative evaluation back to a PETSc Vec */
  ierr = VecGetOwnershipRange(out, &local_range_start, &local_range_end); CHKERRQ(ierr);

  /* Get a double pointer to the data inside the sc_array_t of t8dg_dof_values_t */
  dof_values_ptr = t8dg_dof_get_double_pointer_to_array(appctx->future_local_dofs_derivative);

  for (iter = local_range_start, dof_iter = 0; iter < local_range_end; ++iter, ++dof_iter)
  {
    ierr = VecSetValue(out, iter, (PetscScalar) dof_values_ptr[dof_iter], INSERT_VALUES); CHKERRQ(ierr);
  }
  ierr = VecAssemblyBegin(out); CHKERRQ(ierr);
  ierr = VecAssemblyEnd(out); CHKERRQ(ierr);

  /* Use VecAYPX to scale and add the rest of the MatVec Application to the 'out' vector */
  ierr = VecAYPX (out, impl_euler_coeff, in);
  CHKERRQ (ierr);

  return 0;
}
#endif

#endif
