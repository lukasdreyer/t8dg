#include "t8dg.h"
#include "t8dg_timestepping.h"

#if T8_WITH_PETSC
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
#if T8_WITH_PETSC
typedef struct
{
  PetscScalar         timestep;
  t8dg_time_matrix_application time_derivative_func;
} t8dg_timestepping_impl_euler_ctx;
#endif

#if T8_WITH_PETSC
/* Matrix free Example operator: applies a diagonal matrix */
extern PetscErrorCode MatMult_MF_Operator (Mat, Vec, Vec);

#if 0
/* Matrix-free routine that mimics the application of the matrix resulting from the Implicit/backwards Euler_Method */
extern PetscErrorCode MatMult_MF_Impl_Euler (Mat, Vec, Vec);
#endif

#endif

#if 0
void
t8dg_timestepping_implicit_euler (t8dg_time_matrix_application time_derivative,
                                  t8dg_timestepping_data_t * time_data, t8dg_dof_values_t ** pdof_array, void *user_data)
{
#if T8_WITH_PETSC
  double              time_beginning;
  t8dg_dof_values_t  *dof_beginning;
  t8dg_dof_values_t  *dof_change;
  t8dg_dof_values_t  *dof_new;
  Mat                 A;
  Vec                 u, f;
  t8dg_timestepping_impl_euler_ctx appctx;

  /* Get the current time */
  time_beginning = t8dg_timestepping_data_get_current_time (time_data);

  /* Store the current timestep */
  appctx.timestep = t8dg_timestepping_data_get_time_step (time_data);

  /* Store the time_derivative_function resultung from the PDE */
  appctx.time_derivative_func = time_derivative
    /* Copy the degrees of freedom */
    dof_beginning = t8dg_dof_values_clone (*pdof_array);
  dof_new = t8dg_dof_values_duplicate (dof_beginning);
  dof_change = t8dg_dof_values_duplicate (dof_beginning);

        /********** Setting up and Solving the LS ****************/

  /* Setting up Vectors */
  /* maybe one can use the sc_array inside the dofs and use veccreatempiwitharray() */
  ierr = VecCreate (PETSC_COMM_WORLD, &f);
  CHKERRQ (ierr);
  /* is the size correct or does num_local_elements have to be multiplied by the number of degrees of freedom per element */
  ierr = VecSetSizes (f, (PetscInt) dof_beginning->num_local_elements, (PetscInt) dof_beginning->num_total_elements);
  CHKERRQ (ierr);
  /* Only local values are written into the Vec, there shouldn't be any cross-process value exchange during the assembly */
  /* ierr = VecSetOption(f, VEC_IGNORE_OFF_PROC_ENTRIES, PETSC_TRUE); CHKERRQ(ierr); */
  ierr = VecSetType (f, VECMPI);
  CHKERRQ (ierr);
  ierr = VecDuplicate (f, &u);
  CHKERRQ (ierr);
  /* Storing the degrees of freedom inside the Vecs */

  /* Setting up matrix-free Matrix */
  ierr = MatCreate (PETSC_COMM_WORLD, &A);
  CHKERRQ (ierr);
  ierr = MatSetSizes (A, PETSC_DECIDE, PETSC_DECIDE, PETSC_DETERMINE, PETSC_DETERMINE);
  CHKERRQ (ierr);               //Sizes of process local or global components has to be determined/calculated
  ierr = MatSetType (A, MATSHELL);
  CHKERRQ (ierr);
  ierr = MatSetUp (A);
  CHKERRQ (ierr);
  /* Set a Context with furhter information/data */
  ierr = MatShellSetContext (A, (void *) &appctx);
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
  ierr = KSPSetTolerances (ksp, 1.e-7, PETSC_DEFAULT, PETSC_DEFAULT, PETSC_DEFAULT);
  CHKERRQ (ierr);
  ierr = KSPSetFromOptions (ksp);
  CHKERRQ (ierr);

  /* Solve the Linear System */
  ierr = KSPSolve (ksp, b, x);
  CHKERRQ (ierr);
        /*********************************************************/

  /* Advance a time step */
  t8dg_timestepping_data_set_current_time (time_data, time_beginning + (double) appctx.timestep);

  /* Update the new degrees of freedom */
  t8dg_dof_values_swap (pdof_array, &dof_new);

  /* Free used resources */
  t8dg_dof_values_destroy (&dof_beginning);
  t8dg_dof_values_destroy (&dof_change);
  t8dg_dof_values_destroy (&dof_new);

  /* Free PETSc resources */
  ierr = MatDestroy (&A);
  CHKERRQ (ierr);
  ierr = VecDestroy (&f);
  CHKERRQ (ierr);
  ierr = VecDestroy (&u);
  CHKERRQ (ierr);
  ierr = KSPDestroy (&ksp);
  CHKERRQ (ierr);

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

  /* Initialize PETSc */
  ierr = PetscInitialize (&argc, &args, (char *) 0, help);
  if (ierr)
    return ierr;

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

#if 0
/* Function that mimics the multiplication of the system matrix resulting from the implicit Euler_Method */
#if T8_WITH_PETSC
PetscErrorCode
MatMult_MF_Impl_Euler (Mat A, Vec in, Vec out)
{
  t8dg_timestepping_impl_euler_ctx *appctx;
  t8dg_time_matrix_application time_derivative;
  PetscErrorCode      ierr;
  PetscScalar         impl_euler_coeff;
  /* maybe should copy the 'dofs' of the problem into an petsc vector, otherwise im not sure yet how to describe the global vector partioned over the processes for GMRES */
  /* this is problem in the domain of the creation of the system, not in the operator execution */

  /* Get the matrix-free application context */
  ierr = MatShellGetContext (A, &appctx);
  CHKERRQ (ierr);

  /* Get the size of the timestep */
  impl_euler_coeff = -1 * appctx->timestep;

  /* Get the current function pointer to the time derivative function of the PDE */
  time_derivative = appctx->time_derivative_func;

  /* u_{k+1} - (delta_t} * R(u_{k+1}) = u_{k} is the system to be solved, therefore A*u_{k+1} equals the application u_{k+1} - delta_{t}*R(u_{k+1}) */
  /* therefore, it's just the execution of VecWAXPY(w,a,x,y) => w = a*x + y */
  /* the evaluation of time_derivative(...) needs the dofs_values, but this application needs an Vec array */
  ierr = VecWAXPY (out, impl_euler_coeff, time_derivative (in), in);
  CHKERRQ (ierr);

  return 0;
}
#endif

#endif
