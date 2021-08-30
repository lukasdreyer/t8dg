/** @file t8dg_preconditioner.h */
/* This header hold the structs and functions needed by the matrix preconditioner routines */

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

/* Typedef which resembles an adaption function which is needed in order to adapt a forest */
typedef int         (*t8dg_corase_lvl_adapt_func_t) (t8_forest_t forest, t8_forest_t forest_from, t8_locidx_t itree, t8_locidx_t ielement,
                                                     t8_eclass_scheme_c * ts, int num_elements, t8_element_t * elements[]);

/* Typedef which resembles a pointer to a function describing the time derivation of the degrees of freedom; same as in t8dg_timestepping.h */
typedef void        (*t8dg_time_derivation_matrix_application_t) (t8dg_dof_values_t * src_dof, t8dg_dof_values_t * dest_dof, const double t,
                                                                  const void *application_data);

/* typedef which resembles a pointer to a function describing the application of the problem in concerns of a choosen block-preconditioner */
typedef void        (*t8dg_block_precon_matrix_application_t) (t8dg_dof_values_t * src_dof, t8dg_dof_values_t * dest_dof, const double t,
                                                               const void *application_data, int selection);

/* Struct needed by the matrix application in order to solve for the inverse application of the block preconditioner */
typedef struct
{
  t8dg_dof_values_t **problem_dofs;
  t8dg_dof_values_t  *problem_dofs_derivation;
  void               *problem;
  t8dg_block_precon_matrix_application_t jacobi_preconditioner_advect_diff_application;
  double              current_implicit_coefficient;
  double              current_point_in_time;
  double              timestep;
  size_t              num_local_dofs;
  PetscInt           *local_indexing;
  int                 selection;
} t8dg_precon_block_matrix_ctx_t;

/* Struct that keeps information needed by the Jacobi Preconditioner */
typedef struct
{
  t8dg_precon_block_matrix_ctx_t jacobi_ctx;
  Mat                 M_jacobi;
  Vec                 u, f;
  KSP                 jacobi_preconditioner;
  PC                  pc_jacobi_preconditioner;
  PetscInt           *global_indexing;

} t8dg_block_preconditioner_ctx_t;

/* Struct that keeps the context needed in order to establish a coarser mesh on which the multigrid preconditioner performs */
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

/* Struct that keeps the context needed by the solver of the coarse level which performs a matrix-free calculation on the restricted system */
typedef struct
{
  t8dg_mg_coarse_lvl_t *coarse_lvl;
  void               *user_data;
  t8dg_dof_values_t  *problem_dofs_derivation;
  t8dg_time_matrix_application time_derivative_func;
  double              current_point_in_time;
  double              timestep;
  double              current_coefficient;
} t8dg_coarse_matrix_ctx_t;

/* Struct the holds information needed by restriction and prolongation routines to transfer between fine and coarse level */
typedef struct
{
  t8dg_linear_advection_diffusion_problem_t *problem;
  t8dg_mg_coarse_lvl_t *coarse_lvl;
  size_t              num_local_dofs_coarse_grid;
  size_t              num_local_dofs_fine_grid;
  PetscInt           *fine_lvl_global_indexing;
} t8dg_mg_interpolating_ctx_t;

/* Struct that fulfills a general purpose and holds Matrices, Contexts, etc. for every possible preconditioner */
typedef struct
{
  /* BEGIN Two-Level Multigrid components */
  KSP                 coarse_solver;
  KSP                 smoother;
  PC                  coarse_pc;
  PC                  smoother_pc;
  Mat                 Restriction, Prolongation;
  Mat                 A_coarse;
  t8dg_coarse_matrix_ctx_t cmat_ctx;
  t8dg_mg_interpolating_ctx_t res_prol_ctx;
  t8dg_mg_coarse_lvl_t coarse_lvl;
  /* END Two-Level Multigrid components */

  /* BEGIN Jacobi-Preconditioner components */
  t8dg_block_preconditioner_ctx_t *jacobi_preconditioner_ctx;
  /* END Jacobi-Preconditioner components */

  double              preconditioner_setup_time;
} t8dg_precon_general_preconditioner_t;

/** This function initializes a selected preconditioner
* \param[in, out] pc A pointer to a PETSc \a pc which will be set up with the selected preconditioner
* \param[in] selector An integer describing which preconditioner ought to be initialized and set up
* \param[in, out] general_precon A pointer to a \a t8dg_precon_general_preconditioner_t context holding further information needed by the matrix application regarding of the choice of the preconditioner
* \param[in] problem A pointer to the advect diff problem
* \param[in] problem_dofs A pointer to a pointer which relates to the degrees of freedom of the initial advect diff problem \a problem
* \param[in] time_derivative A function describing the time derivation of the coefficients (in the degrees of freedom \a problem_dofs); resulting from the discretization of the advect diff problem
* \param[in] A A pointer to a PETSc Matrix which resembles the application of the system matrix resulting from the timestepping method
* \param[in] vec_global_index An indexing scheme which describes the relation between the global index in a PETSc Vector and the process-local degrees of freedom */
void                t8dg_precon_initialize_preconditioner (PC * pc, int selector, t8dg_precon_general_preconditioner_t * general_precon,
                                                           void *problem, t8dg_dof_values_t ** problem_dofs,
                                                           t8dg_time_derivation_matrix_application_t time_derivative, Mat * A,
                                                           PetscInt * vec_global_index);

/** This function desctroys a preconditioner and frees allocated memory which was used during the preconditioning routine 
* \param[in] selector An integer describing which preconditioner was used
* \param[in] general_precon A pointer to the preconditioner which ought to be destroyed
*/
void                t8dg_precon_destroy_preconditioner (PC * pc, int selector, t8dg_precon_general_preconditioner_t * general_precon);

/* Initializes the without-preconditioning option */
void                t8dg_precon_init_without_preconditioning (PC * pc);

/* Initializes a block jacobi preconditioner */
void                t8dg_precon_init_jacobi (void *problem, t8dg_dof_values_t ** problem_dofs, Mat * A, PC * pc,
                                             t8dg_block_preconditioner_ctx_t ** precon_jacobi_ctx, PetscInt * vec_global_index,
                                             int selector);

/* Destroys the allocated memory from the block jacobi preconditioner */
void                t8dg_precon_destroy_block_preconditioner (t8dg_block_preconditioner_ctx_t * ctx);

/** Initializes a two level multigrid preconditioner and calls all subroutines, therefore, some parameters need to be passed to the function which will be initilialized right here
* \param[in] problem A void pointer to the initial advection diffusion problem
* \param[in] problem_dofs A pointer to a pointer to t8dg_dof_values_t which holds the current degrees of freedom regarding the \a problem 
* \param[in] time_derivative A function implmenting the matrix-free application to the degrees of freedom \a *problem_dofs describing the time derivation of these coefficients
* \param[in] A_fine A pointer to a PETSc Matrix resembling the system matrix of an implicit timestepping method; the matrix application of the fine level
* \param[in, out] mg_pc A pointer to the preconditioner which preconditions the initial/fine level linear system; it preconditions \a A_fine; This preconditioner will be set up by this function
* \param[in] coarsening_func A \a t8dg_corase_lvl_adapt_func_t function describing the way the coarse level gets constructed/adapted from the initial/fine level problem
* \param[in, out] smoother A pointer to a PETSc KSP resembling the pre- and post-smoothing iterations which should be performed on the fine level; This KSP will be set up by this function
* \param[in, out] smoother_pc A pointer to the preconditioning context of the \a smoother; This preconditioner will be set up by this function
* \param[in, out] Restriction A pointer to a PETSc Matrix resembling the interpolation/restriction from the fine level onto the coarse level; This Matrix will be set up by this function
* \param[in, out] A_coarse A pointer to a PETSc Matrix describing the system matrix of the coarse level problem; it resembles the restriction of the fine level linear system, given \a A_fine; This Matrix will be set up by this function 
* \param[in, out] coarse_solver A pointer to a PETSc KSP which will calculate the solution of the coarse level problem, given by the coarse level matrix application \a A_coarse; This KSP will be set up by this function
* \param[in, out] coarse_pc A pointer to the preconditioning context of \a coarse_solver; This preconditioner will be set up by this function 
* \param[in, out] Prolongation A pointer to a PETSc Matrix resembling the interpolation/prolongation from the coarse level onto the fine level; This Matrix will be set up by this function
* \param[in, out] coarse_lvl A pointer to \a t8dg_mg_coarse_lvl_t context; this struct will hold the constructed coarse level mesh; The members will be filled and set up by this function
* \param[in, out] res_prol_ctx A pointer to \a t8dg_mg_interpolating_ctx_t context; This struct will hold information needed by the \a Restriction and \æ Prolongation routines; The members will be filled and set up by this function
* \param[in, out] cmat_ctx A pointer to \a t8dg_coarse_matrix_ctx_t context; This struct will hold the information needed by the matrix-application \a A_coarse of the coarse level problem; The members will be filled and set up by this function
* \param[in] vec_global_index A pointer to an indexing scheme describing the global position within the PETSc vectors of the fine level of the local degrees of freedom
* */
void                t8dg_precon_init_two_level_mg (void *problem, t8dg_dof_values_t ** problem_dofs,
                                                   t8dg_time_derivation_matrix_application_t time_derivative, Mat * A_fine, PC * mg_pc,
                                                   t8dg_corase_lvl_adapt_func_t coarsening_func, KSP * smoother, PC * smoother_pc,
                                                   Mat * Restriction, Mat * A_coarse, KSP * coarse_solver, PC * coarse_pc,
                                                   Mat * Prolongation, t8dg_mg_coarse_lvl_t * coarse_lvl,
                                                   t8dg_mg_interpolating_ctx_t * res_prol_ctx, t8dg_coarse_matrix_ctx_t * cmat_ctx,
                                                   PetscInt * vec_global_index);

/** This function destroys the allocated data needed by the before called initilization of the two level multigrid preconditioner \a t8dg_precon_init_two_level_mg()
* \param[in] smoother A pointer to the PETSc KSP resembling the pre- and post-smoothing iterations
* \param[in] Restriction A pointer to the PETSc Matrix resembling the interpolation/restriction from the fine level onto the coarse level
* \param[in] A_coarse A pointer to a PETSc Matrix describing the system matrix of the coarse level problem
* \param[in] coarse_solver A pointer to the PETSc KSP which solved the coarse level problem, given by the coarse level matrix application \a A_coarse
* \param[in] Prolongation A pointer the a PETSc Matrix resembling the interpolation/prolongation from the coarse level onto the fine level
* \param[in] coarse_lvl A pointer the \a t8dg_mg_coarse_lvl_t context which hold the coarse level mesh of the multigrid preconditioner
* \param[in] cmat_ctx A pointer to the \a t8dg_coarse_matrix_ctx_t context context which hold the information needed by the coarse matrix application \a A_coarse
 */
void                t8dg_precon_destroy_two_level_mg (KSP * smoother, Mat * Restriction, Mat * A_coarse, KSP * coarse_solver,
                                                      Mat * Prolongation, t8dg_mg_coarse_lvl_t * coarse_lvl,
                                                      t8dg_coarse_matrix_ctx_t * cmat_ctx);

/** A function that writes a t8dg_dof_values_t to a PETSc Vector
* \param[in] dofs A pointer to the degrees of freedom whose element_dofs will be copied into \a p_vec
* \param[in, out] p_vec A pointer to a PETSc Vector which entries will be filled 
* \param[in] indexing A global indexing scheme which provides the corresponding global index in the PETSc Vector to an entry in the (local) t8dg_dof_values_t \a dofs
*/
void                t8dg_precon_write_dof_to_vec (t8dg_dof_values_t * dofs, Vec * p_vec, PetscInt * indexing, size_t num_dofs);

/** A function that writes a PETSc Vector to a t8dg_dof_values_t 
* \param[in] p_vec A pointer to a PETSc Vector which entries will be copied
* \param[in, out] dofs A pointer to the degrees of freedom which will be filled with the entries of the PETSc Vector, specifically, the element_dofs of \a dofs will be overwritten with these entries 
*/
void                t8dg_precon_write_vec_to_dof (Vec * p_vec, t8dg_dof_values_t * dofs, size_t num_dofs);

/** A function that sets up the components needed in order to perform a two level multigrid V-Cycle 
* \param[in] problem A pointer to the initial advection diffusion problem
* \param[in] pdof_array A pointer to a pointer to t8dg_dof_values_t which holds the current degrees of freedom regarding the \a problem 
* \param[in] time_derivative A function implmenting the matrix-free application to the degrees of freedom \a *pdof_array describing the time derivation of these coefficients
* \param[in, out] coarse_lvl_mesh A pointer to a \a t8dg_mg_coarse_lvl_t context which members will be filled by this function
* \param[in] coarsening_func An adaption function describing how the coarse level mesh will be extracted from the initial forest of the \a problem
* \param[in, out] res_prol_ctx A pointer to a \a t8dg_mg_interpolating_ctx_t context which members will be filled by this function
* \param[in, out] cmat_ctx A pointer to a \a t8dg_mg_coarse_lvl_t context which members will be filled with the corresponding coarse level information
* \param[in] fine_forest_indexing A global indexing scheme for the entries in a PETSc Vector regarding the fine/intial level mesh
*/
void
 
 
 
 
 
 
 
 t8dg_mg_set_up_two_lvl_precon (t8dg_linear_advection_diffusion_problem_t * problem, t8dg_dof_values_t ** pdof_array,
                                t8dg_time_matrix_application time_derivative, t8dg_mg_coarse_lvl_t * coarse_lvl_mesh,
                                t8dg_corase_lvl_adapt_func_t coarsening_func, t8dg_mg_interpolating_ctx_t * res_prol_ctx,
                                t8dg_coarse_matrix_ctx_t * cmat_ctx, PetscInt * fine_forest_indexing);

/** A function that constructs the coarse level for the multigrid preconditioner; this function gets called inside \a t8dg_mg_set_up_two_lvl_precon() and is not needed to be called explicitly;
* The receiving parameters just get passed on by \a t8dg_mg_set_up_two_lvl_precon()
* (this function may be changed to static)
* \param[in] problem A pointer to the initial advection diffusion problem
* \param[in] pdof_array A pointer to a pointer to t8dg_dof_values_t which holds the current degrees of freedom regarding the \a problem 
* \param[in, out] coarse_lvl_mesh A pointer to a \a t8dg_mg_coarse_lvl_t context which members will be filled by this function
* \param[in] coarsening_func An adaption function describing how the coarse level mesh will be extracted from the initial forest of the \a problem
*/
void
 
 
 
 
 
 
 
 t8dg_mg_construct_coarse_lvl (t8dg_linear_advection_diffusion_problem_t * problem, t8dg_dof_values_t ** pdof_array,
                               t8dg_mg_coarse_lvl_t * coarse_lvl_mesh, t8dg_corase_lvl_adapt_func_t coarsening_func);

/** This functions creates a matrix-free coarse level application needed by the GMRES solver on the coarse level within the multigrid preconditioner
* \param[in, out] A_coarse A pointer to a PETSc Matrix which resembles the application to a PETSc Vector on the coarse level within the multigrid preconditioner; this matrix is created and assigned with the right properties by this function
* \param[in] coarse_lvl A pointer to a pointer to t8dg_dof_values_t which holds the current degrees of freesom regarding the \a problem 
* \param[in] cmat_ctx A pointer to \a t8dg_coarse_matrix_ctx_t context which holds the information needed for the creation of \a A_coarse
*/
void                t8dg_mg_create_coarse_lvl_matrix (Mat * A_coarse, t8dg_mg_coarse_lvl_t * coarse_lvl,
                                                      t8dg_coarse_matrix_ctx_t * cmat_ctx);

/** This functions creates a matrix-free routine to restrict the problem of the fine level onto the coarse level within the multigrid preconditioner
* \param[in, out] Restrcition A pointer to a PETSc Matrix which is created by this function and is assigned with properties and application to restrict a PETSc Vector of the fine level onto the coarse level in a matrix-free fashion
* \param[in] res_prol_ctx A pointer to \a t8dg_mg_interpolating_ctx_t context which provides information needed to create the \a Restriction matrix
*/
void                t8dg_mg_create_restriction_matrix (Mat * Restriction, t8dg_mg_interpolating_ctx_t * res_prol_ctx);

/** This functions creates a matrix-free routine to prolongate the problem of the coarse level onto the fine level within the multigrid preconditioner
* \param[in, out] Prolongation A pointer to a PETSc Matrix which is created by this function and is assigned with properties and application to prolongate a PETSc Vector of the coarse level onto the fine level in a matrix-free fashion
* \param[in] res_prol_ctx A pointer to \a t8dg_mg_interpolating_ctx_t context which provides information needed to create the \a Restriction matrix
*/
void                t8dg_mg_create_prolongation_matrix (Mat * Prolongation, t8dg_mg_interpolating_ctx_t * res_prol_ctx);

/** Thus function updates the preconditioner within the DIRK methods due to the varying coefficients of the different stages 
* \param[in] selector An integer describing which preconditioner was selected 
* \param[in, out] preconditioner A pointer to the preconditioner which has to be updated 
* \param[in] stage_related_a_coefficient The current a_coefficient (speaking of a Butcher Tableau) of the stage */
void                t8dg_precon_dirk_update_preconditioner (int selector, t8dg_precon_general_preconditioner_t * preconditioner,
                                                            double stage_related_a_coefficient, double stage_current_time);
double              t8dg_precon_get_setup_time (t8dg_precon_general_preconditioner_t * preconditioner);
#endif

T8DG_EXTERN_C_END ();

#endif /* SRC_T8DG_PRECONDITIONER_H_ */
