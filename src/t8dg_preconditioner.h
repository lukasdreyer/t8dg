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

/* Permitted number of maximal multigrid levels to use (initial mesh is included in this number -> a maximum of N-1 coarse levels may be used) */
#define T8DG_PRECON_MAX_MG_LVLS 7

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
  double              block_precon_application_time;
} t8dg_block_preconditioner_ctx_t;

#if 0
typedef struct
{
  int                 mg_coarse_order;
  size_t              num_fine_dofs;
  size_t              num_coarse_dofs;
  t8dg_linear_advection_diffusion_problem_t *problem;
  t8dg_values_t      *initial_dg_values;
  t8dg_values_t      *coarse_dg_values;
  t8dg_dof_values_t  *coarse_lvl_dofs;
  t8dg_dof_values_t  *coarse_lvl_dofs_derivation;
  t8dg_dof_values_t **problem_dofs;
  PetscInt           *coarse_lvl_indexing;
  PetscInt           *fine_lvl_indexing;
} t8dg_p_mg_mat_ctx_t;

typedef struct
{
  Mat                 Restriction;
  Mat                 Prolongation;
  Mat                 Coarse_Mat;
  Mat                 Smoothing_Mat;
  KSP                 smoother;
  KSP                 coarse_solver;
  PC                  smoother_pc;
  PC                  coarse_pc;
  t8dg_p_mg_mat_ctx_t *p_mg_mat_ctx;
} t8dg_p_mg_lvl_ctx_t;
#endif

/* Struct that holds information regarding the advect diff problem; it is needed by the multigrid preconditioner with multiple levels */
typedef struct
{
  void               *problem;
  t8dg_time_matrix_application time_derivative_func;
  double              current_point_in_time;
  double              timestep;
  double              current_coefficient;
  t8dg_dof_values_t **pdof_array;
  double              preconditioner_application_time_coarse_lvl;
  double              preconditioner_application_time_smoothing;
  double              preconditioner_application_time_interpolation;
} t8dg_mg_general_data_t;

/* Struct that keeps information needey by the solver on the coarsesrt level during multigrid preconditioning */
typedef struct
{
  size_t              num_local_dofs_coarsest_lvl;
  PetscInt           *indexing_coarsest_level;
  t8dg_dof_values_t **coarse_lvl_problem_dofs;
  t8dg_dof_values_t  *coarse_lvl_problem_dofs_derivation;
  t8dg_mg_general_data_t *mg_general_data;
} t8dg_mg_lvl_coarse_matrix_ctx_t;

/* Struct that keeps information needed by the restriction and prolongation operators during multigrid preconditioning */
typedef struct
{
  int                 current_lvl;
  int                 num_mg_levels;
  size_t              num_local_dofs[T8DG_PRECON_MAX_MG_LVLS];
  PetscInt           *indexing_scheme[T8DG_PRECON_MAX_MG_LVLS];
  t8_forest_t         coarse_level_forests[T8DG_PRECON_MAX_MG_LVLS];
  t8dg_dof_values_t  *dofs_lvl[T8DG_PRECON_MAX_MG_LVLS];
  t8dg_dof_values_t  *dofs_lvl_derivation[T8DG_PRECON_MAX_MG_LVLS - 2];
  t8dg_adapt_data_t  *adapt_data;
  t8dg_values_t      *dg_values;
  t8dg_local_values_t *local_values_lvl[T8DG_PRECON_MAX_MG_LVLS];
  t8dg_mortar_array_t *mortar_array_lvl[T8DG_PRECON_MAX_MG_LVLS];
  t8dg_mg_general_data_t *mg_general_data;
} t8dg_mg_lvl_interpolation_ctx_t;

/* Struct that keeps all information for multigrid preconditioning with multiple coarse levels */
typedef struct
{
  KSP                 coarse_solver;
  KSP                 smoothers[T8DG_PRECON_MAX_MG_LVLS - 1];
  PC                  coarse_pc;
  PC                  smoother_pcs[T8DG_PRECON_MAX_MG_LVLS - 1];
  Mat                 Restriction_Mats[T8DG_PRECON_MAX_MG_LVLS - 1];
  Mat                 Prolongation_Mats[T8DG_PRECON_MAX_MG_LVLS - 1];
  Mat                 Smoothing_Mats[T8DG_PRECON_MAX_MG_LVLS - 1];
  Mat                 A_coarse_Mat;
  t8dg_mg_lvl_interpolation_ctx_t *interpolation_ctx;
  t8dg_mg_lvl_coarse_matrix_ctx_t *coarse_matrix_ctx;
  t8dg_mg_general_data_t *mg_general_data;
  t8dg_block_preconditioner_ctx_t *precon_jacobi_ctx;
} t8dg_mg_levels_ctx_t;

/* Struct that fulfills a general purpose and holds Matrices, Contexts, etc. for every possible preconditioner */
typedef struct
{
  /* Multiple Level Multigrid preconditioner */
  t8dg_mg_levels_ctx_t *multiple_mg_lvls_ctx;

  /* p-Multigrid preconditioner */
  //t8dg_p_mg_lvl_ctx_t *p_mg_ctx;

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
PetscErrorCode      t8dg_precon_init_without_preconditioning (PC * pc);

/* Initializes a block jacobi preconditioner */
PetscErrorCode      t8dg_precon_init_jacobi (void *problem, t8dg_dof_values_t ** problem_dofs, Mat * A, PC * pc,
                                             t8dg_block_preconditioner_ctx_t ** precon_jacobi_ctx, PetscInt * vec_global_index,
                                             int selector);

/* Destroys the allocated memory from the block jacobi preconditioner */
PetscErrorCode      t8dg_precon_destroy_block_preconditioner (t8dg_block_preconditioner_ctx_t * ctx);

/** A function that writes a t8dg_dof_values_t to a PETSc Vector
* \param[in] dofs A pointer to the degrees of freedom whose element_dofs will be copied into \a p_vec
* \param[in, out] p_vec A pointer to a PETSc Vector which entries will be filled 
* \param[in] indexing A global indexing scheme which provides the corresponding global index in the PETSc Vector to an entry in the (local) t8dg_dof_values_t \a dofs
*/
PetscErrorCode      t8dg_precon_write_dof_to_vec (t8dg_dof_values_t * dofs, Vec * p_vec, PetscInt * indexing, size_t num_dofs);

/** A function that writes a PETSc Vector to a t8dg_dof_values_t 
* \param[in] p_vec A pointer to a PETSc Vector which entries will be copied
* \param[in, out] dofs A pointer to the degrees of freedom which will be filled with the entries of the PETSc Vector, specifically, the element_dofs of \a dofs will be overwritten with these entries 
*/
PetscErrorCode      t8dg_precon_write_vec_to_dof (Vec * p_vec, t8dg_dof_values_t * dofs, size_t num_dofs);

/** This function updates the preconditioner within the DIRK methods due to the varying coefficients of the different stages 
* \param[in] selector An integer describing which preconditioner was selected 
* \param[in, out] preconditioner A pointer to the preconditioner which has to be updated 
* \param[in] stage_related_a_coefficient The current a_coefficient (speaking of a Butcher Tableau) of the stage */
void                t8dg_precon_dirk_update_preconditioner (int selector, t8dg_precon_general_preconditioner_t * preconditioner,
                                                            double stage_related_a_coefficient, double stage_current_time);
double              t8dg_precon_get_setup_time (t8dg_precon_general_preconditioner_t * preconditioner);

/** This function initializes a mutligrid preconditioner with a given number of meshes (-> N-1 coarse level meshes) 
* \param[in] problem A void pointer to the initial advect diff problem
* \param[in] problem_dofs A pointer to a pointer to the initial dofs of the advect diff problem
* \param[in] time_derivative A function pointer describing the time-derivation of the dofs
* \param[in] A_fine A pointer to a PETSc Matrix resembling the global Matrix (on the intial mesh) resulting from a implicit time-stepping method (e.g. Implicit Euler)
* \param[in, out] mg_pc A pointer to the preconditioning context of the KSP solver which is used to slove global sytem dscribed by \a A_fine; This context is constructed by this method
* \param[in] vec_global_index A vector containing the ordering of the dofs of the global system (on the initial mesh)
* \param[in] num_mg_levels The number of levels to use within multigrid preconditioning (-> num_mg_levels -1 coarse levels are constructed)
* \param[in, out] mg_ctx A pointer to a pointer to the multigrid context, which contains the matrices, solvers, indexing_schemes, etc. This context and it's members are filled by this function
* \note The maximum number of permitted multigrid levels is set by the \a T8DG_PRECON_MAX_MG_LVLS macro
*/
PetscErrorCode
 
 
 t8dg_precon_init_multiple_level_mg (void *problem, t8dg_dof_values_t ** problem_dofs,
                                     t8dg_time_derivation_matrix_application_t time_derivative, Mat * A_fine, PC * mg_pc,
                                     PetscInt * vec_global_index, int num_mg_levels, t8dg_mg_levels_ctx_t ** mg_ctx);

/** This function set the operators required (interpolation, coarse level solver, etc.) up
* \param[in] problem A void pointer to the initial advect diff problem
* \param[in] pdof_array A pointer to a pointer to the initial dofs of the advect diff problem
* \param[in] time_derivative A function pointer describing the time-derivation of the dofs
* \param[in] vec_global_index A vector containing the ordering of the dofs of the global system (on the initial mesh)
* \param[in] num_mg_levels The number of levels to use within multigrid preconditioning (-> num_mg_levels -1 coarse levels are constructed)
* \param[in, out] mg_ctx A pointer to a pointer to the multigrid context, which contains the matrices, solvers, indexing_schemes, etc. Additional members of the context are filled by this function
*/
void
 
 
 
 
 
 
 
 t8dg_mg_set_up_multiple_lvl_precon (t8dg_linear_advection_diffusion_problem_t * problem, t8dg_dof_values_t ** pdof_array,
                                     t8dg_time_matrix_application time_derivative, PetscInt * initial_forest_indexing, int num_mg_levels,
                                     t8dg_mg_levels_ctx_t * mg_ctx);

/** This function constructs a new coarse mesh at the level \a num_lvl in the hierachy
* \param[in, out] mg_lvls A pointer to the multigrid context
* \param[in] coarsening_func The coarsening function which is used to construct the new forest
* \param[in] num_lvl the level of the new coarse mesh in the multigrid hierachy 
*/
PetscErrorCode
     t8dg_mg_construct_coarse_level_forest (t8dg_mg_levels_ctx_t * mg_lvls, t8dg_corase_lvl_adapt_func_t coarsening_func, int num_lvl);

/** This functions destroys the multigrid preconditioner and the space it allocated
* \param[in] mg_ctx The reference of the pointer to the multigrid context which is going to be destroyed
*/
PetscErrorCode      t8dg_precon_destroy_mg_levels (t8dg_mg_levels_ctx_t ** mg_ctx);

/** This function constructs the PETSc Matrix on the coarsest level within the multigrid preconditioner
* \param[in, out] mg_lvls A pointer to the multigrid context
*/
PetscErrorCode      t8dg_mg_set_up_coarse_lvl_matrix (t8dg_mg_levels_ctx_t * mg_lvls);

/** This function sets up the prolonagtion operators between the levels, interpolating from the coarse to the fine level
* \param[in, out] mg_lvls A pointer to the multigrid context
*/
PetscErrorCode      t8dg_mg_set_up_prolongation_operators (t8dg_mg_levels_ctx_t * mg_lvls);

/** This function sets up the restriction operators between the levels, interpolating from the fine to the coarse level
* \param[in, out] mg_lvls A pointer to the multigrid context
*/
PetscErrorCode      t8dg_mg_set_up_restriction_operators (t8dg_mg_levels_ctx_t * mg_lvls);

/** This function set up the smoothing operators(PETSc Mats, KSPs, ...) on each multigrid level (except on the coarsest)
* \param[in, out] mg_ctx A pointer to the multigrid context
* \param[in] A_fine A pointer to the fine/initial level matrix, because Smoothing on the finest level is in fact the application of the initial matrix
*/
PetscErrorCode      t8dg_mg_set_up_smoothing_operators (t8dg_mg_levels_ctx_t * mg_ctx, Mat * A_fine);

/** This function assigns the correct smoothers to each level
* \param[in, out] mg_pc A pointer to the preconditioning context of the KSP which solves the initial problem (resulting from an implicit time-stepping method)
* \param[in, out] mg_ctx A pointer to the multigrid context
*/
PetscErrorCode      t8dg_mg_initialize_smoothers (PC * mg_pc, t8dg_mg_levels_ctx_t ** mg_ctx);

/** This function assigns the correct prolongation operators between the levels
* \param[in, out] mg_pc A pointer to the preconditioning context of the KSP which solves the initial problem (resulting from an implicit time-stepping method)
* \param[in, out] mg_ctx A pointer to the multigrid context
*/
PetscErrorCode      t8dg_mg_initialize_prolongations (PC * mg_pc, t8dg_mg_levels_ctx_t ** mg_ctx);

/** This function assigns the correct restriction operators between the levels
* \param[in, out] mg_pc A pointer to the preconditioning context of the KSP which solves the initial problem (resulting from an implicit time-stepping method)
* \param[in, out] mg_ctx A pointer to the multigrid context
*/
PetscErrorCode      t8dg_mg_initialize_restrictions (PC * mg_pc, t8dg_mg_levels_ctx_t ** mg_ctx);

/* These following methods are currently not in use ---> p-multigrid is not yet possible */
#if 0
void
 
 
 
 
 
 
 
 t8dg_precon_init_p_mg (void *problem, t8dg_dof_values_t ** problem_dofs, t8dg_time_derivation_matrix_application_t time_derivative,
                        Mat * A_fine, PC * mg_pc, PetscInt * vec_global_index, int coarse_lvl_order, t8dg_p_mg_lvl_ctx_t ** mg_ctx);

void                t8dg_p_mg_set_up_smoothing_operator (t8dg_p_mg_lvl_ctx_t * mg_ctx, Mat * A_fine);

void                t8dg_p_mg_set_up_coarse_level (t8dg_p_mg_lvl_ctx_t * mg_ctx);

void                t8dg_p_mg_set_up_prolongation_operator (t8dg_p_mg_lvl_ctx_t * mg_ctx);

void                t8dg_p_mg_set_up_restriction_operator (t8dg_p_mg_lvl_ctx_t * mg_ctx);

void
 
 
 
 
 
 
 
 t8dg_p_mg_set_up_coarse_lvl (t8dg_linear_advection_diffusion_problem_t * problem, t8dg_dof_values_t ** problem_dofs, int coarse_lvl_order,
                              t8dg_p_mg_lvl_ctx_t * mg_ctx);

void                t8dg_precon_destroy_p_mg (t8dg_p_mg_lvl_ctx_t ** mg_ctx);
#endif

#endif

T8DG_EXTERN_C_END ();

#endif /* SRC_T8DG_PRECONDITIONER_H_ */
