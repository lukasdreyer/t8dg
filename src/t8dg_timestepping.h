/*
 * timestepping.h
 *
 *  Created on: Mar 14, 2020
 *      Author: lukas
 */

/** @file t8dg_timestepping.h */
#ifndef SRC_TIMESTEPPING_H_
#define SRC_TIMESTEPPING_H_
#include <sc_containers.h>
#include "t8dg.h"
#include "t8dg_dof.h"

#if T8_WITH_PETSC
#include <petscksp.h>
#endif

T8DG_EXTERN_C_BEGIN ();

typedef struct t8dg_timestepping_data t8dg_timestepping_data_t;

typedef void        (*t8dg_time_matrix_application) (t8dg_dof_values_t * src_dof, t8dg_dof_values_t * dest_dof, const double t,
                                                     const void *application_data);

/** This function gets called by \a t8dg_advect_diff_advance_timestep() method and choose either an explicit or an implicit Runge Kutta timestepping method based on the input via the command line; default: is explicit RKV.
* The information which timestepping method was choosen lies within the \a tima_data.
* \param[in] time_derivative A function pointer to a function describing the problem dependent time derivation of the coefficients (degrees of freedom)
* \param[in] time_data A pointer to \a t8dg_timestepping_data_t holding all information regarding time and time-stepping
* \param[in, out] pdof_array A pointer to a pointer to the initial degrees of freedom of the advect diff problem; These degrees of freedom get filled by the timestepping methods with the approximation of the next time step
* \param[in] user_data A void pointer to the initial advect diff problem
*/
void                t8dg_timestepping_choose_impl_expl_method (t8dg_time_matrix_application time_derivative,
                                                               t8dg_timestepping_data_t * time_data, t8dg_dof_values_t ** pdof_array,
                                                               void *user_data);

/* Struct which keeps information regarding the matrix-free application of system matrix resulting from the implicit Euler-Method */
#if T8_WITH_PETSC
typedef struct
{
  PetscScalar         timestep;
  PetscInt           *global_indexing;
  double              next_time_point;
  size_t              num_local_dofs;
  t8dg_time_matrix_application time_derivative_func;
  t8dg_dof_values_t **future_local_dofs;
  t8dg_dof_values_t  *future_local_dofs_derivative;
  t8dg_timestepping_data_t *time_data;
  void               *user_data;
  double              application_time;
  int                 mf_app_count;
} t8dg_timestepping_impl_euler_ctx_t;
#endif

/* Struct that keeps the information needed by the matrix-free application of the DIRK methods */
#if T8_WITH_PETSC
typedef struct
{
  size_t              num_local_dofs;
  double              timestep;
  double              next_time_point;
  double              dirk_current_a_coeff;
  PetscInt           *global_indexing;
  t8dg_time_matrix_application time_derivative_func;
  t8dg_dof_values_t  *future_local_dofs_derivative;
  t8dg_dof_values_t  *future_local_dofs_step;
  t8dg_dof_values_t  *current_local_dofs;
  t8dg_dof_values_t **local_derivation_degrees;
  t8dg_timestepping_data_t *time_data;
  void               *user_data;
} t8dg_timestepping_dirk_ctx_t;
#endif

/** Advances a timestep using the implicit Euler method 
* \param[in] time_derivative The function that describes the derivation of the dofs of the current problem which has to be solved (u_{t} = R(u, t))
* \param[in] time_data       Keeps the information regarding the time, e.g. current time, timetep
* \param[in] pdof_array      Pointer to the dofs of the given problem
* \param[in] user_data       Contains the complete information regarding the problem, in fact it is a pointer to the initial problem 
* The former dofs (from the last time step) \a pdof_array will be directly exchanged with the calculated dofs of the next time step
*/
void
 
 
 
 
 
 
 
 t8dg_timestepping_implicit_euler (t8dg_time_matrix_application time_derivative,
                                   t8dg_timestepping_data_t * time_data, t8dg_dof_values_t ** pdof_array, void *user_data,
                                   int preconditioner_selection);

/** Initializes the apllication context needed by the matrix-free application of the system matrix resulting from the implicit euler method 
* \param[in] time_derivative The function describing the time derivative of the coefficients in the degrees of freedom
* \param[in] time_data A pointer to the time_data of the advect diff problem
* \param[in] pdof_array A pointer to a pointer to the degrees of freedom of the advect diff problem
* \param[in] user_data A void pointer to the advect diff problem
* \param[in, out] appctx A pointer to \a t8dg_timestepping_impl_euler_ctx_t context whose members will be filled by this function
*/
#if T8_WITH_PETSC
void
 
 
 
 
 
 
 
 t8dg_timestepping_init_impl_euler_appctx (t8dg_time_matrix_application time_derivative, t8dg_timestepping_data_t * time_data,
                                           t8dg_dof_values_t ** pdof_array, void *user_data, t8dg_timestepping_impl_euler_ctx_t * appctx);
#endif

/** Destroys the allocated space needed by the matrix application of the implicit euler method 
* \param[in] appctx The appctx which was used by the matrix application of the implicit euler method
*/
#if T8_WITH_PETSC
void                t8dg_timestepping_destroy_impl_euler_appctx (t8dg_timestepping_impl_euler_ctx_t * appctx);
#endif

/** Advances a timestep using a DIRK method  
* \param[in] time_derivative The function that describes the derivation of the dofs of the current problem which has to be solved (u_{t} = R(u, t))
* \param[in] time_data       Keeps the information regarding the time, e.g. current time, timetep
* \param[in] pdof_array      Pointer to the dofs of the given problem
* \param[in] user_data       Contains the complete information regarding the problem, in fact it is a pointer to the initial problem 
* \param[in] num_order_steps The of order/stages of the DIRK methods
* The former dofs (from the last time step) \a pdof_array will be directly exchanged with the calculated dofs of the next time step
* The only DIRK methods implemented are DIRK(2,2) of order 2 using 2 steps and DIRK(3,3) of order 3 using 3 steps. Therefore, \a num_order_stages got to equal 2 or 3
*/
void
 
 
 
 
 
 
 
 t8dg_timestepping_dirk (t8dg_time_matrix_application time_derivative,
                         t8dg_timestepping_data_t * time_data, t8dg_dof_values_t ** pdof_array, void *user_data, int num_order_steps,
                         int preconditioner_selection);

#if T8_WITH_PETSC
/** Allocates and fills members of t8dg_timestepping_dirk_ctx_t which are needed by the application of the DIRK methods 
* \param[in] time_derivative The function that describes the derivation of the dofs of the current problem which has to be solved (u_{t} = R(u, t))
* \param[in] time_data       Keeps the information regarding the time, e.g. current time, timetep
* \param[in] pdof_array      Pointer to the dofs of the given problem
* \param[in] user_data       Contains the complete information regarding the problem, in fact it is a pointer to the initial problem 
* \param[in, out] appctx     The application context of the DIRK methods; this function allocates and fills the members of the context needed in the DIRK routines
*/

void
 
 
 
 
 
 
 
 t8dg_timestepping_init_dirk_appctx (t8dg_time_matrix_application time_derivative, t8dg_timestepping_data_t * time_data,
                                     t8dg_dof_values_t ** pdof_array, void *user_data, t8dg_timestepping_dirk_ctx_t * appctx);
#endif

#if T8_WITH_PETSC
/** Deallocates/Destroys the application context needed by the application of the DIRK methods
* \param[in] appctx The previous initialized \a t8dg_timestepping_dirk_ctx_t (by calling \a t8dg_timestepping_init_dirk_appctx()) which is going to be destroyed
*/
void                t8dg_timestepping_destroy_dirk_appctx (t8dg_timestepping_dirk_ctx_t * appctx);
#endif

/** implements a single step of runge kutta with a-values in the butcher-tableau only on the first minor diagonal
 * The time derivative application is given by time_derivative
 */
void
 
 
 
 
 
 
 
 t8dg_timestepping_runge_kutta_step (t8dg_time_matrix_application time_derivative,
                                     t8dg_timestepping_data_t * time_data, t8dg_dof_values_t ** pdof_array, void *user_data);

void                t8dg_runge_kutta_fill_coefficients (int time_order, double **prk_a, double **prk_b, double **prk_c);

t8dg_timestepping_data_t *t8dg_timestepping_data_new_cfl (int time_order, double start_time, double end_time, double cfl,
                                                          int use_implicit_timestepping, int preconditioner_selection,
                                                          int multigrid_levels);

t8dg_timestepping_data_t *t8dg_timestepping_data_new_constant_timestep (int time_order, double start_time, double end_time, double delta_t,
                                                                        int use_implicit_timestepping, int preconditioner_selection,
                                                                        int multigrid_levels);

void                t8dg_timestepping_data_destroy (t8dg_timestepping_data_t ** ptime_data);

double              t8dg_timestepping_data_get_current_time (const t8dg_timestepping_data_t * time_data);

void                t8dg_timestepping_data_set_current_time (t8dg_timestepping_data_t * time_data, double current_time);

double              t8dg_timestepping_data_get_end_time (const t8dg_timestepping_data_t * time_data);

double              t8dg_timestepping_data_get_cfl (const t8dg_timestepping_data_t * time_data);

int                 t8dg_timestepping_data_get_time_order (const t8dg_timestepping_data_t * time_data);

double              t8dg_timestepping_data_get_time_step (const t8dg_timestepping_data_t * time_data);

void                t8dg_timestepping_data_set_time_step (t8dg_timestepping_data_t * time_data, double delta_t);

int                 t8dg_timestepping_data_is_endtime_reached (const t8dg_timestepping_data_t * time_data);

int                 t8dg_timestepping_data_get_step_number (const t8dg_timestepping_data_t * time_data);

void                t8dg_timestepping_data_increase_step_number (t8dg_timestepping_data_t * time_data);

double              t8dg_timestepping_data_get_time_left (const t8dg_timestepping_data_t * time_data);

int                 t8dg_timestepping_data_get_multigrid_levels (t8dg_timestepping_data_t * time_data);

T8DG_EXTERN_C_END ();

#endif /* SRC_TIMESTEPPING_H_ */
