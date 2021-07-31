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
#include <petscdmda.h>
#include <petscdm.h>
#endif

T8DG_EXTERN_C_BEGIN ();

typedef struct t8dg_timestepping_data t8dg_timestepping_data_t;

typedef void        (*t8dg_time_matrix_application) (t8dg_dof_values_t * src_dof, t8dg_dof_values_t * dest_dof, const double t,
                                                     const void *application_data);

void                t8dg_timestepping_choose_impl_expl_method (t8dg_time_matrix_application time_derivative,
                                                               t8dg_timestepping_data_t * time_data, t8dg_dof_values_t ** pdof_array,
                                                               void *user_data);

/* Struct that provides a context for matrix-free preconditioning */
#if T8_WITH_PETSC
typedef struct
{
  Vec                 matrix_diagonal;
  PetscScalar         timestep;
  PetscScalar         current_b_coeff;
} t8dg_timestepping_precon_jacobi_ctx_t;
#endif

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
} t8dg_timestepping_impl_euler_ctx_t;
#endif

/* Struct that keeps the information needed by the matrix-free application of the DIRK methods */
#if T8_WITH_PETSC
typedef struct
{
  int                 num_local_dofs;
  double              timestep;
  double              next_time_point;
  double              dirk_current_a_coeff;
  PetscInt           *global_indexing;
  t8dg_time_matrix_application time_derivative_func;
  t8dg_dof_values_t  *future_local_dofs_derivative;
  t8dg_dof_values_t  *future_local_dofs_step;
  t8dg_dof_values_t  *current_local_dofs;
  t8dg_dof_values_t  *local_derivation_degrees;
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
                                   t8dg_timestepping_data_t * time_data, t8dg_dof_values_t ** pdof_array, void *user_data);

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
                         t8dg_timestepping_data_t * time_data, t8dg_dof_values_t ** pdof_array, void *user_data, int num_order_steps);

/** implements a single step of runge kutta with a-values in the butcher-tableau only on the first minor diagonal
 * The time derivative application is given by time_derivative
 */
void
 
 
 
 
 
 
 
 t8dg_timestepping_runge_kutta_step (t8dg_time_matrix_application time_derivative,
                                     t8dg_timestepping_data_t * time_data, t8dg_dof_values_t ** pdof_array, void *user_data);

void                t8dg_runge_kutta_fill_coefficients (int time_order, double **prk_a, double **prk_b, double **prk_c);

t8dg_timestepping_data_t *t8dg_timestepping_data_new_cfl (int time_order, double start_time, double end_time, double cfl,
                                                          int use_implicit_timestepping);

t8dg_timestepping_data_t *t8dg_timestepping_data_new_constant_timestep (int time_order, double start_time, double end_time, double delta_t,
                                                                        int use_implicit_timestepping);

void                t8dg_timestepping_data_destroy (t8dg_timestepping_data_t ** ptime_data);

double              t8dg_timestepping_data_get_current_time (t8dg_timestepping_data_t * time_data);

void                t8dg_timestepping_data_set_current_time (t8dg_timestepping_data_t * time_data, double current_time);

double              t8dg_timestepping_data_get_end_time (t8dg_timestepping_data_t * time_data);

double              t8dg_timestepping_data_get_cfl (t8dg_timestepping_data_t * time_data);

int                 t8dg_timestepping_data_get_time_order (t8dg_timestepping_data_t * time_data);

double              t8dg_timestepping_data_get_time_step (t8dg_timestepping_data_t * time_data);

void                t8dg_timestepping_data_set_time_step (t8dg_timestepping_data_t * time_data, double delta_t);

int                 t8dg_timestepping_data_is_endtime_reached (t8dg_timestepping_data_t * time_data);

int                 t8dg_timestepping_data_get_step_number (t8dg_timestepping_data_t * time_data);

void                t8dg_timestepping_data_increase_step_number (t8dg_timestepping_data_t * time_data);

double              t8dg_timestepping_data_get_time_left (t8dg_timestepping_data_t * time_data);

T8DG_EXTERN_C_END ();

#endif /* SRC_TIMESTEPPING_H_ */
