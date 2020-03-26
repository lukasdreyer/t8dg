/*
 * t8dg_advect.h
 *
 *  Created on: Mar 22, 2020
 *      Author: lukas
 */
/** @file t8dg_advect_problem.h */

#ifndef SRC_T8DG_ADVECT_H_
#define SRC_T8DG_ADVECT_H_

typedef struct t8dg_linear_advection_problem t8dg_linear_advection_problem_t;

t8dg_linear_advection_problem_t *t8dg_advect_problem_init_linear_1D (t8_cmesh_t cmesh,
                                                                     t8dg_scalar_function_3d_time_fn u_0,
                                                                     double flow_velocity, int level,
                                                                     int number_LGL_points, double start_time,
                                                                     double end_time, double cfl, int time_order, sc_MPI_Comm comm);

void                t8dg_advect_problem_init_elements_linear_1D (t8dg_linear_advection_problem_t * problem);

void                t8dg_advect_problem_destroy (t8dg_linear_advection_problem_t ** pproblem);

void                t8dg_advect_runge_kutta_step (t8dg_linear_advection_problem_t * problem);

int                 t8dg_advect_problem_endtime_reached (t8dg_linear_advection_problem_t * problem);

void                t8dg_advect_problem_printdof (t8dg_linear_advection_problem_t * problem);

#endif /* SRC_T8DG_ADVECT_H_ */
