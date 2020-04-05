/*
 * t8dg_advect.h
 *
 *  Created on: Mar 22, 2020
 *      Author: lukas
 */
/** @file t8dg_advect_problem.h */

#ifndef SRC_T8DG_ADVECT_H_
#define SRC_T8DG_ADVECT_H_

#include <t8_cmesh.h>
#include <sc.h>
//#include "t8dg_timestepping.h"

typedef struct t8dg_linear_advection_problem t8dg_linear_advection_problem_t;

t8dg_linear_advection_problem_t *t8dg_advect_problem_init_linear_geometry_1D (t8_cmesh_t cmesh,
                                                                              t8dg_scalar_function_3d_time_fn u_0,
                                                                              double flow_velocity, int level,
                                                                              int number_LGL_points, double start_time,
                                                                              double end_time, double cfl, int time_order,
                                                                              sc_MPI_Comm comm);

void                t8dg_advect_problem_init_elements (t8dg_linear_advection_problem_t * problem);

void                t8dg_advect_problem_destroy (t8dg_linear_advection_problem_t ** pproblem);

int                 t8dg_advect_problem_endtime_reached (t8dg_linear_advection_problem_t * problem);

void                t8dg_advect_problem_printdof (t8dg_linear_advection_problem_t * problem);

void                t8dg_advect_write_vtk (t8dg_linear_advection_problem_t * problem);

void                t8dg_advect_problem_advance_timestep (t8dg_linear_advection_problem_t * problem);

//t8dg_timestepping_data_t *t8dg_advect_get_time_data (t8dg_linear_advection_problem_t * problem);

#endif /* SRC_T8DG_ADVECT_H_ */
