/*
 * t8dg_advect.h
 *
 *  Created on: Mar 22, 2020
 *      Author: lukas
 */

#ifndef SRC_T8DG_ADVECT_H_
#define SRC_T8DG_ADVECT_H_

typedef struct t8dg_advect_problem_linear_1D t8dg_advect_problem_linear_1D_t;

t8dg_advect_problem_linear_1D_t *t8dg_advect_problem_init_linear_1D (t8_cmesh_t cmesh,
                                                                     t8dg_scalar_function_3d_time_fn u_0,
                                                                     double flow_velocity, int level,
                                                                     int number_LGL_points, double start_time,
                                                                     double end_time, double cfl, int time_order, sc_MPI_Comm comm);

void                t8dg_advect_problem_init_elements_linear_1D (t8dg_advect_problem_linear_1D_t * problem);

void                t8dg_advect_problem_destroy (t8dg_advect_problem_linear_1D_t ** pproblem);

void                t8dg_advect_evolve (t8dg_advect_problem_linear_1D_t * problem);

int                 t8dg_advect_problem_endtime_reached (t8dg_advect_problem_linear_1D_t * problem);

void                t8dg_advect_problem_printdof (t8dg_advect_problem_linear_1D_t * problem);

#endif /* SRC_T8DG_ADVECT_H_ */
