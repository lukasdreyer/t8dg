/*
 * solver.h
 *
 *  Created on: Mar 11, 2020
 *      Author: lukas
 */

#ifndef SRC_SOLVER_H_
#define SRC_SOLVER_H_

#include <t8_cmesh.h>
#include "t8dg.h"

void                t8dg_advect_solve_1D (t8_cmesh_t cmesh, t8dg_scalar_function_3d_time_fn u_initial,
                                          double flow_velocity, int level, int number_LGL_points, double start_time,
                                          double end_time, double cfl, int time_order, sc_MPI_Comm comm);

#endif /* SRC_SOLVER_H_ */
