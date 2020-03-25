/*
 * solver.h
 *
 *  Created on: Mar 11, 2020
 *      Author: lukas
 */
/** @file t8dg_solver.h */
#ifndef SRC_SOLVER_H_
#define SRC_SOLVER_H_

#include <t8_cmesh.h>
#include "t8dg.h"

/** Solves the (linear) 1D advection problem on a linear geometry on a uniform grid
 *
 * \param [in] cmesh            		The coarse mesh
 * \param [in] u_initial            		time-dependent initial function
 * \param [in] flow_velocity            	scalar flow velocity, TODO: timedependent vector field
 * \param [in] level            		uniform refinement level
 * \param [in] number_LGL_points            	in 1D the number of LGL quadrature vertices also used for functionbasis
 * \param [in] start_time            		start_time of the simulation
 * \param [in] end_time             		end_time of the simulation
 * \param [in] cfl            			cfl number used to determine the timestep delta_t
 * \param [in] time_order            		order of the runge-kutta timestepping
 * \param [in] comm            			MPI Communicator
 */

void                t8dg_advect_solve_1D (t8_cmesh_t cmesh, t8dg_scalar_function_3d_time_fn u_initial,
                                          double flow_velocity, int level, int number_LGL_points, double start_time,
                                          double end_time, double cfl, int time_order, sc_MPI_Comm comm);

#endif /* SRC_SOLVER_H_ */
