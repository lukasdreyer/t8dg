/*
 * solver.h
 *
 *  Created on: Mar 11, 2020
 *      Author: lukas
 */

#ifndef SRC_SOLVER_H_
#define SRC_SOLVER_H_

#include <t8_cmesh.h>
#include <example/common/t8_example_common.h>

#include "t8dg_global.h"





void t8dg_1D_advect_solve (t8_cmesh_t cmesh, t8dg_scalar_function_MAX_DIMd_fn u_0, double flow_velocity,
			   int level, int number_LGL_points, sc_MPI_Comm comm);


#endif /* SRC_SOLVER_H_ */



