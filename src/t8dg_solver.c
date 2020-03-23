/*
 * solver.c
 *
 *  Created on: Mar 11, 2020
 *      Author: lukas
 */
#include <t8.h>
#include <example/common/t8_example_common.h>
#include <t8_forest.h>

#include "t8dg.h"
#include "t8dg_solver.h"
#include "t8dg_advect_problem.h"
/*
static void t8dg_flatten_jacobian_matrix(double *flat_array,t8dg_square_3D_matrix_t jacobian_matrix, int dim){
  int ixdim,iydim;
  for(ixdim=0; ixdim < dim; ixdim++){
    for(iydim = 0; iydim < dim ; iydim++){
      flat_array[ixdim * dim + iydim] = jacobian_matrix[ixdim][iydim];
    }
  }
}
*/



void t8dg_advect_solve_1D (t8_cmesh_t cmesh, t8dg_scalar_function_3d_time_fn u_initial, double flow_velocity,
			   int level, int number_LGL_points,
			   double start_time, double end_time, double cfl, int time_order,
			   sc_MPI_Comm comm)
{
  t8dg_advect_problem_linear_1D_t	*problem;

  t8_debugf("Start Advection Solve\n");

  problem = t8dg_advect_problem_init_linear_1D (cmesh, u_initial, flow_velocity,
  				   level, number_LGL_points,
				   start_time, end_time, cfl, time_order,
				   comm);
  t8dg_advect_problem_init_elements_linear_1D (problem);


  t8dg_advect_problem_printdof(problem);


  /*Timeloop with Rungekutta timestepping: */
  while(!t8dg_advect_problem_endtime_reached(problem)){
      t8dg_advect_evolve(problem);
  }

  t8dg_advect_problem_printdof(problem);

  t8dg_advect_problem_destroy(&problem);
  return;
}



