/*
 * solver.c
 *
 *  Created on: Mar 11, 2020
 *      Author: lukas
 */
#include <t8.h>
//#include <example/common/t8_example_common.h>
//#include <t8_forest.h>

#include "t8dg.h"
#include "t8dg_solver.h"
#include "t8dg_advect_problem.h"

void
t8dg_advect_solve_1D (t8_cmesh_t cmesh, t8dg_scalar_function_3d_time_fn u_initial,
                      double flow_velocity, int level, int number_LGL_points, double start_time,
                      double end_time, double cfl, int time_order, sc_MPI_Comm comm)
{
  t8dg_linear_advection_problem_t *problem;

  t8_debugf ("Start Advection Solve\n");

  problem = t8dg_advect_problem_init_linear_geometry_1D (cmesh, u_initial, flow_velocity,
                                                         level, number_LGL_points, start_time, end_time, cfl, time_order, comm);
  t8dg_advect_problem_init_elements (problem);

  /*Current output */
  t8_global_productionf ("Start Dof values:\n");
  t8dg_advect_problem_printdof (problem);

  t8dg_advect_write_vtk (problem);

  /*Timeloop with Rungekutta timestepping: */
  while (!t8dg_advect_problem_endtime_reached (problem)) {
    t8dg_timestepping_runge_kutta_step (t8dg_advect_time_derivative, t8dg_advect_get_time_data (problem),
                                        t8dg_advect_get_dof_sc_array_pointer (problem), problem);
    t8dg_advect_write_vtk (problem);
  }

  /*Current output */
  t8_global_productionf ("End Dof values:\n");
  t8dg_advect_problem_printdof (problem);

  t8dg_advect_problem_destroy (&problem);
  return;
}
