/*
 * t8dg.c
 *
 *  Created on: Mar 11, 2020
 *      Author: lukas
 */
#include <sc_options.h>
#include <sc.h>

#include <t8_cmesh.h>
#include <t8.h>
#include <t8_vec.h>

#include "t8dg.h"
#include "t8dg_solver.h"

double
u_analytical (const double x[DIM3], double t)
{
//  return 1;
//  return (x[0] > .25 && x[0] < .75) ? 1.0 : 0.0;
  return 0.5 - fabs (0.5 - t8_vec_norm (x));
}

int
main (int argc, char *argv[])
{
  int                 mpiret;
  sc_options_t       *opt;
  char                help[BUFSIZ];
  int                 parsed, helpme;

  int                 level;
  int                 time_order;
  int                 number_LGL_points;
  double              flow_velocity;
  double              cfl;
  double              start_time;
  double              end_time;
  /* brief help message */

  /* long help message */

  snprintf (help, BUFSIZ, "This program solves the linear advection equation on the line. \n");
  mpiret = sc_MPI_Init (&argc, &argv);
  SC_CHECK_MPI (mpiret);

  sc_init (sc_MPI_COMM_WORLD, 1, 1, NULL, SC_LP_ESSENTIAL);
#ifdef T8_ENABLE_DEBUG
  t8_init (SC_LP_DEBUG);
#else
  t8_init (SC_LP_ESSENTIAL);
#endif

  /* initialize command line argument parser */
  opt = sc_options_new (argv[0]);

  sc_options_add_switch (opt, 'h', "help", &helpme, "Display a short help message.");
  sc_options_add_int (opt, 'l', "level", &level, 3, "The uniform refinement level of the mesh. Default: 3");
  sc_options_add_int (opt, 'L', "LGL", &number_LGL_points, 2, "The number of LGL basis points/basisfunctions in 1D. Default: 2");
  sc_options_add_double (opt, 'c', "flow_velocity", &flow_velocity, 1.0, "The flow velocity. Default: 1.0");
  sc_options_add_double (opt, 'C', "CFL", &cfl, 1.0, "The CFL number used to determine the timestep. Default: 1.0");
  sc_options_add_int (opt, 'o', "time_order", &time_order, 2,
                      "The order used for the runge Kutta timestepping (1<= order <=4). Default: 2");
  sc_options_add_double (opt, 't', "start_time", &start_time, 0.0, "The start time of the solve. Default: 0.0");
  sc_options_add_double (opt, 'T', "end_time", &end_time, 1.0, "The end time of the solve. Default: 1.0");

  parsed = sc_options_parse (t8_get_package_id (), SC_LP_ERROR, opt, argc, argv);
  if (helpme) {
    /* display help message and usage */
    t8_global_essentialf ("%s\n", help);
    sc_options_print_usage (t8_get_package_id (), SC_LP_ERROR, opt, NULL);
  }
  else if (parsed >= 0) {
    t8_cmesh_t          cmesh;

    cmesh = t8_cmesh_new_hypercube (T8_ECLASS_LINE, sc_MPI_COMM_WORLD, 0, 0, 1);

#if 0
    cmesh = t8_cmesh_new_periodic_line_more_trees (sc_MPI_COMM_WORLD);
#endif
    /* Computation */
    t8dg_advect_solve_1D (cmesh, u_analytical, flow_velocity,
                          level, number_LGL_points, start_time, end_time, cfl, time_order, sc_MPI_COMM_WORLD);
#if 0
    t8_cmesh_destroy (&cmesh);  /*t8_forest_unref takes care! */
#endif
  }
  else {
    /* wrong usage */
    t8_global_productionf ("\n\tERROR:Wrong usage.\n\n");
    sc_options_print_usage (t8_get_package_id (), SC_LP_ERROR, opt, NULL);
  }

  sc_options_destroy (opt);
  sc_finalize ();
  mpiret = sc_MPI_Finalize ();
  SC_CHECK_MPI (mpiret);
  return 0;
}
