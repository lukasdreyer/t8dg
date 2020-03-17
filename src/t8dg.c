/*
 * t8dg.c
 *
 *  Created on: Mar 11, 2020
 *      Author: lukas
 */
#include <sc_options.h>
#include <t8_cmesh.h>
#include <t8.h>

#include "t8dg.h"
#include "t8dg_solver.hxx"


double u_0(const double x[DIM3]){
  if(x[0]>=0.25 && x[0]<=0.75)return 1 - 4 * abs(x[0]);
  return 0;
}

int
main (int argc, char *argv[])
{
  int 			mpiret;
  sc_options_t		*opt;
  char			help[BUFSIZ];
  int			parsed, helpme;

  t8_eclass_t		eclass = T8_ECLASS_LINE;

  int			level,number_LGL_points;
  double		flow_velocity;


  /* brief help message */

  /* long help message */

  snprintf (help, BUFSIZ,
            "This program solves the linear advection equation on "
            "the line.\n");
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

  sc_options_add_switch (opt, 'h', "help", &helpme,
                         "Display a short help message.");
  sc_options_add_int (opt, 'l', "level", &level, 0,
                      "The uniform refinement level of the mesh.");
  sc_options_add_int (opt, 'L', "Number of LGL points", &number_LGL_points, 2,
                      "The number of LGL basis points/basisfunctions in 1D.");
  sc_options_add_double (opt, 'c', "flow_velocity", &flow_velocity, 1.0,
                      "The flow velocity.");


  parsed =
    sc_options_parse (t8_get_package_id (), SC_LP_ERROR, opt, argc, argv);
  if (helpme) {
    /* display help message and usage */
    t8_global_essentialf ("%s\n", help);
    sc_options_print_usage (t8_get_package_id (), SC_LP_ERROR, opt, NULL);
  }
  else if (parsed >= 0) {
    t8_cmesh_t			cmesh;

    cmesh = t8_cmesh_new_hypercube (eclass, sc_MPI_COMM_WORLD, 0, 0, 1);
    /* Computation */
    t8dg_1D_advect_solve (cmesh, u_0, flow_velocity,
			   level, number_LGL_points, sc_MPI_COMM_WORLD);
#if 0
    t8_cmesh_destroy(&cmesh);/*t8_forest_unref takes care! */
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
