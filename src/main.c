/*
 * solver.c
 *
 *  Created on: Mar 11, 2020
 *      Author: lukas
 */
#include <t8.h>
#include <t8_cmesh.h>
#include <t8_vec.h>

#include "t8dg.h"
#include "t8dg_advect_problem.h"
#include "t8dg_coarse_geometry.h"
#include "t8dg_flux.h"
#include "t8dg_flux_implementation.h"
#include "t8dg_common.h"

#include <sc_options.h>
#include <sc.h>

#include <../example/common/t8_example_common.h>

static              t8dg_scalar_function_3d_time_fn
t8dg_choose_initial_cond_fn (int initial_cond_arg)
{
  switch (initial_cond_arg) {
  case (0):
    return t8_scalar3d_constant_one;
  case (1):
    return t8dg_scalar1d_hat_function;
  case (2):
    return t8_scalar3d_step_function;
  case (3):
    return t8_scalar3d_sinx;
  case (4):
    return t8dg_scalar3d_norm_function;
  case (5):
    return t8dg_scalar2d_hat_function;
  case (6):
    return t8dg_scalar2d_step_function;
  case (7):
    return t8dg_scalar2d_triangle_step_function;
  default:
    return t8_scalar3d_constant_zero;
  }
}

static              t8_cmesh_t
t8dg_choose_cmesh (int icmesh, sc_MPI_Comm comm)
{
  t8_cmesh_t          cmesh;
  switch (icmesh) {
  case 0:
    cmesh = t8_cmesh_new_hypercube (T8_ECLASS_LINE, comm, 0, 0, 1);
    break;
  case 1:
    cmesh = t8_cmesh_new_periodic_line_more_trees (comm);
    break;
  case 2:
    cmesh = t8dg_cmesh_new_periodic_diagonal_line_more_trees (comm);
    break;
  case 3:
    cmesh = t8_cmesh_new_hypercube (T8_ECLASS_QUAD, comm, 0, 0, 1);
    break;
  case 4:
    cmesh = t8dg_cmesh_new_square_more_trees_different_size (comm);
    break;
  case 5:
    cmesh = t8dg_cmesh_new_square_moebius (comm);
    break;
  case 6:
    cmesh = t8dg_cmesh_new_half_moebius_more_trees (comm);
    break;
  case 7:
    cmesh = t8dg_cmesh_new_square_tilted (comm);
    break;

  default:
    T8DG_ABORT ("Not yet implemented");
  }
  return cmesh;
}

static int
t8dg_check_options (int icmesh, int initial_cond_arg,
                    int uniform_level, int refinement_levels,
                    int number_LGL_points, double start_time, double end_time, double cfl, int time_order, int vtk_freq, int adapt_freq,
                    int adapt_arg)
{
  if (!(icmesh >= 0 && icmesh <= 7))
    return 0;
  if (!(initial_cond_arg >= 0 && initial_cond_arg <= 7))
    return 0;
  if (!(uniform_level >= 0 && uniform_level <= 30))
    return 0;
  if (!(refinement_levels >= 0 && uniform_level + refinement_levels <= 30))
    return 0;
  if (!(number_LGL_points >= 1 && number_LGL_points <= MAX_LGL_NUMBER))
    return 0;
  if (!(start_time < end_time))
    return 0;
  if (!(cfl > 0 && cfl <= 1))
    return 0;
  if (!(time_order >= 1 && time_order <= 4))
    return 0;
  if (!(vtk_freq >= 0))
    return 0;
  if (!(adapt_freq >= 0))
    return 0;
  if (!(adapt_arg >= 0 && adapt_arg <= 1))
    return 0;
  return 1;
}

t8_forest_adapt_t
t8dg_choose_adapt_fn (int adapt_arg)
{
  switch (adapt_arg) {
  case 0:
    return t8dg_advect_mass_adapt;
    break;
  case 1:
    return t8dg_advect_gradient_adapt;
    break;

  default:
    T8DG_ABORT ("Wrong adapt fn arg");
    break;
  }
}

/*TODO: which init function creates what, outsource problem description*/
static t8dg_linear_advection_problem_t *
t8dg_advect_problem_init_linear_geometry (int icmesh,
                                          int initial_cond_arg,
                                          double flow_speed,
                                          int uniform_level, int max_level,
                                          int number_LGL_points,
                                          double start_time, double end_time, double cfl, int time_order, int adapt_arg, sc_MPI_Comm comm)
{
  t8dg_scalar_function_3d_time_fn u_initial;
  t8_cmesh_t          cmesh;
  t8dg_coarse_geometry_t *coarse_geometry;

  t8dg_timestepping_data_t *time_data;
  double             *first_tree_vertices;
  int                 dim;
  t8_forest_adapt_t   adapt_fn;

  t8dg_linear_advection_problem_description_t *description;
  description = T8DG_ALLOC (t8dg_linear_advection_problem_description_t, 1);
  description->flux_data = T8_ALLOC (t8dg_linear_flux3D_constant_flux_data_t, 1);       /*TODO: free */
  double             *flow_direction = ((t8dg_linear_flux3D_constant_flux_data_t *) description->flux_data)->flow_direction;
  double             *p_flow_speed = &((t8dg_linear_flux3D_constant_flux_data_t *) description->flux_data)->flow_velocity;

  adapt_fn = t8dg_choose_adapt_fn (adapt_arg);

  if (icmesh <= 2) {
    coarse_geometry = t8dg_coarse_geometry_new_1D_linear ();
    cmesh = t8dg_choose_cmesh (icmesh, comm);
    dim = 1;
    first_tree_vertices = t8_cmesh_get_tree_vertices (cmesh, 0);

    t8_vec_axpyz (first_tree_vertices, first_tree_vertices + 3, flow_direction, -1);
    *p_flow_speed = flow_speed;
    description->velocity_field = t8dg_linear_flux3D_constant_flux_fn;

  }
  else {
    coarse_geometry = t8dg_coarse_geometry_new_2D_linear ();
    cmesh = t8dg_choose_cmesh (icmesh, comm);
    dim = 2;
    flow_direction[0] = 1;
    flow_direction[1] = 0;
    flow_direction[2] = 0;
    *p_flow_speed = flow_speed;
    description->velocity_field = t8dg_linear_flux3D_constant_flux_fn;
  }

  description->initial_condition_fn = t8dg_choose_initial_cond_fn (initial_cond_arg);

  time_data = t8dg_timestepping_data_new (time_order, start_time, end_time, cfl);

  return t8dg_advect_problem_init (cmesh, coarse_geometry, dim, description,
                                   uniform_level, max_level, number_LGL_points, time_data, adapt_fn, comm);
}

/** Solves the linear advection problem on a linear geometry on a uniform grid
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

void
t8dg_advect_solve (int icmesh, int initial_cond_arg,
                   double flow_velocity, int uniform_level, int refinement_levels,
                   int number_LGL_points, double start_time,
                   double end_time, double cfl, int time_order, int vtk_freq, int adapt_freq, int adapt_arg, sc_MPI_Comm comm)
{
  t8dg_linear_advection_problem_t *problem;
  int                 step_number;
  double              total_time, io_time, solve_time;
  double              l2_error, l_inf_error;

  t8dg_debugf ("Start Advection Solve\n");

  total_time = -sc_MPI_Wtime ();

  problem = t8dg_advect_problem_init_linear_geometry (icmesh, initial_cond_arg, flow_velocity,
                                                      uniform_level, uniform_level + refinement_levels,
                                                      number_LGL_points, start_time, end_time, cfl, time_order, adapt_arg, comm);

  t8dg_advect_problem_accumulate_stat (problem, ADVECT_INIT, total_time + sc_MPI_Wtime ());

  t8dg_debugf ("Dof values beginning:\n");
  t8dg_advect_problem_printdof (problem);

  /*Timeloop with Rungekutta timestepping: */
  while (!t8dg_advect_problem_endtime_reached (problem)) {
    t8dg_advect_problem_set_time_step (problem);
    step_number = t8dg_advect_problem_get_stepnumber (problem); /*TODO: could also simply be in this loop */
    if (vtk_freq && step_number % vtk_freq == 0) {
      io_time = -sc_MPI_Wtime ();
      t8dg_advect_write_vtk (problem);
      io_time += sc_MPI_Wtime ();
      t8dg_advect_problem_accumulate_stat (problem, ADVECT_IO, io_time);
    }
    solve_time = -sc_MPI_Wtime ();
    t8dg_advect_problem_advance_timestep (problem);
    solve_time += sc_MPI_Wtime ();
    t8dg_advect_problem_accumulate_stat (problem, ADVECT_SOLVE, solve_time);

    if (adapt_freq && step_number % adapt_freq == adapt_freq - 1) {
      t8dg_advect_problem_adapt (problem, 1);
      t8dg_advect_problem_partition (problem, 1);
    }
  }

  t8dg_advect_write_vtk (problem);

  total_time += sc_MPI_Wtime ();
  t8dg_advect_problem_accumulate_stat (problem, ADVECT_TOTAL, total_time);

  l2_error = t8dg_advect_problem_l2_rel (problem);
  l_inf_error = t8dg_advect_problem_l_infty_rel (problem);

  t8dg_advect_problem_accumulate_stat (problem, ADVECT_ERROR_2, l2_error);
  t8dg_advect_problem_accumulate_stat (problem, ADVECT_ERROR_INF, l_inf_error);

  t8dg_advect_problem_compute_and_print_stats (problem);

  /*Current output */
  t8dg_debugf ("End Dof values:\n");
  t8dg_advect_problem_printdof (problem);

  t8dg_advect_problem_destroy (&problem);
  return;
}

int
main (int argc, char *argv[])
{
  int                 mpiret;
  sc_options_t       *opt;
  char                help[BUFSIZ];
  int                 parsed, helpme;

  int                 initial_cond_arg;
  int                 uniform_level, refinement_levels;
  int                 time_order;
  int                 number_LGL_points;
  int                 vtk_freq;
  int                 adapt_freq;
  int                 adapt_arg;
  int                 icmesh;
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
  t8dg_init (SC_LP_DEBUG);
  t8_init (SC_LP_DEBUG);
#else
  t8dg_init (SC_LP_ESSENTIAL);
  t8_init (SC_LP_ESSENTIAL);
#endif

  /* initialize command line argument parser */
  opt = sc_options_new (argv[0]);

  sc_options_add_switch (opt, 'h', "help", &helpme, "Display a short help message.");
  sc_options_add_int (opt, 'l', "level", &uniform_level, 3, "The uniform initial refinement level of the mesh. Default: 3");
  sc_options_add_int (opt, 'L', "LGL", &number_LGL_points, 2, "The number of LGL basis points/basisfunctions in 1D. Default: 2");
  sc_options_add_double (opt, 'c', "flow_velocity", &flow_velocity, 1.0, "The flow velocity. Default: 1.0");
  sc_options_add_double (opt, 'C', "CFL", &cfl, 1.0, "The CFL number used to determine the timestep. Default: 1.0");
  sc_options_add_int (opt, 'o', "time_order", &time_order, 2,
                      "The order used for the runge Kutta timestepping (1<= order <=4). Default: 2");
  sc_options_add_double (opt, 't', "start_time", &start_time, 0.0, "The start time of the solve. Default: 0.0");
  sc_options_add_double (opt, 'T', "end_time", &end_time, 1.0, "The end time of the solve. Default: 1.0");
  sc_options_add_int (opt, 'i', "initial_cond", &initial_cond_arg, 0, "Choose initial condition function. Default: 0\n"
                      "\t\t0: constant function\n" "\t\t1: 1D hat function\n" "\t\t2: 1D step function\n" "\t\t3: 1D sine function\n"
                      "\t\t4: norm\n" "\t\t5: 2D hat\n" "\t\t6: 2D circle step function\n" "\t\t7: 2D triangle step function\n");
  sc_options_add_int (opt, 'r', "ref_levels", &refinement_levels, 0, "The number of refinement levels(>=0). Default: 0");
  sc_options_add_int (opt, 'v', "vkt_freq", &vtk_freq, 1, "The number of steps until new vtk output. Default: 1\t (0 means no vtk)");
  sc_options_add_int (opt, 'a', "adapt_freq", &adapt_freq, 1, "The number of steps until adapt. Default: 1\t (0 means no adapt)");
  sc_options_add_int (opt, 'm', "cmesh", &icmesh, 0, "Choose cmesh. Default: 0\n" "\t\t0: line 1 tree\n" "\t\t1: line 3 trees\n"
                      "\t\t2: diagonal line more trees\n" "\t\t3: square\n" "\t\t4: square different size trees\n" "\t\t5: square moebius\n"
                      "\t\t6: moebius more tree\n" "\t\t7: parallelogram");
  sc_options_add_int (opt, 'A', "adapt_fn", &adapt_arg, 0,
                      "Choose Adapt Function. Default: 0\n" "\t\t0: mass-crit\n" "\t\t1: gradient-crit\n");

  parsed = sc_options_parse (t8dg_get_package_id (), SC_LP_ERROR, opt, argc, argv);
  if (helpme) {
    /* display help message and usage */
    t8dg_global_essentialf ("%s\n", help);
    sc_options_print_usage (t8dg_get_package_id (), SC_LP_ERROR, opt, NULL);
  }
  else if (parsed >= 0 && t8dg_check_options (icmesh, initial_cond_arg, uniform_level, refinement_levels, number_LGL_points,
                                              start_time, end_time, cfl, time_order, vtk_freq, adapt_freq, adapt_arg)) {

    /* Computation */
    t8dg_advect_solve (icmesh, initial_cond_arg, flow_velocity,
                       uniform_level, refinement_levels, number_LGL_points,
                       start_time, end_time, cfl, time_order, vtk_freq, adapt_freq, adapt_arg, sc_MPI_COMM_WORLD);
  }
  else {
    /* wrong usage */
    t8dg_global_productionf ("\n\tERROR:Wrong usage.\n\n");
    sc_options_print_usage (t8dg_get_package_id (), SC_LP_ERROR, opt, NULL);
  }

  sc_options_destroy (opt);
  sc_finalize ();
  mpiret = sc_MPI_Finalize ();
  SC_CHECK_MPI (mpiret);
  return 0;
}
