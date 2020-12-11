#include <t8.h>
#include <t8_forest.h>
#include <t8_vec.h>
#include <t8_schemes/t8_default_cxx.hxx>
#include <t8_forest/t8_forest_iterate.h>
#include <t8_forest/t8_forest_partition.h>

#include <sc_containers.h>
#include <sc_statistics.h>

#include "t8dg.h"
#include "t8dg_coarse_geometry.h"
#include "t8dg_advect_problem.h"
#include "t8dg_flux.h"
#include "t8dg_global_values.h"
#include "t8dg_local_values.h"
#include "t8dg_values.h"
#include "t8dg_quadrature.h"
#include "t8dg_functionbasis.h"
#include "t8dg_sc_array.h"
#include "t8dg_timestepping.h"
#include "t8dg_mortar.h"
#include "t8dg_geometry.h"
#include "t8dg_flux.h"
#include "t8dg_flux_implementation.h"
#include "t8dg_adapt.h"
#include "t8dg_output.h"

/* Enum for statistics. */
typedef enum
{
  ADVECT_ADAPT = 0,             /* adapt runtime */
  ADVECT_PARTITION,             /* partition runtime */
  ADVECT_PARTITION_DATA,        /* data partitioning runtime */
  ADVECT_BALANCE,               /* balance runtime */
  ADVECT_GHOST,                 /* ghost runtime */
  ADVECT_GHOST_EXCHANGE,        /* ghost exchange runtime */
  ADVECT_GHOST_WAIT,            /* ghost exchange waittime */
  ADVECT_REPLACE,               /* forest_iterate_replace runtime */
  ADVECT_IO,                    /* vtk runtime */
  ADVECT_INIT,                  /* initialization runtime */
  ADVECT_AMR,                   /* AMR runtime (adapt+partition+ghost+balance) including data exchange (partition/ghost) */
  ADVECT_SOLVE,                 /* solver runtime */
  ADVECT_TOTAL,                 /* overall runtime */

  ADVECT_PARTITION_PROCS,       /* number of processes sent to in partition */
  ADVECT_BALANCE_ROUNDS,        /* number of rounds in balance */
  ADVECT_GHOST_SENT,            /* number of ghosts sent to other processes */
  ADVECT_ELEM_AVG,              /* average global number of elements (per time step) */
  ADVECT_ERROR_INF,             /* l_infty error */
  ADVECT_ERROR_2,               /* L_2 error */
  ADVECT_VOL_LOSS,              /* The loss in volume in percent */

  ADVECT_NUM_STATS,             /* The number of statistics that we measure */
  ADVECT_NUM_TIME_STATS = ADVECT_TOTAL + 1      /* The number of time statistics that we only want to have counted once */
} advect_stats_t;

/* Names of statistics that we measure */
const char         *advect_stat_names[ADVECT_NUM_STATS] = {
  "adapt",
  "partition",
  "partition_data",
  "balance",
  "ghost",
  "ghost_exchange",
  "ghost_exchange_wait",
  "replace",
  "vtk_print",
  "init",
  "AMR",
  "solve",
  "total",

  "partition_procs_sent",
  "balance_rounds",
  "ghost_sent",
  "number_elements",

  "l_infty_error",
  "l_2_error",
  "mass_loss_[%]"
};

/** The container for all information needed to solve the advection problem.
 */
struct t8dg_linear_advection_problem
{
  int                 dim;      /**< Dimension of the submanifold */

  int                 uniform_refinement_level; /**< uniform refinement level */
  int                 maximum_refinement_level;

  t8_forest_t         forest;                                   /**< The t8_forest used to adapt the coarse mesh*/

  t8dg_linear_advection_problem_description_t *description;

  t8_forest_adapt_t   adapt_fn;

  t8dg_timestepping_data_t *time_data;

  t8dg_coarse_geometry_t *coarse_geometry;   /**< coarse geometry, that gives the geometry and Jacobian of the geometry from the coarse
						 reference element to the coarse image element*/
  int                 vtk_count;
  sc_MPI_Comm         comm; /**< MPI Communicator */
  sc_statinfo_t       stats[ADVECT_NUM_STATS]; /**< Runtimes and other statistics. */

/* The following sc_arrays contain data, that needs to be available for each processorlocal element. To facilitate partitioning, they are
 * saved in a linear array
 */
  t8dg_values_t      *dg_values;

  /* The dof_values get ghosted */
  t8dg_dof_values_t  *dof_values;              /**< The Value of u at the nodal basis vertices */
  t8dg_dof_values_t  *dof_values_adapt;

  double              total_time;       /* used for measuring the total time the struct is in existence */
};

void
t8dg_advect_problem_compute_and_print_stats (t8dg_linear_advection_problem_t * problem)
{
  sc_stats_compute (problem->comm, ADVECT_NUM_STATS, problem->stats);
  sc_stats_print (t8dg_get_package_id (), SC_LP_ESSENTIAL, ADVECT_NUM_STATS, problem->stats, 1, 1);
}

static void
t8dg_advect_problem_accumulate_stat (t8dg_linear_advection_problem_t * problem, advect_stats_t stat, double value)
{
  T8DG_ASSERT (stat >= 0 && stat < ADVECT_NUM_STATS);
  if (stat < ADVECT_NUM_TIME_STATS) {
    sc_stats_set1 (&problem->stats[stat], problem->stats[stat].sum_values + value, advect_stat_names[stat]);
  }
  else {
    sc_stats_accumulate (&problem->stats[stat], value);
  }
}

t8dg_linear_advection_problem_t *
t8dg_advect_problem_init (t8_cmesh_t cmesh,
                          t8dg_coarse_geometry_t * coarse_geometry,
                          int dim, t8dg_linear_advection_problem_description_t * description,
                          int uniform_level, int max_level, int number_LGL_points, t8dg_timestepping_data_t * time_data,
                          t8_forest_adapt_t adapt_fn, sc_MPI_Comm comm)
{
  t8dg_linear_advection_problem_t *problem;
  t8_scheme_cxx_t    *default_scheme;
  int                 istat;
  double              init_time;
  init_time = -sc_MPI_Wtime ();

  /* allocate problem */
  problem = T8_ALLOC (t8dg_linear_advection_problem_t, 1);

  problem->total_time = -sc_MPI_Wtime ();

  /*initialize stats */
  for (istat = 0; istat < ADVECT_NUM_STATS; istat++) {
    sc_stats_init (&problem->stats[istat], advect_stat_names[istat]);
  }

  default_scheme = t8_scheme_new_default_cxx ();

  t8dg_debugf ("create uniform forest\n");
  problem->forest = t8_forest_new_uniform (cmesh, default_scheme, uniform_level, 1, comm);

  problem->description = description;
  problem->dim = dim;
  problem->uniform_refinement_level = uniform_level;
  problem->maximum_refinement_level = max_level;
  problem->coarse_geometry = coarse_geometry;
  problem->description = description;
  problem->time_data = time_data;
  problem->vtk_count = 0;
  problem->comm = comm;

  problem->adapt_fn = adapt_fn;

  problem->dg_values = t8dg_values_new_LGL_hypercube (dim, number_LGL_points, problem->coarse_geometry, problem->forest);

  problem->dof_values = t8dg_dof_values_new (problem->forest, t8dg_values_get_global_values_array (problem->dg_values));
  problem->dof_values_adapt = NULL;

  t8dg_values_interpolate_scalar_function_3d_time (problem->dg_values, problem->description->initial_condition_fn,
                                                   t8dg_timestepping_data_get_current_time (problem->time_data), problem->dof_values);

  if (problem->maximum_refinement_level > problem->uniform_refinement_level) {
    int                 ilevel;

    for (ilevel = problem->uniform_refinement_level; ilevel < problem->maximum_refinement_level; ilevel++) {
      /* initial adapt */
      t8dg_advect_problem_adapt (problem, 0);
      /* repartition */
      t8dg_advect_problem_partition (problem, 0);
      /* Re initialize the dof_values */
      t8dg_values_interpolate_scalar_function_3d_time (problem->dg_values, problem->description->initial_condition_fn,
                                                       t8dg_timestepping_data_get_current_time (problem->time_data), problem->dof_values);
    }
  }
  t8dg_advect_problem_accumulate_stat (problem, ADVECT_INIT, init_time + sc_MPI_Wtime ());
  return problem;
}

void
t8dg_advect_problem_destroy (t8dg_linear_advection_problem_t ** pproblem)
{
  t8dg_linear_advection_problem_t *problem;
  T8DG_ASSERT (pproblem != NULL);
  problem = *pproblem;
  if (problem == NULL) {
    return;
  }

  t8dg_advect_problem_accumulate_stat (problem, ADVECT_TOTAL, problem->total_time + sc_MPI_Wtime ());

  if (problem->description->flux_data != NULL) {
    T8DG_FREE (problem->description->flux_data);
    T8DG_FREE (problem->description);
  }

  t8dg_timestepping_data_destroy (&(problem->time_data));

  t8dg_coarse_geometry_destroy (&(problem->coarse_geometry));

  t8dg_dof_values_destroy (&problem->dof_values);
  t8dg_values_destroy (&problem->dg_values);
  /* Unref the forest */
  t8_forest_unref (&problem->forest);   /*unrefs coarse mesh as well */
  /* Free the problem and set pointer to NULL */
  T8_FREE (problem);
  *pproblem = NULL;
}

static void
t8dg_advect_time_derivative (t8dg_dof_values_t * dof_values, t8dg_dof_values_t * dof_change, const double t, const void *application_data)
{
  T8DG_ASSERT (application_data != NULL);
  t8dg_linear_advection_problem_t *problem = (t8dg_linear_advection_problem_t *) application_data;
  T8DG_ASSERT (dof_values == problem->dof_values);
  T8DG_ASSERT (t == t8dg_timestepping_data_get_current_time (problem->time_data));

  t8dg_dof_values_t  *dof_flux;
  dof_flux = t8dg_dof_values_duplicate (dof_change);
  double              ghost_exchange_time, ghost_waittime;

  t8dg_values_apply_stiffness_matrix_linear_flux_fn3D (problem->dg_values, problem->description->velocity_field,
                                                       problem->description->flux_data, t, dof_values, dof_change);

  t8_debugf ("A u\n");
  t8dg_dof_values_debug_print (dof_change);

  /*Ghost exchange */
  ghost_exchange_time = -sc_MPI_Wtime ();
  t8dg_dof_values_ghost_exchange (dof_values);
  t8dg_values_ghost_exchange (problem->dg_values);
  ghost_exchange_time += sc_MPI_Wtime ();

  t8dg_advect_problem_accumulate_stat (problem, ADVECT_GHOST_EXCHANGE, ghost_exchange_time);
  ghost_waittime = t8_forest_profile_get_ghostexchange_waittime (problem->forest);
  t8dg_advect_problem_accumulate_stat (problem, ADVECT_GHOST_WAIT, ghost_waittime);
  t8dg_advect_problem_accumulate_stat (problem, ADVECT_AMR, ghost_exchange_time + ghost_waittime);

  t8dg_numerical_linear_flux3D_fn lax_friedrich_NF;
  lax_friedrich_NF = t8dg_linear_numerical_flux3D_lax_friedrich_fn;
  double              velocity_bound = 1;

  t8dg_values_apply_boundary_integrals (problem->dg_values, dof_values, dof_flux, problem->description->velocity_field,
                                        problem->description->flux_data, lax_friedrich_NF, &velocity_bound, t);

  t8dg_dof_values_axpy (-1, dof_flux, dof_change);      /*subtract function */
  t8dg_dof_values_destroy (&dof_flux);

  t8_debugf ("A u  - B u\n");
  t8dg_dof_values_debug_print (dof_change);

  /*apply massinverse */
  t8dg_values_apply_inverse_mass_matrix (problem->dg_values, dof_change, dof_change);

  t8_debugf ("du/dt = M^-1(A u - B u)\n");
  t8dg_dof_values_debug_print (dof_change);

  /*updates dof_values_change */
  /*ghost dof_values */
  /*dudt = Mg + Au */
  /*receive ghosts */
  /*dudt -= Bu */
  /*dudt = M^-1 dudt */
  /*resize dof values to only internal values */
}

static void
t8dg_advect_problem_set_time_step (t8dg_linear_advection_problem_t * problem)
{
  double              delta_t, min_delta_t, flow_velocity, time_left, diam, cfl;
  t8_locidx_t         num_trees, num_elems_in_tree, itree, ielement;
  t8_element_t       *element;
  double             *tree_vertices;

  /* maximum possible delta_t value */
  time_left = t8dg_timestepping_data_get_time_left (problem->time_data);
  min_delta_t = time_left;
  cfl = t8dg_timestepping_data_get_cfl (problem->time_data);

  num_trees = t8_forest_get_num_local_trees (problem->forest);
  for (itree = 0; itree < num_trees; itree++) {
    num_elems_in_tree = t8_forest_get_tree_num_elements (problem->forest, itree);

    tree_vertices = t8_forest_get_tree_vertices (problem->forest, itree);
    for (ielement = 0; ielement < num_elems_in_tree; ielement++) {
      element = t8_forest_get_element_in_tree (problem->forest, itree, ielement);
      /* Compute the minimum diameter */
      diam = t8_forest_element_diam (problem->forest, itree, element, tree_vertices);
      T8_ASSERT (diam > 0);

      flow_velocity = 1;        /*TODO: element_get_flow_velocity function */
      /* Compute minimum necessary time step */
      delta_t = time_left;
      if (flow_velocity > 0) {
        delta_t = cfl * diam / flow_velocity;
      }
      min_delta_t = SC_MIN (delta_t, min_delta_t);
    }
  }
  sc_MPI_Allreduce (&min_delta_t, &delta_t, 1, sc_MPI_DOUBLE, sc_MPI_MIN, problem->comm);
  t8dg_timestepping_data_set_time_step (problem->time_data, delta_t);
}

void
t8dg_advect_problem_advance_timestep (t8dg_linear_advection_problem_t * problem)
{
  double              solve_time = -sc_MPI_Wtime ();
  t8dg_advect_problem_set_time_step (problem);
  t8dg_advect_problem_accumulate_stat (problem, ADVECT_ELEM_AVG, t8_forest_get_global_num_elements (problem->forest));
  t8dg_timestepping_runge_kutta_step (t8dg_advect_time_derivative, problem->time_data, &(problem->dof_values), problem);
  t8dg_timestepping_data_increase_step_number (problem->time_data);
  t8dg_advect_problem_accumulate_stat (problem, ADVECT_SOLVE, solve_time + sc_MPI_Wtime ());
}

void
t8dg_advect_problem_adapt (t8dg_linear_advection_problem_t * problem, int measure_time)
{
  /* Nothing to do */
  if (problem->maximum_refinement_level == problem->uniform_refinement_level)
    return;

  t8_forest_t         forest_adapt;
  double              adapt_time, ghost_time, balance_time, replace_time;
  int                 did_balance = 0, balance_rounds = 0;
  t8_locidx_t         ghosts_sent;

  t8_debugf ("Into advect adapt\n");
  /* Adapt the forest, but keep the old one */
  t8_forest_ref (problem->forest);
  t8_forest_init (&forest_adapt);
  /* Enable profiling to measure the runtime */
  t8_forest_set_profiling (forest_adapt, 1);

  t8dg_adapt_data_t   adapt_data;
  adapt_data.dg_values = problem->dg_values;
  adapt_data.dof_values = problem->dof_values;
  adapt_data.maximum_refinement_level = problem->maximum_refinement_level;
  adapt_data.uniform_refinement_level = problem->uniform_refinement_level;

  /* Set the user data pointer of the new forest */
  t8_forest_set_user_data (forest_adapt, &adapt_data);
  /* Set the adapt function */
  t8_forest_set_adapt (forest_adapt, problem->forest, problem->adapt_fn, 0);
  if (problem->maximum_refinement_level - problem->uniform_refinement_level > 1) {
    /* We also want to balance the forest if there is a possibility of elements
     * with difference in refinement levels greater 1 */
    t8_forest_set_balance (forest_adapt, NULL, 1);
    did_balance = 1;
  }
  /* We also want ghost elements in the new forest */
  t8_forest_set_ghost (forest_adapt, 1, T8_GHOST_FACES);
  /* Commit the forest, adaptation and balance happens here */
  t8_forest_commit (forest_adapt);

  if (measure_time) {
    adapt_time = t8_forest_profile_get_adapt_time (forest_adapt);
    ghost_time = t8_forest_profile_get_ghost_time (forest_adapt, &ghosts_sent);
    if (did_balance) {
      balance_time = t8_forest_profile_get_balance_time (forest_adapt, &balance_rounds);
    }
    t8dg_advect_problem_accumulate_stat (problem, ADVECT_ADAPT, adapt_time);
    t8dg_advect_problem_accumulate_stat (problem, ADVECT_GHOST, ghost_time);
    t8dg_advect_problem_accumulate_stat (problem, ADVECT_GHOST_SENT, ghosts_sent);
    balance_time = 0;
    if (did_balance) {
      t8dg_advect_problem_accumulate_stat (problem, ADVECT_BALANCE, balance_time);
      t8dg_advect_problem_accumulate_stat (problem, ADVECT_BALANCE_ROUNDS, balance_rounds);
    }
  }

  t8dg_values_allocate_adapt (problem->dg_values, forest_adapt);

  problem->dof_values_adapt = t8dg_dof_values_new (forest_adapt, t8dg_values_get_global_values_array (problem->dg_values));

  adapt_data.dof_values_adapt = problem->dof_values_adapt;

  /* We now call iterate_replace in which we interpolate the new element data.
   * It is necessary that the old and new forest only differ by at most one level.
   * We guarantee this by calling adapt non-recursively and calling balance without
   * repartitioning. */
  replace_time = -sc_MPI_Wtime ();
  t8_forest_iterate_replace (forest_adapt, problem->forest, t8dg_adapt_replace);
  replace_time += sc_MPI_Wtime ();
  if (measure_time) {
    t8dg_advect_problem_accumulate_stat (problem, ADVECT_REPLACE, replace_time);
    t8dg_advect_problem_accumulate_stat (problem, ADVECT_AMR, ghost_time + adapt_time + balance_time + replace_time);
  }

  t8dg_values_cleanup_adapt (problem->dg_values);

  /* Free memory for the forest */
  t8_forest_unref (&problem->forest);
  /* Set the forest to the adapted one */
  problem->forest = forest_adapt;
  forest_adapt = NULL;

  t8dg_dof_values_destroy (&problem->dof_values);
  /* Set the phi values to the adapted phi values */
  problem->dof_values = problem->dof_values_adapt;
  problem->dof_values_adapt = NULL;
}

void
t8dg_advect_problem_partition (t8dg_linear_advection_problem_t * problem, int measure_time)
{
  t8_forest_t         forest_partition;
  t8dg_dof_values_t  *dof_values_partition;
  double              partition_time, ghost_time = 0, partition_data_time;
  int                 procs_sent, ghosts_sent;

  t8_forest_ref (problem->forest);
  t8_forest_init (&forest_partition);

  /* Enable profiling to measure the runtime */
  t8_forest_set_profiling (forest_partition, 1);

  t8_forest_set_partition (forest_partition, problem->forest, 0);
  t8_forest_set_ghost (forest_partition, 1, T8_GHOST_FACES);
  t8_forest_commit (forest_partition);

  if (measure_time) {
    partition_time = t8_forest_profile_get_partition_time (forest_partition, &procs_sent);
    ghost_time = t8_forest_profile_get_ghost_time (forest_partition, &ghosts_sent);
    t8dg_advect_problem_accumulate_stat (problem, ADVECT_PARTITION, partition_time);
    t8dg_advect_problem_accumulate_stat (problem, ADVECT_GHOST, ghost_time);
    t8dg_advect_problem_accumulate_stat (problem, ADVECT_GHOST_SENT, ghosts_sent);
    t8dg_advect_problem_accumulate_stat (problem, ADVECT_PARTITION_PROCS, procs_sent);
  }

  t8_debugf ("[ADVECT] partition");

  t8dg_values_partition (problem->dg_values, forest_partition);

  /* Partition dof_values */
  dof_values_partition = t8dg_dof_values_new (forest_partition, t8dg_values_get_global_values_array (problem->dg_values));

  partition_data_time = -sc_MPI_Wtime ();
  t8dg_dof_values_partition (problem->dof_values, dof_values_partition);
  partition_data_time += sc_MPI_Wtime ();

  if (measure_time) {
    t8dg_advect_problem_accumulate_stat (problem, ADVECT_PARTITION_DATA, partition_data_time);
    t8dg_advect_problem_accumulate_stat (problem, ADVECT_AMR, partition_time + ghost_time + partition_data_time);
  }

  t8_debugf (" [ADVECT] Done partition dof_values\n");

  /*destroy old dof values and use partition dof values */
  t8dg_dof_values_destroy (&problem->dof_values);
  problem->dof_values = dof_values_partition;

  t8_forest_unref (&problem->forest);
  problem->forest = forest_partition;
}

double
t8dg_advect_problem_l_infty_rel (t8dg_linear_advection_problem_t * problem)
{
  double              l_infty_rel;
  if (problem->description->analytical_sol_fn != NULL) {
    l_infty_rel = t8dg_values_norm_l_infty_rel (problem->dg_values, problem->dof_values, problem->description->analytical_sol_fn,
                                                t8dg_timestepping_data_get_current_time (problem->time_data), problem->comm);
    t8dg_advect_problem_accumulate_stat (problem, ADVECT_ERROR_INF, l_infty_rel);

  }
  return -1;
}

double
t8dg_advect_problem_l2_rel (t8dg_linear_advection_problem_t * problem)
{
  double              l2_rel;
  if (problem->description->analytical_sol_fn != NULL) {
    l2_rel = t8dg_values_norm_l2_rel (problem->dg_values, problem->dof_values, problem->description->analytical_sol_fn,
                                      t8dg_timestepping_data_get_current_time (problem->time_data), problem->comm);
    t8dg_advect_problem_accumulate_stat (problem, ADVECT_ERROR_2, l2_rel);

  }
  return -1;
}

void
t8dg_advect_problem_write_vtk (t8dg_linear_advection_problem_t * problem)
{
  char                fileprefix[BUFSIZ] = "t8dg_advect";
  double              io_time = -sc_MPI_Wtime ();

  t8dg_output_write_vtk (problem->dof_values, problem->forest, fileprefix, &problem->vtk_count);

  t8dg_advect_problem_accumulate_stat (problem, ADVECT_IO, io_time + sc_MPI_Wtime ());
}

void
t8dg_advect_problem_printdof (t8dg_linear_advection_problem_t * problem)
{
  t8dg_dof_values_debug_print (problem->dof_values);
}

int
t8dg_advect_problem_get_stepnumber (t8dg_linear_advection_problem_t * problem)
{
  return t8dg_timestepping_data_get_step_number (problem->time_data);
}

int
t8dg_advect_problem_endtime_reached (t8dg_linear_advection_problem_t * problem)
{
  return t8dg_timestepping_data_is_endtime_reached (problem->time_data);
}
