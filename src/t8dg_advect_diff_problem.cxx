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
#include "t8dg_advect_diff_problem.h"
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
#include "t8dg_common.h"
#include "t8dg_cmesh.h"

/* Enum for statistics. */
typedef enum
{
  ADVECT_DIFF_ADAPT = 0,        /* adapt runtime */
  ADVECT_DIFF_PARTITION,        /* partition runtime */
  ADVECT_DIFF_BALANCE,          /* balance runtime */
  ADVECT_DIFF_GHOST,            /* ghost runtime */
  ADVECT_DIFF_REPLACE,          /* forest_iterate_replace runtime */
  ADVECT_DIFF_IO,               /* vtk runtime */
  ADVECT_DIFF_INIT,             /* initialization runtime */
  ADVECT_DIFF_AMR,              /* AMR runtime (adapt+partition+ghost+balance) including data exchange (partition/ghost) */
  ADVECT_DIFF_SOLVE,            /* solver runtime */
  ADVECT_DIFF_TOTAL,            /* overall runtime */

  ADVECT_DIFF_ELEM_AVG,         /* average global number of elements (per time step) */
  ADVECT_DIFF_ERROR_INF,        /* l_infty error */
  ADVECT_DIFF_ERROR_2,          /* L_2 error */
  ADVECT_DIFF_VOL_LOSS,         /* The loss in volume in percent */

  ADVECT_DIFF_NUM_STATS,        /* The number of statistics that we measure */
  ADVECT_DIFF_NUM_TIME_STATS = ADVECT_DIFF_TOTAL + 1    /* The number of time statistics that we only want to have counted once */
} advect_diff_stats_t;

/* Names of statistics that we measure */
const char         *advect_diff_stat_names[ADVECT_DIFF_NUM_STATS] = {
  "adapt",
  "partition",
  "balance",
  "ghost",
  "replace",
  "vtk_print",
  "init",
  "AMR",
  "solve",
  "total",

  "number_elements",

  "l_infty_error",
  "l_2_error",
  "mass_loss_[%]"
};

/** The container for all information needed to solve the advection problem.
 */
struct t8dg_linear_advection_diffusion_problem
{
  t8_forest_t         forest;                                   /**< The t8_forest used to adapt the coarse mesh*/
  int                 dim;                                      /**< Dimension of the submanifold */

  t8dg_adapt_data    *adapt_data;

  t8dg_linear_advection_diffusion_problem_description_t *description;

  t8dg_timestepping_data_t *time_data;

  t8dg_vtk_data_t    *vtk_data;

  t8dg_values_t      *dg_values;

  /* The dof_values get ghosted */
  t8dg_dof_values_t  *dof_values;              /**< The Value of u at the nodal basis vertices */
  t8dg_dof_values_t  *dof_values_adapt;

  double              total_time;       /* used for measuring the total time the struct is in existence */
  sc_statinfo_t       stats[ADVECT_DIFF_NUM_STATS]; /**< Runtimes and other statistics. */
  sc_MPI_Comm         comm; /**< MPI Communicator */
};

void
t8dg_advect_diff_problem_compute_and_print_stats (t8dg_linear_advection_diffusion_problem_t * problem)
{
  sc_stats_compute (problem->comm, ADVECT_DIFF_NUM_STATS, problem->stats);
  sc_stats_print (t8dg_get_package_id (), SC_LP_ESSENTIAL, ADVECT_DIFF_NUM_STATS, problem->stats, 1, 1);
}

static void
t8dg_advect_diff_problem_accumulate_stat (t8dg_linear_advection_diffusion_problem_t * problem, advect_diff_stats_t stat, double value)
{
  T8DG_ASSERT (stat >= 0 && stat < ADVECT_DIFF_NUM_STATS);
  if (stat < ADVECT_DIFF_NUM_TIME_STATS) {
    sc_stats_set1 (&problem->stats[stat], problem->stats[stat].sum_values + value, advect_diff_stat_names[stat]);
  }
  else {
    sc_stats_accumulate (&problem->stats[stat], value);
  }
}

t8dg_linear_advection_diffusion_problem_description_t *
t8dg_advect_diff_problem_description_new (int initial_cond_arg, int velocity_field_arg, double flow_velocity, double diffusion_coefficient,
                                          int numerical_flux_arg, int source_sink_arg, int dim)
{
  t8dg_linear_advection_diffusion_problem_description_t *description;
  description = T8DG_ALLOC_ZERO (t8dg_linear_advection_diffusion_problem_description_t, 1);
  description->dim = dim;
  description->initial_condition_fn = t8dg_common_initial_cond_fn (initial_cond_arg);
  description->analytical_sol_fn = t8dg_common_analytic_solution_fn (initial_cond_arg, diffusion_coefficient);
  description->diffusion_coefficient = diffusion_coefficient;
  if (description->analytical_sol_fn == t8dg_scalar3d_sin_product) {
    t8dg_scalar3d_sin_product_data_t *ana_sol_data = T8DG_ALLOC_ZERO (t8dg_scalar3d_sin_product_data_t, 1);
    ana_sol_data->diffusion_coefficient = diffusion_coefficient;
    ana_sol_data->dim = dim;
    description->analytical_sol_data = ana_sol_data;

    t8dg_scalar3d_sin_product_data_t *init_data = T8DG_ALLOC_ZERO (t8dg_scalar3d_sin_product_data_t, 1);
    init_data->diffusion_coefficient = diffusion_coefficient;
    init_data->dim = dim;
    description->initial_condition_data = init_data;
  }
  t8dg_linear_flux3D_constant_flux_data_t *flux_data;
  description->numerical_flux_advection = t8dg_linear_numerical_flux3D_lax_friedrich_fn;
  description->numerical_flux_advection_data = T8DG_ALLOC (double, 1);

  switch (velocity_field_arg) {
  case 0:
    description->velocity_field = t8dg_linear_flux3D_constant_flux_fn;
    flux_data = T8DG_ALLOC_ZERO (t8dg_linear_flux3D_constant_flux_data_t, 1);
    flux_data->flow_direction[0] = 1;
    flux_data->flow_direction[1] = 0;
    flux_data->flow_direction[2] = 0;
    flux_data->flow_velocity = flow_velocity;
    *(double *) description->numerical_flux_advection_data = flux_data->flow_velocity;
    break;

  case 1:
    description->velocity_field = t8dg_rotating_flux_2D_fn;
    flux_data = NULL;
    *(double *) description->numerical_flux_advection_data = 1;

  case 2:
    description->velocity_field = t8dg_spiral_flux_3D_fn;
    flux_data = NULL;
    *(double *) description->numerical_flux_advection_data = 2 * M_PI;

  default:
    flux_data = NULL;
    break;
  }

  description->flux_data = flux_data;
  description->source_sink_fn = NULL;
  description->source_sink_data = NULL;
  /*boundary conditions */

  switch (numerical_flux_arg) {
  case 0:
    description->numerical_flux_diffusion_concentration = t8dg_numerical_flux1D_central;
    description->numerical_flux_diffusion_concentration_data = NULL;
    description->numerical_flux_diffusion_gradient = t8dg_numerical_flux1D_central;
    description->numerical_flux_diffusion_gradient_data = NULL;
    break;
  case 1:
    description->numerical_flux_diffusion_concentration = t8dg_numerical_flux1D_right;
    description->numerical_flux_diffusion_concentration_data = NULL;
    description->numerical_flux_diffusion_gradient = t8dg_numerical_flux1D_left;
    description->numerical_flux_diffusion_gradient_data = NULL;
    break;
  case 2:
    description->numerical_flux_diffusion_concentration = t8dg_numerical_flux1D_left;
    description->numerical_flux_diffusion_concentration_data = NULL;
    description->numerical_flux_diffusion_gradient = t8dg_numerical_flux1D_right;
    description->numerical_flux_diffusion_gradient_data = NULL;
    break;

  default:
    break;
  }

  switch (source_sink_arg) {
  case 0:
    description->source_sink_fn = NULL;
    description->source_sink_data = NULL;
    break;
  case 1:
    description->source_sink_fn = t8dg_cylinder_ring_source_fn;
    description->source_sink_data = NULL;
    break;
  default:
    break;
  }

  return description;
}

void
t8dg_advect_diff_problem_description_destroy (t8dg_linear_advection_diffusion_problem_description_t ** p_description)
{
  t8dg_linear_advection_diffusion_problem_description_t *description = *p_description;
  T8DG_FREE (description->flux_data);
  T8DG_FREE (description->analytical_sol_data);
  T8DG_FREE (description->initial_condition_data);
  T8DG_FREE (description->numerical_flux_advection_data);
  T8DG_FREE (description->numerical_flux_diffusion_concentration_data);
  T8DG_FREE (description->numerical_flux_diffusion_gradient_data);
  T8DG_FREE (description->source_sink_data);
  T8DG_FREE (description);
  *p_description = NULL;
}

t8dg_linear_advection_diffusion_problem_t *
t8dg_advect_diff_problem_init_arguments (int icmesh,
                                         int initial_level,
                                         int number_LGL_points,
                                         int initial_cond_arg,
                                         double flow_speed,
                                         double diffusion_coefficient,
                                         double start_time,
                                         double end_time,
                                         double cfl,
                                         double delta_t,
                                         int time_order,
                                         int min_level,
                                         int max_level,
                                         int adapt_arg,
                                         int adapt_freq, const char *prefix, int vtk_freq, int numerical_flux_arg, int source_sink_arg,
                                         sc_MPI_Comm comm)
{
  int                 dim;
  int                 geometry_arg;
  int                 velocity_field_arg;
  t8_scheme_cxx_t    *default_scheme;
  t8_cmesh_t          cmesh;
  t8_forest_t         forest;

  t8dg_coarse_geometry_t *coarse_geometry;
  t8dg_timestepping_data_t *time_data;
  t8dg_adapt_data_t  *adapt_data;
  t8dg_values_t      *dg_values;
  t8dg_vtk_data_t    *vtk_data;
  t8dg_linear_advection_diffusion_problem_description_t *description;
  double              init_time;
  init_time = -sc_MPI_Wtime ();

  default_scheme = t8_scheme_new_default_cxx ();
  cmesh = t8dg_cmesh_new_arg (icmesh, &dim, &velocity_field_arg, &geometry_arg, comm);
  forest = t8_forest_new_uniform (cmesh, default_scheme, initial_level, 1, comm);

  coarse_geometry = t8dg_coarse_geometry_new_arg (geometry_arg);

  description =
    t8dg_advect_diff_problem_description_new (initial_cond_arg, velocity_field_arg, flow_speed, diffusion_coefficient, numerical_flux_arg,
                                              source_sink_arg, dim);

  dg_values = t8dg_values_new_LGL_hypercube (dim, number_LGL_points, coarse_geometry, forest);

  adapt_data =
    t8dg_adapt_data_new (dg_values, initial_level, min_level, max_level, adapt_arg, adapt_freq, description->source_sink_fn,
                         description->source_sink_data);

  vtk_data = t8dg_output_vtk_data_new (prefix, vtk_freq);
  if (cfl > 0) {
    time_data = t8dg_timestepping_data_new_cfl (time_order, start_time, end_time, cfl);
  }
  else {
    time_data = t8dg_timestepping_data_new_constant_timestep (time_order, start_time, end_time, delta_t);
  }
  init_time += sc_MPI_Wtime ();

  return t8dg_advect_diff_problem_init (forest, description, dg_values, time_data, adapt_data, vtk_data, init_time, comm);
}

t8dg_linear_advection_diffusion_problem_t *
t8dg_advect_diff_problem_init (t8_forest_t forest, t8dg_linear_advection_diffusion_problem_description_t * description,
                               t8dg_values_t * dg_values, t8dg_timestepping_data_t * time_data,
                               t8dg_adapt_data_t * adapt_data, t8dg_vtk_data_t * vtk_data, double init_time, sc_MPI_Comm comm)
{
  t8dg_linear_advection_diffusion_problem_t *problem;
  int                 istat;
  init_time -= sc_MPI_Wtime ();

  /* allocate problem */
  problem = T8_ALLOC (t8dg_linear_advection_diffusion_problem_t, 1);

  problem->total_time = init_time;

  /*initialize stats */
  for (istat = 0; istat < ADVECT_DIFF_NUM_STATS; istat++) {
    sc_stats_init (&problem->stats[istat], advect_diff_stat_names[istat]);
  }

  problem->forest = forest;

  problem->vtk_data = vtk_data;
  problem->adapt_data = adapt_data;
  problem->dim = description->dim;
  problem->description = description;
  problem->time_data = time_data;
  problem->comm = comm;

  problem->dg_values = dg_values;

  problem->dof_values = t8dg_dof_values_new (problem->forest, t8dg_values_get_global_values_array (problem->dg_values));
  problem->dof_values_adapt = NULL;

  t8dg_values_interpolate_scalar_function_3d_time (problem->dg_values, problem->description->initial_condition_fn,
                                                   t8dg_timestepping_data_get_current_time (problem->time_data),
                                                   problem->description->initial_condition_data, problem->dof_values);

  int                 iadapt, adapt_steps;

  adapt_steps =
    SC_MAX (problem->adapt_data->maximum_refinement_level - problem->adapt_data->initial_refinement_level,
            problem->adapt_data->initial_refinement_level - problem->adapt_data->minimum_refinement_level);

  for (iadapt = 0; iadapt < adapt_steps; iadapt++) {
    /* initial adapt */
    t8dg_advect_diff_problem_adapt (problem, 0);
    /* repartition */
    t8dg_advect_diff_problem_partition (problem, 0);
    /* Re initialize the dof_values */
    t8dg_values_interpolate_scalar_function_3d_time (problem->dg_values, problem->description->initial_condition_fn,
                                                     t8dg_timestepping_data_get_current_time (problem->time_data),
                                                     problem->description->initial_condition_data, problem->dof_values);
  }

  t8dg_advect_diff_problem_accumulate_stat (problem, ADVECT_DIFF_INIT, init_time + sc_MPI_Wtime ());
  return problem;
}

void
t8dg_advect_diff_problem_destroy (t8dg_linear_advection_diffusion_problem_t ** pproblem)
{
  t8dg_linear_advection_diffusion_problem_t *problem;
  T8DG_ASSERT (pproblem != NULL);
  problem = *pproblem;
  if (problem == NULL) {
    return;
  }

  t8dg_advect_diff_problem_accumulate_stat (problem, ADVECT_DIFF_TOTAL, problem->total_time + sc_MPI_Wtime ());

  t8dg_advect_diff_problem_compute_and_print_stats (problem);

  t8dg_advect_diff_problem_description_destroy (&problem->description);

  t8dg_timestepping_data_destroy (&(problem->time_data));
  t8dg_adapt_data_destroy (&problem->adapt_data);
  t8dg_output_vtk_data_destroy (&problem->vtk_data);

  t8dg_dof_values_destroy (&problem->dof_values);
  t8dg_values_destroy (&problem->dg_values);
  /* Unref the forest */

  t8_forest_unref (&problem->forest);   /*unrefs coarse mesh as well */
  /* Free the problem and set pointer to NULL */
  T8_FREE (problem);
  *pproblem = NULL;
}

void
t8dg_advect_diff_solve (t8dg_linear_advection_diffusion_problem_t * problem)
{
  int                 step_number, apx_total_steps, modulus;

  t8dg_debugf ("Dof values beginning:\n");
  t8dg_advect_diff_problem_printdof (problem);

  t8dg_advect_diff_problem_set_time_step (problem);
  apx_total_steps = t8dg_advect_diff_problem_get_apx_total_steps (problem);
  modulus = SC_MAX (apx_total_steps / 10, 1);

  /*Timeloop with Rungekutta timestepping: */
  while (!t8dg_advect_diff_problem_endtime_reached (problem)) {
    step_number = t8dg_advect_diff_problem_get_stepnumber (problem);

    if (problem->vtk_data->vtk_freq && step_number % problem->vtk_data->vtk_freq == 0) {

      t8dg_advect_diff_problem_write_vtk (problem);
    }
    if ((step_number + 1) % modulus == 0) {
      t8dg_global_essentialf ("Step %i of apx. %i\n", step_number + 1, apx_total_steps);
    }

    t8dg_advect_diff_problem_advance_timestep (problem);

    if (problem->adapt_data->adapt_freq && (step_number + 1) % problem->adapt_data->adapt_freq == 0) {
      t8dg_advect_diff_problem_adapt (problem, 1);
      t8dg_advect_diff_problem_partition (problem, 1);
    }
  }

  t8dg_advect_diff_problem_write_vtk (problem);

  t8dg_advect_diff_problem_l2_rel (problem);
  t8dg_advect_diff_problem_l_infty_rel (problem);

  t8dg_debugf ("End Dof values:\n");
  t8dg_advect_diff_problem_printdof (problem);
}

static void
t8dg_advect_diff_time_derivative (t8dg_dof_values_t * dof_values, t8dg_dof_values_t * dof_change, const double t,
                                  const void *application_data)
{
  T8DG_ASSERT (application_data != NULL);
  t8dg_linear_advection_diffusion_problem_t *problem = (t8dg_linear_advection_diffusion_problem_t *) application_data;
  T8DG_ASSERT (dof_values == problem->dof_values);
  T8DG_ASSERT (t == t8dg_timestepping_data_get_current_time (problem->time_data));

  double              sqrt_diffusion_coefficient;
  t8dg_dof_values_t  *bu;
  t8dg_dof_values_t  *gradient_component;
  t8dg_dof_values_t  *gradient_component_stiffness;
  t8dg_dof_values_t  *gradient_component_boundary;
  t8dg_dof_values_t  *dof_summand;
  t8dg_dof_values_t  *dof_sum;

  sqrt_diffusion_coefficient = sqrt (problem->description->diffusion_coefficient);

  gradient_component = t8dg_dof_values_duplicate (dof_values);
  gradient_component_stiffness = t8dg_dof_values_duplicate (dof_values);
  gradient_component_boundary = t8dg_dof_values_duplicate (dof_values);
  dof_summand = t8dg_dof_values_duplicate (dof_values);
  dof_sum = t8dg_dof_values_duplicate (dof_values);

  bu = t8dg_dof_values_clone (dof_values);
  t8dg_dof_values_ax (bu, sqrt_diffusion_coefficient);

  t8dg_dof_values_set_zero (dof_sum);

  if (sqrt_diffusion_coefficient > 0) {
    int                 icomp;
    for (icomp = 0; icomp < problem->dim; icomp++) {
      t8dg_values_apply_component_stiffness_matrix_dof (problem->dg_values, icomp, bu, gradient_component_stiffness);

      t8dg_values_apply_component_boundary_integrals (problem->dg_values, bu, gradient_component_boundary, icomp,
                                                      problem->description->numerical_flux_diffusion_gradient,
                                                      problem->description->numerical_flux_diffusion_gradient_data, t);

      t8dg_dof_values_axpyz (-1, gradient_component_stiffness, gradient_component_boundary, gradient_component);

      t8dg_values_apply_inverse_mass_matrix (problem->dg_values, gradient_component, gradient_component);

      t8dg_dof_values_ax (gradient_component, sqrt_diffusion_coefficient);

      t8dg_values_apply_component_boundary_integrals (problem->dg_values, gradient_component, dof_summand, icomp,
                                                      problem->description->numerical_flux_diffusion_concentration,
                                                      problem->description->numerical_flux_diffusion_concentration_data, t);
      t8dg_dof_values_add (dof_sum, dof_summand);

      t8dg_values_apply_component_stiffness_matrix_dof (problem->dg_values, icomp, gradient_component, dof_summand);
      t8dg_dof_values_subtract (dof_sum, dof_summand);
    }

  }
  t8dg_values_apply_stiffness_matrix_linear_flux_fn3D (problem->dg_values, problem->description->velocity_field,
                                                       problem->description->flux_data, t, dof_values, dof_summand);

  t8dg_dof_values_add (dof_sum, dof_summand);

  t8dg_debugf ("stiffnessmatrix:\n");
  t8dg_dof_values_debug_print (dof_summand);

  t8dg_values_apply_boundary_integrals (problem->dg_values, dof_values, dof_summand, problem->description->velocity_field,
                                        problem->description->flux_data, problem->description->numerical_flux_advection,
                                        problem->description->numerical_flux_advection_data, t);

  t8dg_dof_values_subtract (dof_sum, dof_summand);      /*subtract function */
  t8dg_debugf ("boundary_matrix:\n");
  t8dg_dof_values_debug_print (dof_summand);

  if (problem->description->source_sink_fn != NULL) {
    t8dg_values_interpolate_scalar_function_3d_time (problem->dg_values, problem->description->source_sink_fn, t,
                                                     problem->description->source_sink_data, dof_summand);
    t8dg_values_apply_mass_matrix (problem->dg_values, dof_summand, dof_summand);
    t8dg_dof_values_add (dof_sum, dof_summand);
  }

  /*apply massinverse */
  t8dg_values_apply_inverse_mass_matrix (problem->dg_values, dof_sum, dof_change);

  t8_debugf ("du/dt = ");
  t8dg_dof_values_debug_print (dof_change);

  t8dg_dof_values_destroy (&gradient_component);
  t8dg_dof_values_destroy (&gradient_component_stiffness);
  t8dg_dof_values_destroy (&gradient_component_boundary);
  t8dg_dof_values_destroy (&bu);
  t8dg_dof_values_destroy (&dof_sum);
  t8dg_dof_values_destroy (&dof_summand);
}

int
t8dg_advect_diff_problem_get_apx_total_steps (t8dg_linear_advection_diffusion_problem_t * problem)
{
  double              timestep = t8dg_timestepping_data_get_time_step (problem->time_data);
  double              timelenght = t8dg_timestepping_data_get_time_left (problem->time_data);
  return (int) timelenght / timestep;
}

void
t8dg_advect_diff_problem_set_time_step (t8dg_linear_advection_diffusion_problem_t * problem)
{
  double              delta_t, min_delta_t, flow_velocity, time_left, diam, cfl;
  t8_locidx_t         num_trees, num_elems_in_tree, itree, ielement;
  t8_element_t       *element;

  /* maximum possible delta_t value */
  time_left = t8dg_timestepping_data_get_time_left (problem->time_data);
  min_delta_t = time_left;
  cfl = t8dg_timestepping_data_get_cfl (problem->time_data);

  if (cfl <= 0) {
    /*No CFL criterion given, constant timestep used except in last step */
    delta_t = SC_MAX (t8dg_timestepping_data_get_time_step (problem->time_data), min_delta_t);
    t8dg_timestepping_data_set_time_step (problem->time_data, delta_t);
    return;
  }

  num_trees = t8_forest_get_num_local_trees (problem->forest);
  for (itree = 0; itree < num_trees; itree++) {
    num_elems_in_tree = t8_forest_get_tree_num_elements (problem->forest, itree);

    for (ielement = 0; ielement < num_elems_in_tree; ielement++) {
      element = t8_forest_get_element_in_tree (problem->forest, itree, ielement);
      /* Compute the minimum diameter */
      diam = t8_forest_element_diam (problem->forest, itree, element);
      T8_ASSERT (diam > 0);
      if (problem->description->diffusion_coefficient > 0) {
        delta_t = cfl * diam * diam;
      }
      else if (problem->description->velocity_field == t8dg_linear_flux3D_constant_flux_fn) {
        flow_velocity = ((t8dg_linear_flux3D_constant_flux_data_t *) problem->description->flux_data)->flow_velocity;   /*TODO: element_get_flow_velocity function */
        delta_t = cfl * diam / flow_velocity;
      }
      else if (problem->description->velocity_field == t8dg_rotating_flux_2D_fn) {
        delta_t = cfl * diam;
      }
      else if (problem->description->velocity_field == t8dg_spiral_flux_3D_fn) {
        delta_t = cfl * diam / (2 * M_PI);
      }
      else {
        T8DG_ABORT ("Not implemented \n ");
      }
      min_delta_t = SC_MIN (delta_t, min_delta_t);
    }
  }
  sc_MPI_Allreduce (&min_delta_t, &delta_t, 1, sc_MPI_DOUBLE, sc_MPI_MIN, problem->comm);
  t8dg_timestepping_data_set_time_step (problem->time_data, delta_t);
}

void
t8dg_advect_diff_problem_advance_timestep (t8dg_linear_advection_diffusion_problem_t * problem)
{
  double              solve_time = -sc_MPI_Wtime ();
  t8dg_advect_diff_problem_set_time_step (problem);
  t8dg_advect_diff_problem_accumulate_stat (problem, ADVECT_DIFF_ELEM_AVG, t8_forest_get_global_num_elements (problem->forest));
  t8dg_timestepping_runge_kutta_step (t8dg_advect_diff_time_derivative, problem->time_data, &(problem->dof_values), problem);
  t8dg_timestepping_data_increase_step_number (problem->time_data);
  t8dg_advect_diff_problem_accumulate_stat (problem, ADVECT_DIFF_SOLVE, solve_time + sc_MPI_Wtime ());
}

void
t8dg_advect_diff_problem_adapt (t8dg_linear_advection_diffusion_problem_t * problem, int measure_time)
{
  /* Nothing to do */
  if (problem->adapt_data->maximum_refinement_level == problem->adapt_data->minimum_refinement_level)
    return;

  t8_forest_t         forest_adapt;
  double              adapt_time, ghost_time, balance_time = 0, replace_time;
  int                 did_balance = 0, balance_rounds = 0;
  t8_locidx_t         ghosts_sent;

  t8_debugf ("Into advect adapt\n");
  /* Adapt the forest, but keep the old one */
  t8_forest_ref (problem->forest);
  t8_forest_init (&forest_adapt);
  /* Enable profiling to measure the runtime */
  t8_forest_set_profiling (forest_adapt, 1);

  problem->adapt_data->dof_values = problem->dof_values;
  t8dg_adapt_data_interpolate_source_fn (problem->adapt_data);

  /* Set the user data pointer of the new forest */
  t8_forest_set_user_data (forest_adapt, problem->adapt_data);
  /* Set the adapt function */
  t8_forest_set_adapt (forest_adapt, problem->forest, problem->adapt_data->adapt_fn, 0);
  if (problem->adapt_data->maximum_refinement_level - problem->adapt_data->minimum_refinement_level > 1) {
    /* We also want to balance the forest if there is a possibility of elements
     * with difference in refinement levels greater 1 */
    t8_forest_set_balance (forest_adapt, NULL, 1);
    did_balance = 1;
  }
  /* We also want ghost elements in the new forest */
  t8_forest_set_ghost (forest_adapt, 1, T8_GHOST_FACES);
  /* Commit the forest, adaptation and balance happens here */
  t8_forest_commit (forest_adapt);

  t8dg_dof_values_destroy (&problem->adapt_data->source_sink_dof);

  if (measure_time) {
    adapt_time = t8_forest_profile_get_adapt_time (forest_adapt);
    ghost_time = t8_forest_profile_get_ghost_time (forest_adapt, &ghosts_sent);
    if (did_balance) {
      balance_time = t8_forest_profile_get_balance_time (forest_adapt, &balance_rounds);
    }
    t8dg_advect_diff_problem_accumulate_stat (problem, ADVECT_DIFF_ADAPT, adapt_time);
    t8dg_advect_diff_problem_accumulate_stat (problem, ADVECT_DIFF_GHOST, ghost_time);
    if (did_balance) {
      t8dg_advect_diff_problem_accumulate_stat (problem, ADVECT_DIFF_BALANCE, balance_time);
    }
  }

  t8dg_values_allocate_adapt (problem->dg_values, forest_adapt);

  problem->dof_values_adapt = t8dg_dof_values_new (forest_adapt, t8dg_values_get_global_values_array (problem->dg_values));

  problem->adapt_data->dof_values_adapt = problem->dof_values_adapt;

  /* We now call iterate_replace in which we interpolate the new element data.
   * It is necessary that the old and new forest only differ by at most one level.
   * We guarantee this by calling adapt non-recursively and calling balance without
   * repartitioning. */
  replace_time = -sc_MPI_Wtime ();
  t8_forest_iterate_replace (forest_adapt, problem->forest, t8dg_adapt_replace);
  replace_time += sc_MPI_Wtime ();
  if (measure_time) {
    t8dg_advect_diff_problem_accumulate_stat (problem, ADVECT_DIFF_REPLACE, replace_time);
    t8dg_advect_diff_problem_accumulate_stat (problem, ADVECT_DIFF_AMR, ghost_time + adapt_time + balance_time + replace_time);
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
t8dg_advect_diff_problem_partition (t8dg_linear_advection_diffusion_problem_t * problem, int measure_time)
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

  t8_debugf ("[ADVECT] partition");
  partition_data_time = -sc_MPI_Wtime ();
  t8dg_values_partition (problem->dg_values, forest_partition);

  /* Partition dof_values */
  dof_values_partition = t8dg_dof_values_new (forest_partition, t8dg_values_get_global_values_array (problem->dg_values));

  t8dg_dof_values_partition (problem->dof_values, dof_values_partition);
  partition_data_time += sc_MPI_Wtime ();

  if (measure_time) {
    partition_time = t8_forest_profile_get_partition_time (forest_partition, &procs_sent);
    ghost_time = t8_forest_profile_get_ghost_time (forest_partition, &ghosts_sent);
    t8dg_advect_diff_problem_accumulate_stat (problem, ADVECT_DIFF_PARTITION, partition_time + partition_data_time);
    t8dg_advect_diff_problem_accumulate_stat (problem, ADVECT_DIFF_GHOST, ghost_time);
    t8dg_advect_diff_problem_accumulate_stat (problem, ADVECT_DIFF_AMR, partition_time + ghost_time + partition_data_time);
  }

  t8_debugf (" [ADVECT] Done partition dof_values\n");

  /*destroy old dof values and use partition dof values */
  t8dg_dof_values_destroy (&problem->dof_values);
  problem->dof_values = dof_values_partition;

  t8_forest_unref (&problem->forest);
  problem->forest = forest_partition;
}

double
t8dg_advect_diff_problem_l_infty_rel (t8dg_linear_advection_diffusion_problem_t * problem)
{
  double              l_infty_rel;
  if (problem->description->analytical_sol_fn != NULL) {
    l_infty_rel = t8dg_values_norm_l_infty_rel (problem->dg_values, problem->dof_values, problem->description->analytical_sol_fn,
                                                t8dg_timestepping_data_get_current_time (problem->time_data),
                                                problem->description->analytical_sol_data, problem->comm);
    t8dg_advect_diff_problem_accumulate_stat (problem, ADVECT_DIFF_ERROR_INF, l_infty_rel);

  }
  return -1;
}

double
t8dg_advect_diff_problem_l2_rel (t8dg_linear_advection_diffusion_problem_t * problem)
{
  double              l2_rel;
  if (problem->description->analytical_sol_fn != NULL) {
    l2_rel = t8dg_values_norm_l2_rel (problem->dg_values, problem->dof_values, problem->description->analytical_sol_fn,
                                      t8dg_timestepping_data_get_current_time (problem->time_data),
                                      problem->description->analytical_sol_data, problem->comm);
    t8dg_advect_diff_problem_accumulate_stat (problem, ADVECT_DIFF_ERROR_2, l2_rel);

  }
  return -1;
}

void
t8dg_advect_diff_problem_write_vtk (t8dg_linear_advection_diffusion_problem_t * problem)
{
  double              io_time = -sc_MPI_Wtime ();

  t8dg_output_write_vtk (problem->dof_values, problem->vtk_data);

  t8dg_advect_diff_problem_accumulate_stat (problem, ADVECT_DIFF_IO, io_time + sc_MPI_Wtime ());
}

void
t8dg_advect_diff_problem_printdof (t8dg_linear_advection_diffusion_problem_t * problem)
{
  t8dg_dof_values_debug_print (problem->dof_values);
}

int
t8dg_advect_diff_problem_get_stepnumber (t8dg_linear_advection_diffusion_problem_t * problem)
{
  return t8dg_timestepping_data_get_step_number (problem->time_data);
}

int
t8dg_advect_diff_problem_endtime_reached (t8dg_linear_advection_diffusion_problem_t * problem)
{
  return t8dg_timestepping_data_is_endtime_reached (problem->time_data);
}
