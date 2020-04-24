/*
 * t8dg_advect.c
 *
 *  Created on: Mar 22, 2020
 *      Author: lukas
 */

#include <t8.h>
#include <t8_forest.h>
#include <t8_vec.h>
#include <t8_schemes/t8_default_cxx.hxx>
#include <t8_forest_vtk.h>
#include <t8_vtk.h>
#include <t8_forest/t8_forest_iterate.h>
#include <t8_forest/t8_forest_partition.h>

#include <sc_containers.h>
#include <sc_statistics.h>

#include "t8dg.h"
#include "t8dg_coarse_geometry.h"
#include "t8dg_advect_problem.h"
#include "t8dg_flux.h"
#include "t8dg_global_precomputed_values.h"
#include "t8dg_local_precomputed_values.h"
#include "t8dg_precomputed_values.h"
#include "t8dg_quadrature.h"
#include "t8dg_functionbasis.h"
#include "t8dg_sc_array.h"
#include "t8dg_timestepping.h"
#include "t8dg_mortar.h"

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

typedef struct t8dg_linear_advection_problem_description
{
  t8dg_scalar_function_3d_time_fn initial_condition_fn;           /**< Initial condition function */
  t8dg_linear_flux_t *flux;
  t8dg_linear_numerical_flux_fn numerical_flux_fn;             /**< Approximation to the Riemann problem */
  t8dg_scalar_function_3d_time_fn source_sink_fn;

} t8dg_linear_advection_problem_description_t;

/** The container for all information needed to solve the advection problem.
 */
struct t8dg_linear_advection_problem
{
  int                 dim;      /**< Dimension of the submanifold */

  int                 uniform_refinement_level; /**< uniform refinement level */
  int                 maximum_refinement_level;

  t8_forest_t         forest;                                   /**< The t8_forest used to adapt the coarse mesh*/

  t8dg_linear_advection_problem_description_t description;
  t8dg_timestepping_data_t *time_data;
  t8dg_global_precomputed_values_t *global_values;
  t8dg_coarse_geometry_t *coarse_geometry;   /**< coarse geometry, that gives the geometry and Jacobian of the geometry from the coarse
						 reference element to the coarse image element*/
  int                 vtk_count;
  sc_MPI_Comm         comm; /**< MPI Communicator */
  sc_statinfo_t       stats[ADVECT_NUM_STATS]; /**< Runtimes and other statistics. */

/* The following sc_arrays contain data, that needs to be available for each processorlocal element. To facilitate partitioning, they are
 * saved in a linear array
 */
  t8dg_local_precomputed_values_t *local_values;
  t8dg_local_precomputed_values_t *local_values_adapt;

  /* The dof_values get ghosted */
  sc_array_t         *dof_values;       /**< The Value of u at the nodal basis vertices */
  sc_array_t         *dof_values_adapt;

  /* those need to be recalculated for each time step, remain processor local */
  sc_array_t         *face_mortar[MAX_FACES];                   /**< contains pointer to face_mortars, so that fluxes need only be calculated once */

};

/* Getter functions for the sc_array_contained in advect_linear_problem.
 * TODO: setter functions
 *   */

static double      *
t8dg_advect_problem_get_element_dof_values (const t8dg_linear_advection_problem_t * problem, t8_locidx_t idata)
{
  return ((double *) t8_sc_array_index_locidx (problem->dof_values, idata));
}

/*  get functions for structs at element and faces: */

static t8dg_mortar_t *
t8dg_advect_element_get_face_mortar (const t8dg_linear_advection_problem_t * problem, t8_locidx_t idata, int iface)
{
  T8_ASSERT (iface >= 0 && iface < MAX_FACES);
  return *((t8dg_mortar_t **) t8_sc_array_index_locidx (problem->face_mortar[iface], idata));
}

static void
t8dg_advect_element_set_face_mortar (const t8dg_linear_advection_problem_t * problem, t8_locidx_t idata, int iface, t8dg_mortar_t * mortar)
{
  T8_ASSERT (iface >= 0 && iface < MAX_FACES);
  *((t8dg_mortar_t **) t8_sc_array_index_locidx (problem->face_mortar[iface], idata)) = mortar;
}

/* Given an allocated mortar, set it for both elements adjacent to the face. */
static void
t8dg_advect_element_set_face_mortar_both (const t8dg_linear_advection_problem_t * problem, t8dg_mortar_t * mortar)
{
  t8_locidx_t         idata;
  int                 iface;

  t8dg_mortar_get_idata_iface (mortar, &idata, &iface, 0);
  if (idata < t8_forest_get_num_element (problem->forest)) {
    t8dg_advect_element_set_face_mortar (problem, idata, iface, mortar);
  }
  t8dg_mortar_get_idata_iface (mortar, &idata, &iface, 1);
  if (idata < t8_forest_get_num_element (problem->forest)) {
    t8dg_advect_element_set_face_mortar (problem, idata, iface, mortar);
  }
}

/* Given an allocated mortar, set it for both elements adjacent to the face. */
static void
t8dg_advect_element_set_face_mortar_both_to_NULL (const t8dg_linear_advection_problem_t * problem, t8dg_mortar_t * mortar)
{
  t8_locidx_t         idata;
  int                 iface;

  t8dg_mortar_get_idata_iface (mortar, &idata, &iface, 0);
  if (idata < t8_forest_get_num_element (problem->forest)) {
    t8dg_advect_element_set_face_mortar (problem, idata, iface, NULL);
  }
  t8dg_mortar_get_idata_iface (mortar, &idata, &iface, 1);
  if (idata < t8_forest_get_num_element (problem->forest)) {
    t8dg_advect_element_set_face_mortar (problem, idata, iface, NULL);
  }
}

t8dg_timestepping_data_t *
t8dg_advect_get_time_data (t8dg_linear_advection_problem_t * problem)
{
  return problem->time_data;
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

void
t8dg_advect_problem_printdof (t8dg_linear_advection_problem_t * problem)
{
  t8dg_sc_array_block_double_debug_print (problem->dof_values);
}

void
t8dg_advect_problem_compute_and_print_stats (t8dg_linear_advection_problem_t * problem)
{
  sc_stats_compute (problem->comm, ADVECT_NUM_STATS, problem->stats);
  sc_stats_print (t8dg_get_package_id (), SC_LP_ESSENTIAL, ADVECT_NUM_STATS, problem->stats, 1, 1);
}

void
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

void
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

struct initial_fn_data
{
  const t8dg_linear_advection_problem_t *problem;
  t8_eclass_scheme_c *scheme;
  t8_element_t       *element;
  t8_locidx_t         itree;
  t8dg_scalar_function_3d_time_fn function;
};
static double
t8dg_advect_problem_evaluate_function_on_reference_element (const double reference_vertex[3], void *scalar_fn_data)
{
  double              coarse_vertex[DIM3];
  double              image_vertex[DIM3];
  double              start_time;
  struct initial_fn_data *data;
  const t8dg_linear_advection_problem_t *problem;

  data = (struct initial_fn_data *) scalar_fn_data;
  problem = data->problem;
  start_time = t8dg_timestepping_data_get_current_time (problem->time_data);
  /* transform into coarse reference element */
  t8dg_local_precomputed_values_fine_to_coarse_geometry (reference_vertex, coarse_vertex, data->scheme, data->element);

  /* tree vertices are application data for linear geometry */
  problem->coarse_geometry->geometry (coarse_vertex, image_vertex, problem->forest, data->itree);

  /* apply initial condition function at image vertex and start time */
  return data->function (image_vertex, start_time);

}

static void
t8dg_element_set_dofs_initial (t8dg_linear_advection_problem_t * problem, t8_locidx_t itree, t8_eclass_scheme_c * scheme,
                               t8_locidx_t ielement, t8_element_t * element)
{
  t8_locidx_t         idata;

  idata = t8dg_itree_ielement_to_idata (problem->forest, itree, ielement);

  sc_array_t         *element_dof_view = t8dg_sc_array_block_double_new_view (problem->dof_values, idata);

  struct initial_fn_data data = { problem, scheme, element, itree, problem->description.initial_condition_fn };

  t8dg_functionbasis_interpolate_scalar_fn (t8dg_global_precomputed_values_get_functionbasis (problem->global_values),
                                            t8dg_advect_problem_evaluate_function_on_reference_element, &data, element_dof_view);
  sc_array_destroy (element_dof_view);
}

double
t8dg_advect_problem_l_infty_rel (const t8dg_linear_advection_problem_t * problem, t8dg_scalar_function_3d_time_fn analytical_sol)
{
  t8_locidx_t         num_elements, num_trees;
  t8_locidx_t         ielement, itree, idata;
  sc_array_t         *elem_ana_sol;
  sc_array_t         *elem_dof_val;
  sc_array_t         *elem_error;
  double              error = 0, global_error;
  double              ana_norm = 0, global_ana_norm;

  struct initial_fn_data data = { problem, t8_forest_get_eclass_scheme (problem->forest, T8_ECLASS_LINE), NULL, 0,
    problem->description.initial_condition_fn
  };

  elem_ana_sol = sc_array_new_count (sizeof (double), t8dg_global_precomputed_values_get_num_dof (problem->global_values));
  elem_error = t8dg_sc_array_duplicate (elem_ana_sol);

  num_trees = t8_forest_get_num_local_trees (problem->forest);
  for (itree = 0, idata = 0; itree < num_trees; itree++) {
    data.itree = itree;
    num_elements = t8_forest_get_tree_num_elements (problem->forest, itree);
    for (ielement = 0; ielement < num_elements; ielement++, idata++) {
      data.element = t8_forest_get_element_in_tree (problem->forest, itree, ielement);
      elem_dof_val = t8dg_sc_array_block_double_new_view (problem->dof_values, idata);
      t8dg_functionbasis_interpolate_scalar_fn (t8dg_global_precomputed_values_get_functionbasis (problem->global_values),
                                                t8dg_advect_problem_evaluate_function_on_reference_element, &data, elem_ana_sol);
      t8dg_sc_array_block_double_axpyz (-1, elem_ana_sol, elem_dof_val, elem_error);
      error = SC_MAX (error, t8dg_precomputed_values_element_norm_infty (elem_error));
      ana_norm = SC_MAX (ana_norm, t8dg_precomputed_values_element_norm_infty (elem_ana_sol));
      sc_array_destroy (elem_dof_val);
    }
  }
  /* Compute the maximum of the error among all processes */
  sc_MPI_Allreduce (&ana_norm, &global_ana_norm, 1, sc_MPI_DOUBLE, sc_MPI_MAX, problem->comm);
  sc_MPI_Allreduce (&error, &global_error, 1, sc_MPI_DOUBLE, sc_MPI_MAX, problem->comm);

  sc_array_destroy (elem_ana_sol);
  sc_array_destroy (elem_error);

  /* Return the relative error, that is the l_infty error divided by
   * the l_infty norm of the analytical solution */
  return global_error / global_ana_norm;
}

double
t8dg_advect_problem_l2_rel (const t8dg_linear_advection_problem_t * problem, t8dg_scalar_function_3d_time_fn analytical_sol)
{
  t8_locidx_t         num_elements, num_trees;
  t8_locidx_t         ielement, itree, idata;
  sc_array_t         *elem_ana_sol;
  sc_array_t         *elem_dof_val;
  sc_array_t         *elem_error;
  double              error = 0, global_error;
  double              ana_norm = 0, global_ana_norm;

  struct initial_fn_data data = { problem, t8_forest_get_eclass_scheme (problem->forest, T8_ECLASS_LINE), NULL, 0,
    problem->description.initial_condition_fn
  };

  elem_ana_sol = sc_array_new_count (sizeof (double), t8dg_global_precomputed_values_get_num_dof (problem->global_values));
  elem_error = t8dg_sc_array_duplicate (elem_ana_sol);

  num_trees = t8_forest_get_num_local_trees (problem->forest);
  for (itree = 0, idata = 0; itree < num_trees; itree++) {
    data.itree = itree;
    num_elements = t8_forest_get_tree_num_elements (problem->forest, itree);
    for (ielement = 0; ielement < num_elements; ielement++, idata++) {
      data.element = t8_forest_get_element_in_tree (problem->forest, itree, ielement);
      elem_dof_val = t8dg_sc_array_block_double_new_view (problem->dof_values, idata);
      t8dg_functionbasis_interpolate_scalar_fn (t8dg_global_precomputed_values_get_functionbasis (problem->global_values),
                                                t8dg_advect_problem_evaluate_function_on_reference_element, &data, elem_ana_sol);
      t8dg_sc_array_block_double_axpyz (-1, elem_ana_sol, elem_dof_val, elem_error);
      error += t8dg_precomputed_values_element_norm_l2_squared (elem_error, problem->global_values, problem->local_values, idata),
        ana_norm += t8dg_precomputed_values_element_norm_l2_squared (elem_ana_sol, problem->global_values, problem->local_values, idata);
      sc_array_destroy (elem_dof_val);
    }
  }
  /* Compute the sum of the error among all processes */
  sc_MPI_Allreduce (&ana_norm, &global_ana_norm, 1, sc_MPI_DOUBLE, sc_MPI_SUM, problem->comm);
  sc_MPI_Allreduce (&error, &global_error, 1, sc_MPI_DOUBLE, sc_MPI_SUM, problem->comm);

  sc_array_destroy (elem_ana_sol);
  sc_array_destroy (elem_error);

  global_ana_norm = sqrt (global_ana_norm);
  global_error = sqrt (global_error);

  /* Return the relative error, that is the l_2 error divided by
   * the l_2 norm of the analytical solution */
  return global_error / global_ana_norm;
}

static void
t8dg_advect_problem_mortars_fill (t8dg_linear_advection_problem_t * problem)
{
  t8_locidx_t         itree, ielement, idata;
  t8_locidx_t         num_trees, num_elems_in_tree;
  int                 iface, num_faces;
  t8dg_mortar_t      *mortar;

  num_trees = t8_forest_get_num_local_trees (problem->forest);
  for (itree = 0, idata = 0; itree < num_trees; itree++) {

    num_elems_in_tree = t8_forest_get_tree_num_elements (problem->forest, itree);
    for (ielement = 0; ielement < num_elems_in_tree; ielement++, idata++) {

      num_faces = t8dg_global_precomputed_values_get_num_faces (problem->global_values);
      for (iface = 0; iface < num_faces; iface++) {
        mortar = t8dg_advect_element_get_face_mortar (problem, idata, iface);

        /*for all faces check wether mortar is already allocated/computed */
        if (mortar == NULL) {
          mortar = t8dg_mortar_new (problem->forest, itree, ielement, iface);
          t8dg_advect_element_set_face_mortar_both (problem, mortar);
        }
        if (!(t8dg_mortar_is_valid (mortar))) {
          t8dg_mortar_fill (mortar, problem->dof_values, problem->time_data,
                            problem->global_values, problem->local_values,
                            problem->description.flux, problem->description.numerical_flux_fn);
        }
      }
    }
  }
}

static void
t8dg_advect_problem_mortars_new (t8dg_linear_advection_problem_t * problem)
{
  int                 iface;
  t8_locidx_t         num_elements, num_trees, num_elems_in_tree;
  t8_locidx_t         itree, ielement, idata;
  num_elements = t8_forest_get_num_element (problem->forest);

  for (iface = 0; iface < t8dg_global_precomputed_values_get_num_faces (problem->global_values); iface++) {
    problem->face_mortar[iface] = sc_array_new_count (sizeof (t8dg_mortar_t *), num_elements);
    num_trees = t8_forest_get_num_local_trees (problem->forest);

    for (itree = 0, idata = 0; itree < num_trees; itree++) {
      num_elems_in_tree = t8_forest_get_tree_num_elements (problem->forest, itree);
      for (ielement = 0; ielement < num_elems_in_tree; ielement++, idata++) {
        t8dg_advect_element_set_face_mortar (problem, idata, iface, NULL);
      }
    }
  }
}

static void
t8dg_advect_problem_mortars_invalidate (t8dg_linear_advection_problem_t * problem)
{
  t8_locidx_t         itree, ielement, idata;
  t8_locidx_t         num_trees, num_elems_in_tree;
  int                 iface, num_faces;
  num_trees = t8_forest_get_num_local_trees (problem->forest);

  t8dg_mortar_t      *mortar;

  for (itree = 0, idata = 0; itree < num_trees; itree++) {

    num_elems_in_tree = t8_forest_get_tree_num_elements (problem->forest, itree);

    for (ielement = 0; ielement < num_elems_in_tree; ielement++, idata++) {
      num_faces = t8dg_global_precomputed_values_get_num_faces (problem->global_values);
      for (iface = 0; iface < num_faces; iface++) {
        mortar = t8dg_advect_element_get_face_mortar (problem, idata, iface);
        if (mortar != NULL) {
          t8dg_mortar_invalidate (mortar);
        }
      }
    }
  }
}

static void
t8dg_advect_problem_mortars_destroy (t8dg_linear_advection_problem_t * problem)
{

  t8_locidx_t         itree, ielement, idata, iface;
  t8_locidx_t         num_trees, num_elems_in_tree, num_faces;

  t8dg_mortar_t      *mortar;

  num_trees = t8_forest_get_num_local_trees (problem->forest);

  for (itree = 0, idata = 0; itree < num_trees; itree++) {

    num_elems_in_tree = t8_forest_get_tree_num_elements (problem->forest, itree);

    for (ielement = 0; ielement < num_elems_in_tree; ielement++, idata++) {
      num_faces = t8dg_global_precomputed_values_get_num_faces (problem->global_values);
      for (iface = 0; iface < num_faces; iface++) {
        mortar = t8dg_advect_element_get_face_mortar (problem, idata, iface);
        if (mortar != NULL) {
          t8dg_advect_element_set_face_mortar_both_to_NULL (problem, mortar);
          t8dg_mortar_destroy (&mortar);
        }
      }

    }
  }
  for (iface = 0; iface < t8dg_global_precomputed_values_get_num_faces (problem->global_values); iface++) {
    sc_array_destroy (problem->face_mortar[iface]);
  }
}

static t8dg_linear_advection_problem_t *
t8dg_advect_problem_init (t8_cmesh_t cmesh,
                          t8dg_coarse_geometry_t * coarse_geometry,
                          int dim,
                          t8dg_scalar_function_3d_time_fn u_initial,
                          double flow_speed,
                          int uniform_level, int max_level,
                          int number_LGL_points, double start_time, double end_time, double cfl, int time_order, sc_MPI_Comm comm)
{
  t8dg_linear_advection_problem_t *problem;
  t8_scheme_cxx_t    *default_scheme;
  int                 num_elements;
  int                 istat;
  /* allocate problem */
  problem = T8_ALLOC (t8dg_linear_advection_problem_t, 1);

  default_scheme = t8_scheme_new_default_cxx ();
  t8_debugf ("create uniform forest\n");
  problem->forest = t8_forest_new_uniform (cmesh, default_scheme, uniform_level, 1, comm);

  problem->dim = dim;
  problem->uniform_refinement_level = uniform_level;
  problem->maximum_refinement_level = max_level;

  problem->coarse_geometry = coarse_geometry;

  double              tangential_vector[3] = { 1, 0, 0 };       /*TODO: make dependent on cmesh! */
  problem->description.initial_condition_fn = u_initial;
  problem->description.flux = t8dg_linear_flux_new_1D_linear_geometry (tangential_vector, flow_speed);
  problem->description.numerical_flux_fn = t8dg_linear_numerical_flux_upwind_1D;

  problem->time_data = t8dg_timestepping_data_new (time_order, start_time, end_time, cfl);
//  t8dg_timestepping_data_set_time_step (problem->time_data, cfl * pow (2, -uniform_level));

  problem->vtk_count = 0;
  problem->comm = comm;

  for (istat = 0; istat < ADVECT_NUM_STATS; istat++) {
    sc_stats_init (&problem->stats[istat], advect_stat_names[istat]);
  }

  t8_debugf ("precompute global values\n");
  /* these allocate memory: */
  problem->global_values = t8dg_global_precomputed_values_new_1D_LGL (number_LGL_points);

  t8dg_quadrature_t  *quadrature = t8dg_global_precomputed_values_get_quadrature (problem->global_values);

  t8_debugf ("precompute local values\n");
  num_elements = t8_forest_get_num_element (problem->forest);

  problem->local_values = t8dg_local_precomputed_values_new (quadrature, num_elements);
  problem->local_values_adapt = NULL;

  /* the dof_values need to be ghosted. */
  problem->dof_values =
    sc_array_new_count (sizeof (double) * t8dg_global_precomputed_values_get_num_dof (problem->global_values),
                        num_elements + t8_forest_get_num_ghosts (problem->forest));

  problem->dof_values_adapt = NULL;

  t8dg_advect_problem_mortars_new (problem);
  /*for each element and face a pointer to a mortar */

  t8_debugf ("start elementinit\n");

  t8dg_advect_problem_init_elements (problem);

  if (problem->maximum_refinement_level > problem->uniform_refinement_level) {
    int                 ilevel;

    for (ilevel = problem->uniform_refinement_level; ilevel < problem->maximum_refinement_level; ilevel++) {
      /* initial adapt */
      t8dg_advect_problem_adapt (problem, 0);
      /* repartition */
      t8dg_advect_problem_partition (problem, 0);
      /* Re initialize the elements */
      t8dg_advect_problem_init_elements (problem);
    }
  }

  return problem;
}

/*TODO: which init function creates what, outsource problem description*/
t8dg_linear_advection_problem_t *
t8dg_advect_problem_init_linear_geometry_1D (int icmesh,
                                             t8dg_scalar_function_3d_time_fn u_initial,
                                             double flow_speed,
                                             int uniform_level, int max_level,
                                             int number_LGL_points,
                                             double start_time, double end_time, double cfl, int time_order, sc_MPI_Comm comm)
{

  t8_cmesh_t          cmesh;
  switch (icmesh) {
  case 0:
    cmesh = t8_cmesh_new_hypercube (T8_ECLASS_LINE, sc_MPI_COMM_WORLD, 0, 0, 1);
    break;
  case 1:
    cmesh = t8_cmesh_new_periodic_line_more_trees (sc_MPI_COMM_WORLD);
    break;
  default:
    T8DG_ABORT ("Not yet implemented");
  }

  t8dg_coarse_geometry_t *coarse_geometry = t8dg_coarse_geometry_new_1D_linear ();
  return t8dg_advect_problem_init (cmesh, coarse_geometry, 1, u_initial, flow_speed, uniform_level, max_level,
                                   number_LGL_points, start_time, end_time, cfl, time_order, comm);
}

void
t8dg_advect_problem_init_elements (t8dg_linear_advection_problem_t * problem)
{
  t8_locidx_t         itree, ielement, idata;
  t8_locidx_t         num_trees, num_elems_in_tree;
  t8_element_t       *element;

  t8_eclass_scheme_c *scheme;

  t8dg_quadrature_t  *quadrature = t8dg_global_precomputed_values_get_quadrature (problem->global_values);

  t8_debugf ("Start element init \n");
  num_trees = t8_forest_get_num_local_trees (problem->forest);

  for (itree = 0, idata = 0; itree < num_trees; itree++) {
    scheme = t8_forest_get_eclass_scheme (problem->forest, t8_forest_get_tree_class (problem->forest, itree));

    num_elems_in_tree = t8_forest_get_tree_num_elements (problem->forest, itree);
    for (ielement = 0; ielement < num_elems_in_tree; ielement++, idata++) {
      element = t8_forest_get_element_in_tree (problem->forest, itree, ielement);

      t8dg_local_precomputed_values_set_element (problem->local_values, problem->forest, itree, scheme, ielement, quadrature);

      t8dg_element_set_dofs_initial (problem, itree, scheme, ielement, element);
    }
  }
  t8_debugf ("End element init \n");
}

void
t8dg_advect_problem_destroy (t8dg_linear_advection_problem_t ** pproblem)
{
  t8dg_linear_advection_problem_t *problem;

  T8_ASSERT (pproblem != NULL);
  problem = *pproblem;
  if (problem == NULL) {
    return;
  }

  t8dg_linear_flux_destroy (&(problem->description.flux));
  t8dg_timestepping_data_destroy (&(problem->time_data));

  t8dg_local_precomputed_values_destroy (&(problem->local_values));
  t8dg_advect_problem_mortars_destroy (problem);
  sc_array_destroy (problem->dof_values);
  t8dg_global_precomputed_values_destroy (&problem->global_values);
  t8dg_coarse_geometry_destroy (&(problem->coarse_geometry));

  /* Unref the forest */
  t8_forest_unref (&problem->forest);   /*unrefs coarse mesh as well */
  /* Free the problem and set pointer to NULL */
  T8_FREE (problem);
  *pproblem = NULL;
}

static void
t8dg_advect_problem_apply_stiffness_matrix (t8dg_linear_advection_problem_t * problem, sc_array_t * src_dof, sc_array_t * dest_dof)
{
  t8_locidx_t         itree, ielement, idata;
  t8_locidx_t         num_trees, num_elems_in_tree, num_dof;
  t8dg_quad_idx_t     num_quad_vertices;

  sc_array_t         *element_quad_values;
  sc_array_t         *element_dof_values;
  sc_array_t         *element_dof_derivative_values;
  sc_array_t         *element_res_dof_values;

  num_quad_vertices = t8dg_global_precomputed_values_get_num_elem_quad (problem->global_values);
  num_dof = t8dg_global_precomputed_values_get_num_dof (problem->global_values);

  num_trees = t8_forest_get_num_local_trees (problem->forest);
  for (itree = 0, idata = 0; itree < num_trees; itree++) {

    num_elems_in_tree = t8_forest_get_tree_num_elements (problem->forest, itree);

    for (ielement = 0; ielement < num_elems_in_tree; ielement++, idata++) {
      element_dof_values = t8dg_sc_array_block_double_new_view (src_dof, idata);
      element_quad_values = sc_array_new_count (sizeof (double), num_quad_vertices);
      element_dof_derivative_values = sc_array_new_count (sizeof (double), num_dof);
      element_res_dof_values = t8dg_sc_array_block_double_new_view (dest_dof, idata);

      t8dg_global_precomputed_values_transform_element_dof_to_element_quad (problem->global_values, element_dof_values,
                                                                            element_quad_values);
      t8_debugf ("element_quad_values\n");
      t8dg_sc_array_block_double_debug_print (element_quad_values);

      t8dg_local_precomputed_values_element_multiply_trafo_quad_weight (problem->local_values, element_quad_values, idata);
      t8_debugf ("element_quad_values*qtw\n");
      t8dg_sc_array_block_double_debug_print (element_quad_values);
      t8dg_flux_element_multiply_flux_value (problem->description.flux, element_quad_values,
                                             t8dg_timestepping_data_get_current_time (problem->time_data), problem->local_values,
                                             problem->forest, itree, ielement,
                                             t8dg_global_precomputed_values_get_quadrature (problem->global_values),
                                             problem->coarse_geometry);
      t8_debugf ("element_quad_values*qtw*flux_value\n");
      t8dg_sc_array_block_double_debug_print (element_quad_values);
      t8dg_global_precomputed_values_transform_element_quad_to_element_dof (problem->global_values, element_quad_values,
                                                                            element_dof_derivative_values);
      t8dg_global_precomputed_values_element_apply_derivative_matrix_transpose (problem->global_values, element_dof_derivative_values,
                                                                                element_res_dof_values);

      sc_array_destroy (element_dof_values);
      sc_array_destroy (element_quad_values);
      sc_array_destroy (element_dof_derivative_values);
      sc_array_destroy (element_res_dof_values);
    }
  }
}

static void
t8dg_advect_problem_apply_inverse_mass_matrix (t8dg_linear_advection_problem_t * problem, sc_array_t * dof_array,
                                               sc_array_t * result_dof_array)
{
  /*
   * src_dof = dest_dof is allowed
   **/
  t8_locidx_t         itree, ielement, idata;
  t8_locidx_t         num_trees, num_elems_in_tree;
  sc_array_t         *element_dof_array;
  sc_array_t         *element_result_dof_array;

  num_trees = t8_forest_get_num_local_trees (problem->forest);

  for (itree = 0, idata = 0; itree < num_trees; itree++) {

    num_elems_in_tree = t8_forest_get_tree_num_elements (problem->forest, itree);

    for (ielement = 0; ielement < num_elems_in_tree; ielement++, idata++) {
      element_dof_array = t8dg_sc_array_block_double_new_view (dof_array, idata);
      element_result_dof_array = t8dg_sc_array_block_double_new_view (result_dof_array, idata);
      t8dg_precomputed_values_apply_element_inverse_mass_matrix (problem->global_values, problem->local_values, idata,
                                                                 element_dof_array, element_result_dof_array);
      sc_array_destroy (element_dof_array);
      sc_array_destroy (element_result_dof_array);
    }
  }
}

static void
t8dg_advect_problem_apply_boundary_integrals (t8dg_linear_advection_problem_t * problem, sc_array_t * dest)
{
  t8_locidx_t         itree, ielement, idata;
  t8_locidx_t         num_trees, num_elems_in_tree;
  t8dg_quad_idx_t     num_elem_quad;
  int                 iface, num_faces;
  t8dg_mortar_t      *mortar;

  sc_array_t         *element_quad_flux;
  sc_array_t         *element_dest;

  double              alpha;

  num_faces = t8dg_global_precomputed_values_get_num_faces (problem->global_values);
  num_elem_quad = t8dg_global_precomputed_values_get_num_elem_quad (problem->global_values);

  t8dg_sc_array_block_double_set_zero (dest);

  num_trees = t8_forest_get_num_local_trees (problem->forest);
  for (itree = 0, idata = 0; itree < num_trees; itree++) {

    num_elems_in_tree = t8_forest_get_tree_num_elements (problem->forest, itree);

    for (ielement = 0; ielement < num_elems_in_tree; ielement++, idata++) {
      element_dest = t8dg_sc_array_block_double_new_view (dest, idata);

      element_quad_flux = sc_array_new_count (sizeof (double), num_elem_quad);

      for (iface = 0; iface < num_faces; iface++) {
        mortar = t8dg_advect_element_get_face_mortar (problem, idata, iface);
        t8dg_global_precomputed_values_transform_face_quad_to_element_dof
          (problem->global_values, iface, t8dg_mortar_get_flux (mortar), element_quad_flux);

        t8_debugf ("test face_quad_element,id:%i,if:%i\n", idata, iface);
        t8dg_sc_array_block_double_debug_print (element_quad_flux);
        /*decide wether to subtract or add the calculated fluxes! TODO:check */
        alpha = -t8dg_mortar_get_side (mortar, idata);

        t8dg_sc_array_block_double_axpy (alpha, element_quad_flux, element_dest);
      }
      t8_debugf ("test element_dest\n");
      t8dg_sc_array_block_double_debug_print (element_dest);

      sc_array_destroy (element_quad_flux);
      sc_array_destroy (element_dest);

    }
  }

}

static void
t8dg_advect_time_derivative (const sc_array_t * dof_values, sc_array_t * dof_change, const double t, const void *application_data)
{
  T8DG_ASSERT (application_data != NULL);
  t8dg_linear_advection_problem_t *problem = (t8dg_linear_advection_problem_t *) application_data;
  T8DG_ASSERT (dof_values == problem->dof_values);
  T8DG_ASSERT (t == t8dg_timestepping_data_get_current_time (problem->time_data));

  sc_array_t         *dof_flux;
  dof_flux = t8dg_sc_array_duplicate (dof_change);
  double              ghost_exchange_time, ghost_waittime;

  t8dg_advect_problem_apply_stiffness_matrix (problem, problem->dof_values, dof_change);

  t8_debugf ("A u\n");
  t8dg_sc_array_block_double_debug_print (dof_change);

  /*Ghost exchange */
  ghost_exchange_time = -sc_MPI_Wtime ();
  t8_forest_ghost_exchange_data (problem->forest, problem->dof_values);
  ghost_exchange_time += sc_MPI_Wtime ();
  t8dg_advect_problem_accumulate_stat (problem, ADVECT_GHOST_EXCHANGE, ghost_exchange_time);
  ghost_waittime = t8_forest_profile_get_ghostexchange_waittime (problem->forest);
  t8dg_advect_problem_accumulate_stat (problem, ADVECT_GHOST_WAIT, ghost_waittime);
  t8dg_advect_problem_accumulate_stat (problem, ADVECT_AMR, ghost_exchange_time + ghost_waittime);

  t8_debugf ("fill mortars\n");
  t8dg_advect_problem_mortars_fill (problem);
  t8_debugf ("mortars filled\n");
  t8dg_advect_problem_apply_boundary_integrals (problem, dof_flux);
  t8dg_advect_problem_mortars_invalidate (problem);

  t8dg_sc_array_block_double_axpy (-1, dof_flux, dof_change);
  sc_array_destroy (dof_flux);

  t8_debugf ("A u  - B u\n");
  t8dg_sc_array_block_double_debug_print (dof_change);

  /*apply massinverse */
  t8dg_advect_problem_apply_inverse_mass_matrix (problem, dof_change, dof_change);

  t8_debugf ("du/dt = M^-1(A u - B u)\n");
  t8dg_sc_array_block_double_debug_print (dof_change);

  /*updates dof_values_change */
  /*ghost dof_values */
  /*dudt = Mg + Au */
  /*receive ghosts */
  /*dudt -= Bu */
  /*dudt = M^-1 dudt */
  /*resize dof values to only internal values */
}

void
t8dg_advect_problem_advance_timestep (t8dg_linear_advection_problem_t * problem)
{
  t8dg_advect_problem_accumulate_stat (problem, ADVECT_ELEM_AVG, t8_forest_get_global_num_elements (problem->forest));
  t8dg_timestepping_runge_kutta_step (t8dg_advect_time_derivative, t8dg_advect_get_time_data (problem), &(problem->dof_values), problem);
  t8dg_timestepping_data_increase_step_number (problem->time_data);
}

void
t8dg_advect_write_vtk (t8dg_linear_advection_problem_t * problem)
{
  double             *dof_array;
  t8_locidx_t         num_local_elements, idata;
  t8_vtk_data_field_t vtk_data;
  char                fileprefix[BUFSIZ];
  double             *dof_values;
  double              average;
  int                 idof, number_of_dof;

  number_of_dof = t8dg_global_precomputed_values_get_num_dof (problem->global_values);

  num_local_elements = t8_forest_get_num_element (problem->forest);
  dof_array = T8_ALLOC_ZERO (double, num_local_elements);

  for (idata = 0; idata < num_local_elements; idata++) {
    dof_values = t8dg_advect_problem_get_element_dof_values (problem, idata);
    average = 0;
    for (idof = 0; idof < number_of_dof; idof++) {
      average += dof_values[idof];
    }
    average /= number_of_dof;
    dof_array[idata] = average;
  }

  /* Write meta data for vtk */
  snprintf (vtk_data.description, BUFSIZ, "Num. Solution");
  vtk_data.type = T8_VTK_SCALAR;
  vtk_data.data = dof_array;
  /* Write filename */
  snprintf (fileprefix, BUFSIZ, "t8dg_advection_%03i", problem->vtk_count);
  /* Write vtk files */
  if (t8_forest_vtk_write_file (problem->forest, fileprefix, 1, 1, 1, 1, 0, 1, &vtk_data)) {
    t8_debugf ("[Advect] Wrote pvtu to files %s\n", fileprefix);
  }
  else {
    t8_errorf ("[Advect] Error writing to files %s\n", fileprefix);
  }
  /* clean-up */
  T8_FREE (dof_array);
  problem->vtk_count++;
}

#if 0
static int
t8dg_advect_test_adapt (t8_forest_t forest,
                        t8_forest_t forest_from,
                        t8_locidx_t which_tree,
                        t8_locidx_t lelement_id, t8_eclass_scheme_c * ts, int num_elements, t8_element_t * elements[])
{
  t8dg_linear_advection_problem_t *problem;
  t8_locidx_t         first_idata;
  double             *dof_values;
  int                 level;

  first_idata = t8dg_itree_ielement_to_idata (forest_from, which_tree, lelement_id);
  problem = (t8dg_linear_advection_problem_t *) t8_forest_get_user_data (forest);

  level = ts->t8_element_level (elements[0]);
  if (level == problem->maximum_refinement_level && num_elements == 1) {
    /* It is not possible to refine this level */
    return 0;
  }

  dof_values = t8dg_advect_problem_get_element_dof_values (problem, first_idata);
  if (dof_values[0] < 0.2) {
    return -(num_elements > 1 && level > problem->uniform_refinement_level);    //could discard second check
  }
  else if (dof_values[0] > 0.8) {
    return level < problem->maximum_refinement_level;
  }
  return 0;
}
#endif

static int
t8dg_advect_gradient_adapt (t8_forest_t forest,
                            t8_forest_t forest_from,
                            t8_locidx_t which_tree,
                            t8_locidx_t lelement_id, t8_eclass_scheme_c * ts, int num_elements, t8_element_t * elements[])
{
  t8dg_linear_advection_problem_t *problem;
  t8_locidx_t         first_idata;
  double             *dof_values;
  int                 level;
  double              diam;
  double             *tree_vertices;

  double              gradient_threshold_refine = 0.6;
  double              gradient_threshold_coarsen = 0.2;

  int                 num_dof;

  first_idata = t8dg_itree_ielement_to_idata (forest_from, which_tree, lelement_id);
  problem = (t8dg_linear_advection_problem_t *) t8_forest_get_user_data (forest);

  T8DG_CHECK_ABORT (problem->dim == 1, "Not yet implemented");

  num_dof = t8dg_global_precomputed_values_get_num_dof (problem->global_values);

  level = ts->t8_element_level (elements[0]);
  if (level == problem->maximum_refinement_level && num_elements == 1) {
    /* It is not possible to refine this level */
    return 0;
  }
  tree_vertices = t8_forest_get_tree_vertices (forest_from, which_tree);

  if (num_elements == 1) {
    dof_values = t8dg_advect_problem_get_element_dof_values (problem, first_idata);
    diam = t8_forest_element_diam (forest_from, which_tree, elements[0], tree_vertices);

    double              gradient = fabs (dof_values[0] - dof_values[num_dof - 1]) / diam;
    return gradient > gradient_threshold_refine;
  }
  else {
    dof_values = t8dg_advect_problem_get_element_dof_values (problem, first_idata);
    diam = t8_forest_element_diam (forest_from, which_tree, elements[0], tree_vertices);

    double              gradient_left = fabs (dof_values[0] - dof_values[num_dof - 1]) / diam;

    dof_values = t8dg_advect_problem_get_element_dof_values (problem, first_idata + 1);
    diam = t8_forest_element_diam (forest_from, which_tree, elements[1], tree_vertices);

    double              gradient_right = fabs (dof_values[0] - dof_values[num_dof - 1]) / diam;
    if (gradient_left > gradient_threshold_refine && level < problem->maximum_refinement_level)
      return 1;
    return -(gradient_left < gradient_threshold_coarsen && gradient_right < gradient_threshold_coarsen
             && level > problem->uniform_refinement_level);
  }
  return 0;
}

static void
t8dg_advect_test_replace (t8_forest_t forest_old,
                          t8_forest_t forest_new,
                          t8_locidx_t which_tree,
                          t8_eclass_scheme_c * ts,
                          int num_elems_old, t8_locidx_t first_ielem_old, int num_elems_new, t8_locidx_t first_ielem_new)
{
  t8dg_linear_advection_problem_t *problem;
  t8_locidx_t         first_idata_old, first_idata_new;

  int                 ichild, num_children;

  sc_array_t         *element_dof_parent;
  sc_array_t         *element_dof_child[MAX_SUBELEMENTS];

  t8dg_quadrature_t  *quadrature;

  problem = (t8dg_linear_advection_problem_t *) t8_forest_get_user_data (forest_new);
  T8DG_ASSERT (forest_old == problem->forest);
  T8DG_CHECK_ABORT (t8dg_global_precomputed_values_get_dim (problem->global_values) == 1, "Not yet implemented");

  quadrature = t8dg_global_precomputed_values_get_quadrature (problem->global_values);

  first_idata_old = t8dg_itree_ielement_to_idata (forest_old, which_tree, first_ielem_old);
  first_idata_new = t8dg_itree_ielement_to_idata (forest_new, which_tree, first_ielem_new);

  num_children = 2;             /*TODO: scheme dependent */

  if (num_elems_old == num_elems_new && num_elems_old == 1) {
    t8dg_local_precomputed_values_copy_element_values (problem->local_values, first_idata_old,
                                                       problem->local_values_adapt, first_idata_new);
    t8dg_sc_array_copy_only_at_indices (problem->dof_values, first_idata_old, problem->dof_values_adapt, first_idata_new);
  }
  else if (num_elems_old == 1) {
    element_dof_parent = t8dg_sc_array_block_double_new_view (problem->dof_values, first_idata_old);
    for (ichild = 0; ichild < num_children; ichild++) {
      element_dof_child[ichild] = t8dg_sc_array_block_double_new_view (problem->dof_values_adapt, first_idata_new + ichild);
      t8dg_global_precomputed_values_transform_element_dof_to_child_dof (problem->global_values, element_dof_parent,
                                                                         element_dof_child[ichild], ichild);
      t8dg_local_precomputed_values_set_element (problem->local_values_adapt, forest_new, which_tree, ts, first_ielem_new + ichild,
                                                 quadrature);
      sc_array_destroy (element_dof_child[ichild]);
    }
    sc_array_destroy (element_dof_parent);
  }
  else {
    /* Needed before! we calculate the dof_values */
    t8dg_local_precomputed_values_set_element (problem->local_values_adapt, forest_new, which_tree, ts, first_ielem_new, quadrature);
    element_dof_parent = t8dg_sc_array_block_double_new_view (problem->dof_values_adapt, first_idata_new);
    t8dg_sc_array_block_double_set_zero (element_dof_parent);
    for (ichild = 0; ichild < num_children; ichild++) {
      element_dof_child[ichild] = t8dg_sc_array_block_double_new_view (problem->dof_values, first_idata_old + ichild);

    }
    t8dg_precomputed_values_transform_child_dof_to_parent_dof (problem->global_values, element_dof_child, element_dof_parent, num_children,
                                                               problem->local_values, problem->local_values_adapt,
                                                               first_idata_old, first_idata_new);

    for (ichild = 0; ichild < num_children; ichild++) {
      sc_array_destroy (element_dof_child[ichild]);
    }
    sc_array_destroy (element_dof_parent);

  }
}

void
t8dg_advect_problem_adapt (t8dg_linear_advection_problem_t * problem, int measure_time)
{
  /* Nothing to do */
  if (problem->maximum_refinement_level == problem->uniform_refinement_level)
    return;

  t8_locidx_t         num_elems_p_ghosts, num_elems;
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

  /* Set the user data pointer of the new forest */
  t8_forest_set_user_data (forest_adapt, problem);
  /* Set the adapt function */
  t8_forest_set_adapt (forest_adapt, problem->forest, t8dg_advect_gradient_adapt, 0);
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
    if (did_balance) {
      t8dg_advect_problem_accumulate_stat (problem, ADVECT_BALANCE, balance_time);
      t8dg_advect_problem_accumulate_stat (problem, ADVECT_BALANCE_ROUNDS, balance_rounds);
    }
  }

  /* Allocate new memory for the element_data of the advected forest */
  num_elems = t8_forest_get_num_element (forest_adapt);
  num_elems_p_ghosts = num_elems + t8_forest_get_num_ghosts (forest_adapt);

  problem->local_values_adapt =
    t8dg_local_precomputed_values_new (t8dg_global_precomputed_values_get_quadrature (problem->global_values), num_elems);
  problem->dof_values_adapt =
    sc_array_new_count (t8dg_global_precomputed_values_get_num_dof (problem->global_values) * sizeof (double), num_elems_p_ghosts);

  /* We now call iterate_replace in which we interpolate the new element data.
   * It is necessary that the old and new forest only differ by at most one level.
   * We guarantee this by calling adapt non-recursively and calling balance without
   * repartitioning. */
  replace_time = -sc_MPI_Wtime ();
  t8_forest_iterate_replace (forest_adapt, problem->forest, t8dg_advect_test_replace);
  replace_time += sc_MPI_Wtime ();
  if (measure_time) {
    t8dg_advect_problem_accumulate_stat (problem, ADVECT_REPLACE, replace_time);
    t8dg_advect_problem_accumulate_stat (problem, ADVECT_AMR, ghost_time + adapt_time + balance_time + replace_time);
  }

  /* clean the old element data */
  t8dg_advect_problem_mortars_destroy (problem);
  t8dg_local_precomputed_values_destroy (&problem->local_values);
  sc_array_destroy (problem->dof_values);
  /* Free memory for the forest */
  t8_forest_unref (&problem->forest);
  /* Set the forest to the adapted one */
  problem->forest = forest_adapt;
  forest_adapt = NULL;
  /* Set the elem data to the adapted elem data */
  problem->local_values = problem->local_values_adapt;
  problem->local_values_adapt = NULL;
  /* Set the phi values to the adapted phi values */
  problem->dof_values = problem->dof_values_adapt;
  problem->dof_values_adapt = NULL;
  /*Create new mortar arrays */
  t8dg_advect_problem_mortars_new (problem);
}

void
t8dg_advect_problem_partition (t8dg_linear_advection_problem_t * problem, int measure_time)
{
  t8_forest_t         forest_partition;
  t8dg_local_precomputed_values_t *local_values_partition;
  sc_array_t         *dof_values_partition;
  sc_array_t         *dof_values_local_view;
  sc_array_t         *dof_values_partition_local_view;
  t8_locidx_t         num_local_elems_new, num_local_elems_old, num_ghosts_new;

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

  num_local_elems_old = t8_forest_get_num_element (problem->forest);
  num_local_elems_new = t8_forest_get_num_element (forest_partition);

  num_ghosts_new = t8_forest_get_num_ghosts (forest_partition);

  t8_debugf ("[ADVECT] partition with: num_old:%i, num_new:%i, ghost_new:%i\n", num_local_elems_old, num_local_elems_new, num_ghosts_new);

  /* Partition local precomputed values */
  local_values_partition = t8dg_local_precomputed_values_new (t8dg_global_precomputed_values_get_quadrature (problem->global_values),
                                                              num_local_elems_new);
  t8dg_local_precomputed_values_partition (problem->forest, forest_partition, problem->local_values, local_values_partition);

  t8dg_local_precomputed_values_destroy (&problem->local_values);
  problem->local_values = local_values_partition;

  t8_debugf ("[ADVECT] Done partition local_data\n");

  /* Partition dof_values */
  dof_values_partition = sc_array_new_count (t8dg_global_precomputed_values_get_num_dof (problem->global_values) * sizeof (double),
                                             num_local_elems_new + num_ghosts_new);

  dof_values_local_view = sc_array_new_view (problem->dof_values, 0, num_local_elems_old);
  dof_values_partition_local_view = sc_array_new_view (dof_values_partition, 0, num_local_elems_new);

  partition_data_time = -sc_MPI_Wtime ();
  t8_forest_partition_data (problem->forest, forest_partition, dof_values_local_view, dof_values_partition_local_view);
  partition_data_time += sc_MPI_Wtime ();

  if (measure_time) {
    t8dg_advect_problem_accumulate_stat (problem, ADVECT_PARTITION_DATA, partition_data_time);
    t8dg_advect_problem_accumulate_stat (problem, ADVECT_AMR, partition_time + ghost_time + partition_data_time);
  }

  t8_debugf (" [ADVECT] Done partition dof_values\n");

  /*destroy views */
  sc_array_destroy (dof_values_local_view);
  sc_array_destroy (dof_values_partition_local_view);

  /*destroy old dof values and use partition dof values */
  sc_array_destroy (problem->dof_values);
  problem->dof_values = dof_values_partition;

  t8_debugf (" [ADVECT] begin mortars destroy\n");
  t8dg_advect_problem_mortars_destroy (problem);
  t8_forest_unref (&problem->forest);
  problem->forest = forest_partition;
  t8dg_advect_problem_mortars_new (problem);
  t8_debugf (" [ADVECT] Done partition\n");
}
