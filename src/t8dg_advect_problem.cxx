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
#include "t8dg_geometry.h"

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

  struct t8dg_linear_advection_problem_description
  {
    t8dg_scalar_function_3d_time_fn initial_condition_fn;           /**< Initial condition function */
    t8dg_flux_t        *flux;
    t8dg_scalar_function_3d_time_fn source_sink_fn;
    t8dg_scalar_function_3d_time_fn analytical_sol_fn;           /**< Analytical solution function */
  } description;

  t8dg_timestepping_data_t *time_data;
  t8dg_global_precomputed_values_t *global_values;
  t8dg_coarse_geometry_t *coarse_geometry;   /**< coarse geometry, that gives the geometry and Jacobian of the geometry from the coarse
						 reference element to the coarse image element*/
  int                 vtk_count;
  sc_MPI_Comm         comm; /**< MPI Communicator */
  sc_statinfo_t       stats[ADVECT_NUM_STATS]; /**< Runtimes and other statistics. */

  int                 max_num_element_values, max_num_faces, max_num_face_values;
/* The following sc_arrays contain data, that needs to be available for each processorlocal element. To facilitate partitioning, they are
 * saved in a linear array
 */
  t8dg_local_precomputed_values_t *local_values;
  t8dg_local_precomputed_values_t *local_values_adapt;

  /* The dof_values get ghosted */
  sc_array_t         *dof_values;       /**< The Value of u at the nodal basis vertices */
  sc_array_t         *dof_values_adapt;

  /* those need to be recalculated for each time step, remain processor local */
  t8dg_mortar_array_t *face_mortars;                           /**< contains pointer to face_mortars, so that fluxes need only be calculated once */

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

double
t8dg_advect_problem_l_infty_rel (const t8dg_linear_advection_problem_t * problem)
{
  t8_locidx_t         num_elements, num_trees;
  t8_locidx_t         ielement, itree, idata;
  sc_array_t         *elem_ana_sol;
  sc_array_t         *elem_dof_val;
  sc_array_t         *elem_error;
  double              error = 0, global_error;
  double              ana_norm = 0, global_ana_norm;

  t8dg_geometry_transformation_data_t geometry_data = { problem->coarse_geometry, problem->forest, 0, 0 };

  t8dg_precomputed_values_fn_evaluation_data_t evaluation_data = { &geometry_data, problem->description.analytical_sol_fn,
    t8dg_timestepping_data_get_current_time (problem->time_data)
  };

  elem_ana_sol = sc_array_new_count (sizeof (double), t8dg_global_precomputed_values_get_num_dof (problem->global_values));
  elem_error = t8dg_sc_array_duplicate (elem_ana_sol);

  num_trees = t8_forest_get_num_local_trees (problem->forest);
  for (itree = 0, idata = 0; itree < num_trees; itree++) {
    geometry_data.itree = itree;
    num_elements = t8_forest_get_tree_num_elements (problem->forest, itree);
    for (ielement = 0; ielement < num_elements; ielement++, idata++) {
      geometry_data.ielement = ielement;
      elem_dof_val = t8dg_sc_array_block_double_new_view (problem->dof_values, idata);
      t8dg_functionbasis_interpolate_scalar_fn (t8dg_global_precomputed_values_get_functionbasis (problem->global_values),
                                                t8dg_precomputed_values_transform_reference_vertex_and_evaluate, &evaluation_data,
                                                elem_ana_sol);
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
t8dg_advect_problem_l2_rel (const t8dg_linear_advection_problem_t * problem)
{
  t8_locidx_t         num_elements, num_trees;
  t8_locidx_t         ielement, itree, idata;
  sc_array_t         *elem_ana_sol;
  sc_array_t         *elem_dof_val;
  sc_array_t         *elem_error;
  double              error = 0, global_error;
  double              ana_norm = 0, global_ana_norm;

  t8dg_geometry_transformation_data_t geometry_data = { problem->coarse_geometry, problem->forest, 0, 0 };

  t8dg_precomputed_values_fn_evaluation_data_t evaluation_data = { &geometry_data, problem->description.analytical_sol_fn,
    t8dg_timestepping_data_get_current_time (problem->time_data)
  };

  elem_ana_sol = sc_array_new_count (sizeof (double), t8dg_global_precomputed_values_get_num_dof (problem->global_values));
  elem_error = t8dg_sc_array_duplicate (elem_ana_sol);

  num_trees = t8_forest_get_num_local_trees (problem->forest);
  for (itree = 0, idata = 0; itree < num_trees; itree++) {
    geometry_data.itree = itree;
    num_elements = t8_forest_get_tree_num_elements (problem->forest, itree);
    for (ielement = 0; ielement < num_elements; ielement++, idata++) {
      geometry_data.ielement = ielement;
      elem_dof_val = t8dg_sc_array_block_double_new_view (problem->dof_values, idata);
      t8dg_functionbasis_interpolate_scalar_fn (t8dg_global_precomputed_values_get_functionbasis (problem->global_values),
                                                t8dg_precomputed_values_transform_reference_vertex_and_evaluate, &evaluation_data,
                                                elem_ana_sol);
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

t8dg_linear_advection_problem_t *
t8dg_advect_problem_init (t8_cmesh_t cmesh,
                          t8dg_coarse_geometry_t * coarse_geometry,
                          int dim,
                          t8dg_scalar_function_3d_time_fn u_initial,
                          t8dg_flux_t * flux,
                          int uniform_level, int max_level, int number_LGL_points, t8dg_timestepping_data_t * time_data, sc_MPI_Comm comm)
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

  problem->description.initial_condition_fn = u_initial;
  problem->description.analytical_sol_fn = u_initial;   /*Assumes that the solution is a whole number of revolutions around the periodic domain */
  problem->description.flux = flux;

  problem->time_data = time_data;

  problem->vtk_count = 0;
  problem->comm = comm;

  for (istat = 0; istat < ADVECT_NUM_STATS; istat++) {
    sc_stats_init (&problem->stats[istat], advect_stat_names[istat]);
  }

  t8_debugf ("precompute global values\n");
  /* these allocate memory: */
  problem->global_values = t8dg_global_precomputed_values_new_hypercube_LGL (dim, number_LGL_points);

  t8_debugf ("precompute local values\n");
  num_elements = t8_forest_get_num_element (problem->forest);

  problem->max_num_element_values = t8dg_global_precomputed_values_get_num_dof (problem->global_values);
  problem->max_num_faces = t8dg_global_precomputed_values_get_num_faces (problem->global_values);
  problem->max_num_face_values = t8dg_global_precomputed_values_get_max_num_facevalues (problem->global_values);
  problem->local_values =
    t8dg_local_precomputed_values_new (num_elements, problem->dim, problem->max_num_element_values, problem->max_num_faces,
                                       problem->max_num_face_values);
  problem->local_values_adapt = NULL;

  /* the dof_values need to be ghosted. */
  problem->dof_values =
    sc_array_new_count (sizeof (double) * t8dg_global_precomputed_values_get_num_dof (problem->global_values),
                        num_elements + t8_forest_get_num_ghosts (problem->forest));

  problem->dof_values_adapt = NULL;

  problem->face_mortars =
    t8dg_mortar_array_new_empty (problem->forest, t8dg_global_precomputed_values_get_num_faces (problem->global_values));

  t8_debugf ("start element init\n");

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

void
t8dg_advect_problem_init_elements (t8dg_linear_advection_problem_t * problem)
{
  t8_locidx_t         itree = 0, ielement, idata;
  t8_locidx_t         num_trees, num_elems_in_tree;
  sc_array_t         *element_dof_view;

  t8dg_geometry_transformation_data_t geometry_data = { problem->coarse_geometry, problem->forest, 0, 0 };

  t8dg_precomputed_values_fn_evaluation_data_t evaluation_data = { &geometry_data, problem->description.initial_condition_fn,
    t8dg_timestepping_data_get_current_time (problem->time_data)
  };

  t8_debugf ("Start element init \n");
  num_trees = t8_forest_get_num_local_trees (problem->forest);

  for (itree = 0, idata = 0; itree < num_trees; itree++) {
    geometry_data.itree = itree;

    num_elems_in_tree = t8_forest_get_tree_num_elements (problem->forest, itree);
    for (ielement = 0; ielement < num_elems_in_tree; ielement++, idata++) {
      geometry_data.ielement = ielement;
      t8dg_local_precomputed_values_set_element (problem->local_values, &geometry_data, problem->global_values);

      element_dof_view = t8dg_sc_array_block_double_new_view (problem->dof_values, idata);

      t8dg_functionbasis_interpolate_scalar_fn (t8dg_global_precomputed_values_get_functionbasis (problem->global_values),
                                                t8dg_precomputed_values_transform_reference_vertex_and_evaluate, &evaluation_data,
                                                element_dof_view);
      sc_array_destroy (element_dof_view);
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

  t8dg_flux_destroy (&(problem->description.flux));
  t8dg_timestepping_data_destroy (&(problem->time_data));

  t8dg_local_precomputed_values_destroy (&(problem->local_values));
  t8dg_mortar_array_destroy (&problem->face_mortars);
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
  int                 num_quad_vertices;
  int                 idim;

  sc_array_t         *element_quad_values;
  sc_array_t         *element_flux_quad_values;
  sc_array_t         *element_dof_values;
  sc_array_t         *element_res_dof_values;
  sc_array_t         *element_res_summand_dof_values;

  t8dg_geometry_transformation_data_t geometry_data = { problem->coarse_geometry, problem->forest, 0, 0 };

  num_quad_vertices = t8dg_global_precomputed_values_get_num_elem_quad (problem->global_values);
  num_dof = t8dg_global_precomputed_values_get_num_dof (problem->global_values);

  num_trees = t8_forest_get_num_local_trees (problem->forest);
  for (itree = 0, idata = 0; itree < num_trees; itree++) {
    geometry_data.itree = itree;

    num_elems_in_tree = t8_forest_get_tree_num_elements (problem->forest, itree);
    for (ielement = 0; ielement < num_elems_in_tree; ielement++, idata++) {
      geometry_data.ielement = ielement;
      element_dof_values = t8dg_sc_array_block_double_new_view (src_dof, idata);
      element_quad_values = sc_array_new_count (sizeof (double), num_quad_vertices);
      element_flux_quad_values = sc_array_new_count (sizeof (double), num_quad_vertices);
      element_res_summand_dof_values = sc_array_new_count (sizeof (double), num_dof);
      element_res_dof_values = t8dg_sc_array_block_double_new_view (dest_dof, idata);

      t8dg_local_precomputed_values_element_multiply_trafo_quad_weight (problem->local_values, idata, element_dof_values,
                                                                        element_quad_values);
//      t8_debugf ("element_quad_values*qtw\n");
//      t8dg_sc_array_block_double_debug_print (element_quad_values);

      t8dg_sc_array_block_double_set_zero (element_res_dof_values);
      for (idim = 0; idim < t8dg_global_precomputed_values_get_dim (problem->global_values); idim++) {

        t8dg_local_precomputed_values_element_multiply_flux_value (problem->local_values, problem->description.flux, &geometry_data,
                                                                   t8dg_global_precomputed_values_get_quadrature (problem->global_values),
                                                                   t8dg_timestepping_data_get_current_time (problem->time_data),
                                                                   idim, element_quad_values, element_flux_quad_values);
//        t8_debugf ("element_dof_derivative_values\n");
//        t8dg_sc_array_block_double_debug_print (element_dof_derivative_values);
        t8dg_global_precomputed_values_element_apply_derivative_matrix_transpose (problem->global_values, idim,
                                                                                  element_flux_quad_values, element_res_summand_dof_values);
//        t8_debugf ("element_res_summand_dof_values\n");
//        t8dg_sc_array_block_double_debug_print (element_res_summand_dof_values);
        t8dg_sc_array_block_double_axpy (1, element_res_summand_dof_values, element_res_dof_values);
      }

      sc_array_destroy (element_dof_values);
      sc_array_destroy (element_quad_values);
      sc_array_destroy (element_flux_quad_values);
      sc_array_destroy (element_res_dof_values);
      sc_array_destroy (element_res_summand_dof_values);
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
t8dg_advect_problem_apply_boundary_integrals (t8dg_linear_advection_problem_t * problem, sc_array_t * result_dof)
{
  t8_locidx_t         itree, ielement, idata;
  t8_locidx_t         num_trees, num_elems_in_tree;
  sc_array_t         *element_result_dof;

  num_trees = t8_forest_get_num_local_trees (problem->forest);
  for (itree = 0, idata = 0; itree < num_trees; itree++) {

    num_elems_in_tree = t8_forest_get_tree_num_elements (problem->forest, itree);

    for (ielement = 0; ielement < num_elems_in_tree; ielement++, idata++) {
      element_result_dof = t8dg_sc_array_block_double_new_view (result_dof, idata);
      t8dg_precomputed_values_apply_element_boundary_integral (problem->global_values, problem->local_values, problem->face_mortars, idata,
                                                               element_result_dof);
      sc_array_destroy (element_result_dof);
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
  t8dg_geometry_transformation_data_t geometry_data = { problem->coarse_geometry, problem->forest, 0, 0 };
  t8dg_mortar_fill_data_t mortar_fill_data = { problem->global_values, problem->local_values, &geometry_data, problem->description.flux,
    problem->dof_values, t
  };

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
  t8dg_mortar_array_fill (problem->face_mortars, &mortar_fill_data);
  t8_debugf ("mortars filled\n");
  t8dg_advect_problem_apply_boundary_integrals (problem, dof_flux);
  t8dg_mortar_array_invalidate_all (problem->face_mortars);

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

  problem = (t8dg_linear_advection_problem_t *) t8_forest_get_user_data (forest_new);
  T8DG_ASSERT (forest_old == problem->forest);
  T8DG_CHECK_ABORT (t8dg_global_precomputed_values_get_dim (problem->global_values) == 1, "Not yet implemented");

  t8dg_geometry_transformation_data_t geometry_data = { problem->coarse_geometry, forest_new, which_tree, 0 };

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
      geometry_data.ielement = first_ielem_new + ichild;
      t8dg_local_precomputed_values_set_element (problem->local_values_adapt, &geometry_data, problem->global_values);
      sc_array_destroy (element_dof_child[ichild]);
    }
    sc_array_destroy (element_dof_parent);
  }
  else {
    /* Needed before! we calculate the dof_values */
    geometry_data.ielement = first_ielem_new;
    t8dg_local_precomputed_values_set_element (problem->local_values_adapt, &geometry_data, problem->global_values);
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
    t8dg_local_precomputed_values_new (num_elems, problem->dim, problem->max_num_element_values, problem->max_num_faces,
                                       problem->max_num_face_values);
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
  t8dg_mortar_array_destroy (&problem->face_mortars);
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
  problem->face_mortars =
    t8dg_mortar_array_new_empty (problem->forest, t8dg_global_precomputed_values_get_num_faces (problem->global_values));
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
  local_values_partition =
    t8dg_local_precomputed_values_new (num_local_elems_new, problem->dim, problem->max_num_element_values, problem->max_num_faces,
                                       problem->max_num_face_values);
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
  t8dg_mortar_array_destroy (&problem->face_mortars);
  t8_forest_unref (&problem->forest);
  problem->forest = forest_partition;
  problem->face_mortars =
    t8dg_mortar_array_new_empty (problem->forest, t8dg_global_precomputed_values_get_num_faces (problem->global_values));
  t8_debugf (" [ADVECT] Done partition\n");
}
