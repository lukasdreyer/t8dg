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
#include <t8_forest/t8_forest_ghost.h>

#include <sc_containers.h>

#include "t8dg.h"
#include "t8dg_coarse_geometry.h"
#include "t8dg_advect_problem.h"
#include "t8dg_flux.h"
#include "t8dg_global_precomputed_values.h"
#include "t8dg_local_precomputed_values.h"
#include "t8dg_quadrature.h"
#include "t8dg_functionbasis.h"
#include "t8dg_sc_array.h"
#include "t8dg_timestepping.h"
#include "t8dg_mortar.h"

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

/* The following sc_arrays contain data, that needs to be available for each processorlocal element. To facilitate partitioning, they are
 * saved in a linear array
 */
  t8dg_local_precomputed_values_t *local_values;
  t8dg_local_precomputed_values_t *local_values_adapt;

  /* The dof_values get ghosted */
  sc_array_t         *element_dof_values;       /**< The Value of u at the nodal basis vertices */
  sc_array_t         *element_dof_values_adapt;

  /* those need to be recalculated for each time step, remain processor local */
  sc_array_t         *face_mortar[MAX_FACES];                   /**< contains pointer to face_mortars, so that fluxes need only be calculated once */

};

/* Getter functions for the sc_array_contained in advect_linear_problem.
 * TODO: setter functions
 *   */

static double      *
t8dg_advect_element_get_element_dof_values (const t8dg_linear_advection_problem_t * problem, t8_locidx_t idata)
{
  return ((double *) t8_sc_array_index_locidx (problem->element_dof_values, idata));
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
t8dg_advect_problem_endtime_reached (t8dg_linear_advection_problem_t * problem)
{
  return t8dg_timestepping_data_is_endtime_reached (problem->time_data);
}

void
t8dg_advect_problem_printdof (t8dg_linear_advection_problem_t * problem)
{
  t8dg_sc_array_block_double_print (problem->element_dof_values);
}

static void
t8dg_element_set_dofs_initial (t8dg_linear_advection_problem_t * problem, t8_locidx_t itree, t8_eclass_scheme_c * scheme,
                               t8_locidx_t ielement, t8_element_t * element)
{
  int                 idof;
  double             *element_dof_values;
  double              reference_vertex[DIM3];
  double              coarse_vertex[DIM3];
  double              image_vertex[DIM3];
  double              start_time = t8dg_timestepping_data_get_current_time (problem->time_data);
  t8_locidx_t         idata;

  idata = t8dg_itree_ielement_to_idata (problem->forest, itree, ielement);

  t8dg_functionbasis_t *functionbasis = t8dg_global_precomputed_values_get_functionbasis (problem->global_values);

  element_dof_values = t8dg_advect_element_get_element_dof_values (problem, idata);

  for (idof = 0; idof < t8dg_functionbasis_get_num_dof (functionbasis); idof++) {
    /*TODO: make available for general functionbasis */
    /* get_basisfunction_nodal vertex */
    t8dg_functionbasis_get_vertex (reference_vertex, functionbasis, idof);
    t8_debugf ("reference_vertex\n");
    t8dg_vec_print (reference_vertex);
    /* transform into coarse reference element */
    t8dg_local_precomputed_values_fine_to_coarse_geometry (reference_vertex, coarse_vertex, scheme, element);
    t8_debugf ("coarse_vertex\n");
    t8dg_vec_print (coarse_vertex);

    /* tree vertices are application data for linear geometry */
    problem->coarse_geometry->geometry (coarse_vertex, image_vertex, problem->forest, itree);
    t8_debugf ("image_vertex\n");
    t8dg_vec_print (image_vertex);

    /* apply initial condition function at image vertex and start time */
    element_dof_values[idof] = problem->description.initial_condition_fn (image_vertex, start_time);
  }
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
          t8dg_mortar_fill (mortar, problem->element_dof_values, problem->time_data,
                            problem->global_values, problem->local_values,
                            problem->description.flux, problem->description.numerical_flux_fn);
        }
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
                          int level,
                          int number_LGL_points, double start_time, double end_time, double cfl, int time_order, sc_MPI_Comm comm)
{
  t8dg_linear_advection_problem_t *problem;
  t8_scheme_cxx_t    *default_scheme;
  int                 iface;
  int                 num_elements;
  /* allocate problem */
  problem = T8_ALLOC (t8dg_linear_advection_problem_t, 1);

  default_scheme = t8_scheme_new_default_cxx ();
  t8_debugf ("create uniform forest\n");
  problem->forest = t8_forest_new_uniform (cmesh, default_scheme, level, 1, comm);

  problem->dim = dim;
  problem->uniform_refinement_level = level;

  problem->coarse_geometry = coarse_geometry;

  double              tangential_vector[3] = { 1, 0, 0 };       /*TODO: make dependent on cmesh! */
  problem->description.initial_condition_fn = u_initial;
  problem->description.flux = t8dg_linear_flux_new_1D_linear_geometry (tangential_vector, flow_speed);
  problem->description.numerical_flux_fn = t8dg_linear_numerical_flux_upwind_1D;

  problem->time_data = t8dg_timestepping_data_new (time_order, start_time, end_time, cfl);
  t8dg_timestepping_data_set_time_step (problem->time_data, cfl * pow (2, -level));     /* TODO: make dependent on cfl number and element diameter */

  problem->vtk_count = 0;
  problem->comm = comm;

  t8_debugf ("precompute global values\n");
  /* these allocate memory: */
  problem->global_values = t8dg_global_precomputed_values_new_1D_LGL (number_LGL_points);

  t8dg_quadrature_t  *quadrature = t8dg_global_precomputed_values_get_quadrature (problem->global_values);

  num_elements = t8_forest_get_num_element (problem->forest);

  problem->local_values = t8dg_local_precomputed_values_new (quadrature, num_elements);
  problem->local_values_adapt = NULL;

  /*currently no ghost, since serial, but generally the dof_values need to be ghosted. */
  problem->element_dof_values =
    sc_array_new_count (sizeof (double) * t8dg_global_precomputed_values_get_num_dof (problem->global_values),
                        num_elements + t8_forest_get_num_ghosts (problem->forest));

  problem->element_dof_values_adapt = NULL;

  /*for each element and face a pointer to a mortar */
  for (iface = 0; iface < t8dg_quadrature_get_num_faces (quadrature); iface++) {
    problem->face_mortar[iface] = sc_array_new_count (sizeof (t8dg_mortar_t *), num_elements);
  }

  t8_debugf ("finished problem init\n");

  return problem;
}

/*TODO: which init function creates what, outsource problem description*/
t8dg_linear_advection_problem_t *
t8dg_advect_problem_init_linear_geometry_1D (t8_cmesh_t cmesh,
                                             t8dg_scalar_function_3d_time_fn u_initial,
                                             double flow_speed,
                                             int level,
                                             int number_LGL_points,
                                             double start_time, double end_time, double cfl, int time_order, sc_MPI_Comm comm)
{

  t8dg_coarse_geometry_t *coarse_geometry = t8dg_coarse_geometry_new_1D_linear ();
  return t8dg_advect_problem_init (cmesh, coarse_geometry, 1, u_initial, flow_speed, level, number_LGL_points, start_time, end_time, cfl,
                                   time_order, comm);
}

void
t8dg_advect_problem_init_elements (t8dg_linear_advection_problem_t * problem)
{
  t8_locidx_t         itree, ielement, idata;
  t8_locidx_t         num_trees, num_elems_in_tree;
  t8_element_t       *element;

  int                 iface;

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

      for (iface = 0; iface < MAX_FACES; iface++) {
        t8dg_advect_element_set_face_mortar (problem, idata, iface, NULL);
      }
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
  sc_array_destroy (problem->element_dof_values);

  /* Free the arrays */
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
      t8dg_local_precomputed_values_element_multiply_trafo_quad_weight (problem->local_values, element_quad_values, idata);
      t8dg_flux_element_multiply_flux_value (problem->description.flux, element_quad_values,
                                             t8dg_timestepping_data_get_current_time (problem->time_data), problem->local_values,
                                             problem->forest, itree, ielement,
                                             t8dg_global_precomputed_values_get_quadrature (problem->global_values),
                                             problem->coarse_geometry);
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
t8dg_advect_problem_apply_inverse_mass_matrix_inplace (t8dg_linear_advection_problem_t * problem, sc_array_t * dof_array)
{
  /* TODO: include Vandermonde matrix.
   * Not inplace, but allow src_dof = dest_dof
   **/
  t8_locidx_t         itree, ielement, idata;
  t8_locidx_t         num_trees, num_elems_in_tree;
  sc_array_t         *element_dof_array;

  num_trees = t8_forest_get_num_local_trees (problem->forest);

  for (itree = 0, idata = 0; itree < num_trees; itree++) {

    num_elems_in_tree = t8_forest_get_tree_num_elements (problem->forest, itree);

    for (ielement = 0; ielement < num_elems_in_tree; ielement++, idata++) {
      element_dof_array = t8dg_sc_array_block_double_new_view (dof_array, idata);
      t8dg_local_precomputed_values_element_divide_trafo_quad_weight (problem->local_values, element_dof_array, idata);
      sc_array_destroy (element_dof_array);
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
  t8_debugf ("start calculating time derivate\n");
  T8DG_ASSERT (application_data != NULL);
  t8dg_linear_advection_problem_t *problem = (t8dg_linear_advection_problem_t *) application_data;
  T8DG_ASSERT (dof_values == problem->element_dof_values);
  T8DG_ASSERT (t == t8dg_timestepping_data_get_current_time (problem->time_data));

  sc_array_t         *dof_flux;
  dof_flux = t8dg_sc_array_duplicate (dof_change);
  t8_debugf ("test time derivate\n");

  t8dg_advect_problem_apply_stiffness_matrix (problem, problem->element_dof_values, dof_change);

  t8_debugf ("A u\n");
  t8dg_sc_array_block_double_debug_print (dof_change);

  /*Ghost exchange */
  t8_forest_ghost_exchange_data (problem->forest, problem->element_dof_values);

  t8dg_advect_problem_mortars_fill (problem);
  t8dg_advect_problem_apply_boundary_integrals (problem, dof_flux);
  t8dg_advect_problem_mortars_invalidate (problem);

  t8dg_sc_array_block_double_axpy (-1, dof_flux, dof_change);
  sc_array_destroy (dof_flux);

  t8_debugf ("A u  - B u\n");
  t8dg_sc_array_block_double_debug_print (dof_change);

  /*apply massinverse */
  t8dg_advect_problem_apply_inverse_mass_matrix_inplace (problem, dof_change);

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
  t8dg_timestepping_runge_kutta_step (t8dg_advect_time_derivative, t8dg_advect_get_time_data (problem),
                                      &(problem->element_dof_values), problem);
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
    dof_values = t8dg_advect_element_get_element_dof_values (problem, idata);
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

}

static void
t8dg_advect_problem_adapt (t8dg_linear_advection_problem_t * problem)
{
  /* Nothing to do */
  if (problem->maximum_refinement_level - problem->uniform_refinement_level == 0)
    return;

  t8_locidx_t         num_elems_p_ghosts, num_elems;
  t8_locidx_t         ghost_sent;
  t8_forest_t         forest_adapt;

  /* Adapt the forest, but keep the old one */
  t8_forest_ref (problem->forest);
  t8_forest_init (&forest_adapt);
  /* Set the user data pointer of the new forest */
  t8_forest_set_user_data (forest_adapt, problem);
  /* Set the adapt function */
  t8_forest_set_adapt (forest_adapt, problem->forest, t8dg_advect_test_adapt, 0);
  if (problem->maximum_refinement_level - problem->uniform_refinement_level > 1) {
    /* We also want to balance the forest if there is a possibility of elements
     * with difference in refinement levels greater 1 */
    t8_forest_set_balance (forest_adapt, NULL, 1);
//    did_balance = 1;
  }
  /* We also want ghost elements in the new forest */
  t8_forest_set_ghost (forest_adapt, 1, T8_GHOST_FACES);
  /* Commit the forest, adaptation and balance happens here */
  t8_forest_commit (forest_adapt);

  /* Allocate new memory for the element_data of the advected forest */
  num_elems = t8_forest_get_num_element (forest_adapt);
  num_elems_p_ghosts = num_elems + t8_forest_get_num_ghosts (forest_adapt);
#if 0
  problem->element_data_adapt = sc_array_new_count (sizeof (t8_advect_element_data_t), num_elems);
  problem->phi_values_adapt = sc_array_new_count ((problem->dummy_op ? 2 : 1) * sizeof (double), num_elems_p_ghosts);
#endif

  /* We now call iterate_replace in which we interpolate the new element data.
   * It is necessary that the old and new forest only differ by at most one level.
   * We guarantee this by calling adapt non-recursively and calling balance without
   * repartitioning. */
  t8_forest_iterate_replace (problem->forest_adapt, problem->forest, t8_advect_replace);

  /* clean the old element data */
  t8_advect_problem_elements_destroy (problem);
  sc_array_destroy (problem->element_data);
  sc_array_destroy (problem->phi_values);
  /* Free memory for the forest */
  t8_forest_unref (&problem->forest);
  /* Set the forest to the adapted one */
  problem->forest = problem->forest_adapt;
  problem->forest_adapt = NULL;
  /* Set the elem data to the adapted elem data */
  problem->element_data = problem->element_data_adapt;
  problem->element_data_adapt = NULL;
  /* Set the phi values to the adapted phi values */
  problem->phi_values = problem->phi_values_adapt;
  problem->phi_values_adapt = NULL;
}
#endif
