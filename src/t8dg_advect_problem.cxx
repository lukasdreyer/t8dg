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

#include <sc_containers.h>

#include "t8dg.h"
#include "t8dg_geometry.h"
#include "t8dg_advect_problem.h"
#include "t8dg_numerical_flux.h"
#include "t8dg_LGL.h"
#include "t8dg_sc_array.h"
#include "t8dg_timestepping.h"

/** The container for all information needed to solve the advection problem.
 */
struct t8dg_linear_advection_problem
{
  int                 dim;      /**< Dimension of the submanifold */

/*TODO: add source term functions*/
  t8dg_scalar_function_3d_time_fn initial_condition_fn;         /**< Initial condition function */
#if 0
  t8dg_scalar_function_3d_time_fn source_sink_fn;
#endif
  t8dg_linear_flux_velocity_3D_time_fn flux_velocity_fn;        /**< Divergence free vector field */
  t8dg_linear_numerical_flux_1D_fn numerical_flux_fn;           /**< Approximation to the Riemann problem */

  t8_forest_t         forest;                                   /**< The t8_forest used to adapt the coarse mesh*/

/* The following sc_arrays contain data, that needs to be available for each processorlocal element. To facilitate partitioning, they are
 * saved in a linear array
 */
  sc_array_t         *element_fine_to_coarse_geometry_data;     /**< For each element, the data to determine the linear function from the fine reference
								element into the coarse reference element */
  sc_array_t         *element_trafo_quad_weight;                /**< For each element, an array of length num_quad_points with the combined weight
								of quadrature and transformation */
  sc_array_t         *element_transformed_gradient_tangential_vectors;  /**< For each element, (d_x F^-1)(tau_d), needed to calculate the gradient */

  /* The dof_values get ghosted */
  sc_array_t         *element_dof_values;       /**< The Value of u at the nodal basis vertices */
  sc_array_t         *element_dof_change;       /**< time-derivative used in the rungekutta steps to advance the solution */
  sc_array_t         *element_dof_values_new;   /**< partial sum of the runge-kutta timestepping, at the end
  of a rk timestep use to overwrite dof_values */

  /* To avoid another dimension when flattening the array, and since the number of Faces is bound, have an sc_array for each faceindex of length
   * of the processorlocal elements*/
  sc_array_t         *face_trafo_quad_weight[MAX_FACES];        /**< for each element, an array of length num_quad_points with the combined weight
								of quadrature and transformation */
  sc_array_t         *face_normal_vectors[MAX_FACES];           /**< For each face and quadrature point the 3-dim normal vector in image space*/

  /* those need to be recalculated for each time step, remain processor local */
  sc_array_t         *face_mortar[MAX_FACES];                   /**< contains pointer to face_mortars, so that fluxes need only be calculated once */

  int                 uniform_refinement_level; /**< uniform refinement level */

  int                 time_order;/**< time order of the Runge-kutta timestepping*/

  double              delta_t;  /**< time step */
  double              t;        /**< current time */
  double              T;        /**< end time */
  double              cfl;      /**< cfl number*/

  t8dg_LGL_quadrature_t *quadrature;            /**< LGL quadrature, based on same vertex-set as functionbasis*/
  t8dg_LGL_functionbasis_t *functionbasis;      /**< LGL functionbasis, based on same vertex-set as quadrature */
  t8dg_coarse_geometry_3D_t *coarse_geometry;   /**< coarse geometry, that gives the geometry and Jacobian of the geometry from the coarse
						 reference element to the coarse image element*/

  sc_MPI_Comm         comm; /**< MPI Communicator */
};

/* Additional functions in t8dg_sc_array that return a view on the values for an element  */

static double      *
t8dg_advect_element_get_element_dof_values (const t8dg_linear_advection_problem_t * problem, t8_locidx_t idata)
{
  return ((double *)
          t8_sc_array_index_locidx (problem->element_dof_values, idata));
}

static double      *
t8dg_advect_element_get_face_quad_trafo_weights (const t8dg_linear_advection_problem_t * problem, t8_locidx_t idata, int faceindex)
{
  return ((double *)
          t8_sc_array_index_locidx (problem->face_trafo_quad_weight[faceindex], idata));
}

static double      *
t8dg_advect_element_get_element_quad_trafo_weights (const t8dg_linear_advection_problem_t * problem, t8_locidx_t idata)
{
  return ((double *)
          t8_sc_array_index_locidx (problem->element_trafo_quad_weight, idata));
}

/*Assumes only one tangential vector per quad point, TODO: add getter function for 3D vectors in sc_arrays*/
static double      *
t8dg_element_quad_get_transformed_gradient_tangential_vector (t8dg_linear_advection_problem_t * problem,
                                                              t8_locidx_t idata, t8dg_quad_idx_t iquad)
{
  T8_ASSERT (problem->dim == 1);
  T8_ASSERT (iquad >= 0 && (size_t) iquad < problem->element_transformed_gradient_tangential_vectors->elem_size / (3 * sizeof (double)));
  return ((double *) t8_sc_array_index_locidx (problem->element_transformed_gradient_tangential_vectors, idata)) + 3 * iquad;
}

static double      *
t8dg_element_quad_get_normal_vector (t8dg_linear_advection_problem_t * problem, t8_locidx_t idata, int iface, t8dg_quad_idx_t iquad)
{
  T8_ASSERT (iface >= 0 && iface < MAX_FACES);
  T8_ASSERT (iquad >= 0 && (size_t) iquad < problem->face_normal_vectors[iface]->elem_size / (3 * sizeof (double)));
  return ((double *) t8_sc_array_index_locidx (problem->face_normal_vectors[iface], idata)) + 3 * iquad;
}

/*  get functions for structs at element and faces: */

static t8dg_element_fine_to_coarse_geometry_data_t *
t8dg_advect_element_get_fine_to_coarse_geometry_data (const t8dg_linear_advection_problem_t * problem, t8_locidx_t idata)
{
  return ((t8dg_element_fine_to_coarse_geometry_data_t *)
          t8_sc_array_index_locidx (problem->element_fine_to_coarse_geometry_data, idata));
}

#if 0
static t8dg_mortar_t *t8dg_advect_element_get_face_mortar ();
#endif

static void
t8dg_element_set_dofs_initial (t8dg_linear_advection_problem_t * problem, t8_locidx_t idata, double *tree_vertices)
{
  int                 idof;
  double             *element_dof_values;
  double              reference_vertex[DIM3];
  double              coarse_vertex[DIM3];
  double              image_vertex[DIM3];

  t8dg_element_fine_to_coarse_geometry_data_t *element_geometry_data;

  element_geometry_data = t8dg_advect_element_get_fine_to_coarse_geometry_data (problem, idata);
  element_dof_values = t8dg_advect_element_get_element_dof_values (problem, idata);

  for (idof = 0; idof < t8dg_LGL_functionbasis_get_num_dof (problem->functionbasis); idof++) {
    /* get_basisfunction_nodal vertex */
    t8dg_LGL_functionbasis_get_vertex (reference_vertex, problem->functionbasis, idof);

    /* transform into coarse reference element */
    t8dg_fine_to_coarse_geometry (coarse_vertex, reference_vertex, element_geometry_data);

    /* tree vertices are application data for linear geometry */
    problem->coarse_geometry->geometry (image_vertex, coarse_vertex, tree_vertices);

    /* apply initial condition function at image vertex and start time */
    element_dof_values[idof] = problem->initial_condition_fn (image_vertex, problem->t);
  }
}

static void
t8dg_element_set_precalculated_values_1D_linear (t8dg_linear_advection_problem_t * problem, t8_locidx_t idata, double *tree_vertices)
{
  T8_ASSERT (problem->dim == 1);
  double              coarse_1D_tangential_vector[DIM3];        /* vector between endvertices of the line in the image */

  double              gram_det, h, norm_tangential_vector;
  int                 iquad, iface;
  int                 num_faces, num_elem_quad, num_face_quad;

  t8dg_element_fine_to_coarse_geometry_data_t *geometry_data;

  /*pointer on the values to fill */
  double             *element_quad_trafo;       /*size: num_elem_quad */
  double             *face_quad_trafo;  /*size: num_face_quad */
  double             *transformed_gradient_tangential_vector;   /*size: 3, gets evaluated for each element quad point */
  double             *normal_vector;    /*size: 3, gets evaluated for each element quad point */

  num_faces = t8dg_LGL_quadrature_get_num_faces (problem->quadrature);
  num_elem_quad = t8dg_LGL_quadrature_get_num_element_vertices (problem->quadrature);

  /*tangential vector of the coarse element */
  t8_vec_axpyz (tree_vertices, tree_vertices + DIM3, coarse_1D_tangential_vector, -1);

  geometry_data = t8dg_advect_element_get_fine_to_coarse_geometry_data (problem, idata);
  /*for more complicated geometries these values differ for different quadrature points */
  h = geometry_data->scaling_factor;    /*TODO: implement */
  norm_tangential_vector = t8_vec_norm (coarse_1D_tangential_vector);
  gram_det = h * norm_tangential_vector;

  element_quad_trafo = t8dg_advect_element_get_element_quad_trafo_weights (problem, idata);
  for (iquad = 0; iquad < num_elem_quad; iquad++) {
    element_quad_trafo[iquad] = gram_det * t8dg_LGL_quadrature_get_element_weight (problem->quadrature, iquad);

    transformed_gradient_tangential_vector = t8dg_element_quad_get_transformed_gradient_tangential_vector (problem, idata, iquad);
    t8_vec_axb (coarse_1D_tangential_vector, transformed_gradient_tangential_vector,
                1. / (h * norm_tangential_vector * norm_tangential_vector), 0);
  }

  for (iface = 0; iface < num_faces; iface++) {
    face_quad_trafo = t8dg_advect_element_get_face_quad_trafo_weights (problem, idata, iface);
    num_face_quad = t8dg_LGL_quadrature_get_num_face_vertices (problem->quadrature, iface);
    for (iquad = 0; iquad < num_face_quad; iquad++) {
      /*for 1D elements the faceintegrals are just the value at the facequadrature point */
      face_quad_trafo[iquad] = 1;
      normal_vector = t8dg_element_quad_get_normal_vector (problem, idata, iface, iquad);
      /*scale the element tangential_vector to a unit vector in the right direction */
      switch (iface) {
      case 0:
        t8_vec_axb (coarse_1D_tangential_vector, normal_vector, -1. / norm_tangential_vector, 0);
        break;
      case 1:
        t8_vec_axb (coarse_1D_tangential_vector, normal_vector, 1. / norm_tangential_vector, 0);
        break;
      default:
        T8DG_ASSERT (0);
        break;
      }
    }
  }
}

t8dg_linear_advection_problem_t *
t8dg_advect_problem_init_linear_1D (t8_cmesh_t cmesh, t8dg_scalar_function_3d_time_fn u_initial, double flow_velocity,
                                    int level, int number_LGL_points,
                                    double start_time, double end_time, double cfl, int time_order, sc_MPI_Comm comm)
{
  t8dg_linear_advection_problem_t *problem;
  t8_scheme_cxx_t    *default_scheme;
  int                 iface;
  int                 num_elements;
  /* allocate problem */
  problem = T8_ALLOC (t8dg_linear_advection_problem_t, 1);

  problem->dim = 1;
  problem->uniform_refinement_level = level;
  problem->comm = comm;

  problem->initial_condition_fn = u_initial;

  problem->time_order = time_order;
  problem->T = end_time;
  problem->t = start_time;
  problem->cfl = cfl;
  problem->delta_t = cfl * 0.1; /* TODO: make dependent on cfl number and element diameter */

  t8_debugf ("start LGL construction\n");
  /* these allocate memory: */
  t8dg_LGL_quadrature_and_functionbasis_new_1D (&problem->quadrature, &problem->functionbasis, number_LGL_points);
  problem->coarse_geometry = t8dg_coarse_geometry_new_1D_linear ();

  problem->numerical_flux_fn = t8dg_upwind_flux_1D;

  default_scheme = t8_scheme_new_default_cxx ();
  t8_debugf ("create uniform forest\n");

  problem->forest = t8_forest_new_uniform (cmesh, default_scheme, level, 1, comm);

  num_elements = t8_forest_get_num_element (problem->forest);

  t8_debugf ("start creating sc_arrays\n");

  /*coarse geometry data for each local element */
  problem->element_fine_to_coarse_geometry_data = sc_array_new_count (sizeof (t8dg_element_fine_to_coarse_geometry_data_t), num_elements);

  /*for each element an array of double values */
  problem->element_trafo_quad_weight =
    sc_array_new_count (sizeof (double) * t8dg_LGL_quadrature_get_num_element_vertices (problem->quadrature), num_elements);

  for (iface = 0; iface < t8dg_LGL_quadrature_get_num_faces (problem->quadrature); iface++) {
    problem->face_trafo_quad_weight[iface] =
      sc_array_new_count (sizeof (double) * t8dg_LGL_quadrature_get_num_face_vertices (problem->quadrature, iface), num_elements);
    /*for each element and face a pointer to a mortar */
    problem->face_mortar[iface] = sc_array_new_count (sizeof (t8dg_mortar_t *), num_elements);

    problem->face_normal_vectors[iface] =
      sc_array_new_count (sizeof (double) * 3 * t8dg_LGL_quadrature_get_num_element_vertices (problem->quadrature), num_elements);

  }

  problem->element_transformed_gradient_tangential_vectors =
    sc_array_new_count (sizeof (double) * 3 * problem->dim * t8dg_LGL_quadrature_get_num_element_vertices (problem->quadrature),
                        num_elements);

  /*currently no ghost, since serial, but generally the dof_values need to be ghosted. */
  problem->element_dof_values =
    sc_array_new_count (sizeof (double) * t8dg_LGL_functionbasis_get_num_dof (problem->functionbasis),
                        num_elements + t8_forest_get_num_ghosts (problem->forest));

  problem->element_dof_change =
    sc_array_new_count (sizeof (double) * t8dg_LGL_functionbasis_get_num_dof (problem->functionbasis), num_elements);

  problem->element_dof_values_new =
    sc_array_new_count (sizeof (double) * t8dg_LGL_functionbasis_get_num_dof (problem->functionbasis), num_elements);

  t8_debugf ("finished creating sc_arrays\n");

  return problem;
}

void
t8dg_advect_problem_destroy (t8dg_linear_advection_problem_t ** pproblem)
{
  t8dg_linear_advection_problem_t *problem;
  int                 iface;

  T8_ASSERT (pproblem != NULL);
  problem = *pproblem;
  if (problem == NULL) {
    return;
  }

  /* Free the arrays */
  sc_array_destroy (problem->element_fine_to_coarse_geometry_data);
  sc_array_destroy (problem->element_dof_values);
  sc_array_destroy (problem->element_dof_change);
  sc_array_destroy (problem->element_dof_values_new);
  sc_array_destroy (problem->element_trafo_quad_weight);
  sc_array_destroy (problem->element_transformed_gradient_tangential_vectors);
  for (iface = 0; iface < t8dg_LGL_quadrature_get_num_faces (problem->quadrature); iface++) {
    sc_array_destroy (problem->face_trafo_quad_weight[iface]);
    sc_array_destroy (problem->face_mortar[iface]);
    sc_array_destroy (problem->face_normal_vectors[iface]);
  }

  t8dg_LGL_quadrature_and_functionbasis_destroy (&problem->quadrature, &problem->functionbasis);
  t8dg_coarse_geometry_destroy (&(problem->coarse_geometry));

  /* Unref the forest */
  t8_forest_unref (&problem->forest);   /*unrefs coarse mesh as well */
  /* Free the problem and set pointer to NULL */
  T8_FREE (problem);
  *pproblem = NULL;
}

void
t8dg_advect_problem_init_elements_linear_1D (t8dg_linear_advection_problem_t * problem)
{
  t8_locidx_t         itree, ielement, idata;
  t8_locidx_t         num_trees, num_elems_in_tree;
  t8_element_t       *element;

  t8dg_element_fine_to_coarse_geometry_data_t *geometry_data;

  t8_eclass_scheme_c *scheme;
  double             *tree_vertices;

  t8_debugf ("Start element init \n");
  num_trees = t8_forest_get_num_local_trees (problem->forest);

  for (itree = 0, idata = 0; itree < num_trees; itree++) {
    scheme = t8_forest_get_eclass_scheme (problem->forest, t8_forest_get_tree_class (problem->forest, itree));

    num_elems_in_tree = t8_forest_get_tree_num_elements (problem->forest, itree);

    tree_vertices = t8_forest_get_tree_vertices (problem->forest, itree);

    for (ielement = 0; ielement < num_elems_in_tree; ielement++, idata++) {
      element = t8_forest_get_element_in_tree (problem->forest, itree, ielement);

      geometry_data = t8dg_advect_element_get_fine_to_coarse_geometry_data (problem, idata);

      /*precompute values for the idata_th processorlocal element */

      t8dg_element_set_geometry_data (geometry_data, element, idata, scheme);

      t8dg_element_set_dofs_initial (problem, idata, tree_vertices);

      t8dg_element_set_precalculated_values_1D_linear (problem, idata, tree_vertices);
    }
  }
  t8_debugf ("End element init \n");
}

static void
t8dg_advect_calculate_time_derivative (t8dg_linear_advection_problem_t * problem)
{
  unsigned            i;
  for (i = 0; i < problem->element_dof_change->elem_count * problem->element_dof_change->elem_size / sizeof (double); i++) {
    ((double *) problem->element_dof_change->array)[i] = (2. / problem->t) * ((double *) problem->element_dof_values->array)[i];
  }

  /*updates dof_values_change */
  /*ghost dof_values */
  /*dudt = Mg + Au */
  /*receive ghosts */
  /*dudt -= Bu */
  /*dudt = M^-1 dudt */
  /*resize dof values to only internal values */
}

void
t8dg_advect_runge_kutta_step (t8dg_linear_advection_problem_t * problem)
{
  int                 istep;
  double             *rk_a, *rk_b, *rk_c;
  sc_array_t         *element_dof_beginning;
  double              time_beginning;

  /*Set already during construction? */
  if (problem->t + problem->delta_t > problem->T) {
    problem->delta_t = problem->T - problem->t;
  }

  t8dg_runge_kutta_fill_coefficients (problem->time_order, &rk_a, &rk_b, &rk_c);

  element_dof_beginning = t8dg_sc_array_clone (problem->element_dof_values);
  time_beginning = problem->t;

  t8dg_advect_calculate_time_derivative (problem);
   /**/
    t8dg_sc_array_block_double_zaxpy (problem->element_dof_values_new, rk_b[0] * problem->delta_t, problem->element_dof_change,
                                      problem->element_dof_values);

  for (istep = 0; istep < problem->time_order - 1; istep++) {
    /*calculate the y-value for which the derivative needs to be evaluated
     * since a has only values on the first minor diagonal, only the k from the step before and the original y is needed*/
    t8dg_sc_array_block_double_zaxpy (problem->element_dof_values, rk_a[istep] * problem->delta_t, problem->element_dof_change,
                                      element_dof_beginning);
    /* calculate the derivative at the step time and y value */

    problem->t = time_beginning + rk_c[istep] * problem->delta_t;
    t8dg_advect_calculate_time_derivative (problem);
     /**/
      /*add weighted summand to result */
      t8dg_sc_array_block_double_axpy (rk_b[istep + 1] * problem->delta_t, problem->element_dof_change, problem->element_dof_values_new);
  }
  problem->t = time_beginning + problem->delta_t;
  t8dg_sc_array_swap (&problem->element_dof_values, &problem->element_dof_values_new);

  sc_array_destroy (element_dof_beginning);

}

int
t8dg_advect_problem_endtime_reached (t8dg_linear_advection_problem_t * problem)
{
  return problem->t >= problem->T;
}

void
t8dg_advect_problem_printdof (t8dg_linear_advection_problem_t * problem)
{
  t8dg_sc_array_block_double_print (problem->element_dof_values);
}

#if 0
void
t8dg_apply_mass_matrix (sc_array_t * dest, const sc_array_t * src, const void *application_data)
{
  T8_ASSERT (application_data != NULL);
  const sc_array_t   *element_dof_values;
  sc_array_t         *element_dof_values_new;
  int                 idata;
  int                 num_data;

  t8dg_linear_advection_problem_t *problem = (t8dg_linear_advection_problem_t *) application_data;
  for (idata = 0; idata < num_data; idata++) {
    element_dof_values = t8dg_sc_array_new_double_element_view (src);
    element_dof_values_new = t8dg_sc_array_new_double_element_view (dest);
    t8dg_apply_trafo_quad_weights (element_dof_values_new, element_dof_values, idata, problem);
    sc_array_destroy (element_dof_values);
    sc_array_destroy (element_dof_values_new);
  }
}

void
t8dg_apply_stiffness_matrix (sc_array_t * dest, const sc_array_t * src, const void *application_data)
{
  T8_ASSERT (application_data != NULL);
  const sc_array_t   *element_dof_values;
  sc_array_t         *element_dof_values_new;
  int                 idata;
  int                 num_data;

  t8dg_linear_advection_problem_t *problem = (t8dg_linear_advection_problem_t *) application_data;
  for (idata = 0; idata < num_data; idata++) {
    element_dof_values = t8dg_sc_array_new_double_element_view (src);
    element_dof_values_new = t8dg_sc_array_new_double_element_view (dest);

    sc_array_destroy (element_dof_values);
    sc_array_destroy (element_dof_values_new);
  }
}

void                t8dg_apply_boundary_fluxes (sc_array_t * dest, const sc_array_t * src, const void *application_data);
#endif
