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

struct t8dg_advect_problem_linear_1D
{
  int                 dim = 1;

/*TODO: add source term functions*/
  t8dg_scalar_function_3d_time_fn initial_condition_fn;
  t8dg_scalar_function_3d_time_fn source_sink_fn;

  t8dg_linear_flux_velocity_3D_time_fn flux_velocity_fn;
  t8dg_linear_numerical_flux_1D_fn numerical_flux_fn;

  t8_forest_t         forest;

  /*Once per element */
  sc_array_t         *element_fine_to_coarse_geometry_data;     /*those get partitioned */

  sc_array_t         *element_trafo_quad_weight;        /* Q */
  sc_array_t         *face_trafo_quad_weight[MAX_FACES];        /* FQ */

  sc_array_t         *face_normal_vectors[MAX_FACES];   /* For each face and quadrature point a 3-dim vector */
  sc_array_t         *element_transformed_gradient_tangential_vectors;  /* For each */

  sc_array_t         *element_dof_values;       /*those get ghosted */
  sc_array_t         *element_dof_values_new;   /*those are only needed locally to save the result of a rungekutta timestep */

  sc_array_t         *face_mortar[MAX_FACES];   /*those need to be recalculated for each time step, remain processor local */

  int                 uniform_refinement_level; /*uniform refinement level */

  int                 time_order;

  double              delta_t;  /*time step */
  double              t;        /*current time */
  double              T;        /*end time */
  double              cfl;

  t8dg_LGL_quadrature_t *quadrature;
  t8dg_LGL_functionbasis_t *functionbasis;
  t8dg_coarse_geometry_3D_t *coarse_geometry;
  t8dg_time_matrix_application evolution_matrix;

  sc_MPI_Comm         comm;
};

/*Additional functions in t8dg_sc_array that return a view on the values for an element  */

static double      *
t8dg_advect_element_get_element_dof_values (const t8dg_advect_problem_linear_1D_t * problem, t8_locidx_t ielement)
{
  return ((double *)
          t8_sc_array_index_locidx (problem->element_dof_values, ielement));
}

static double      *
t8dg_advect_element_get_face_quad_trafo_weights (const t8dg_advect_problem_linear_1D_t * problem, t8_locidx_t ielement, int faceindex)
{
  return ((double *)
          t8_sc_array_index_locidx (problem->face_trafo_quad_weight[faceindex], ielement));
}

static double      *
t8dg_advect_element_get_element_quad_trafo_weights (const t8dg_advect_problem_linear_1D_t * problem, t8_locidx_t ielement)
{
  return ((double *)
          t8_sc_array_index_locidx (problem->element_trafo_quad_weight, ielement));
}

/*Assumes only one tangential vector per quad point*/
static double      *
t8dg_element_quad_get_transformed_gradient_tangential_vector (t8dg_advect_problem_linear_1D_t * problem,
                                                              t8_locidx_t ielement, t8dg_quad_idx_t iquad)
{
  T8_ASSERT (iquad >= 0 && (size_t) iquad < problem->element_transformed_gradient_tangential_vectors->elem_size / (3 * sizeof (double)));
  return ((double *) t8_sc_array_index_locidx (problem->element_transformed_gradient_tangential_vectors, ielement)) + 3 * iquad;
}

static double      *
t8dg_element_quad_get_normal_vector (t8dg_advect_problem_linear_1D_t * problem, t8_locidx_t ielement, int iface, t8dg_quad_idx_t iquad)
{
  T8_ASSERT (iface >= 0 && iface < MAX_FACES);
  T8_ASSERT (iquad >= 0 && (size_t) iquad < problem->face_normal_vectors[iface]->elem_size / (3 * sizeof (double)));
  return ((double *) t8_sc_array_index_locidx (problem->face_normal_vectors[iface], ielement)) + 3 * iquad;
}

/*  get functions for structs at element and faces: */

static t8dg_element_fine_to_coarse_geometry_data_t *
t8dg_advect_element_get_fine_to_coarse_geometry_data (const t8dg_advect_problem_linear_1D_t * problem, t8_locidx_t ielement)
{
  return ((t8dg_element_fine_to_coarse_geometry_data_t *)
          t8_sc_array_index_locidx (problem->element_fine_to_coarse_geometry_data, ielement));
}

#if 0
static t8dg_mortar_t *t8dg_advect_element_get_face_mortar ();
#endif

static void
t8dg_element_set_dofs_initial (t8dg_advect_problem_linear_1D_t * problem, t8_locidx_t ielement, double *tree_vertices)
{
  int                 idof;
  double             *element_dof_values;
  double              reference_vertex[DIM3];
  double              coarse_vertex[DIM3];
  double              image_vertex[DIM3];

  t8dg_element_fine_to_coarse_geometry_data_t *element_geometry_data;

  element_geometry_data = t8dg_advect_element_get_fine_to_coarse_geometry_data (problem, ielement);
  element_dof_values = t8dg_advect_element_get_element_dof_values (problem, ielement);

  for (idof = 0; idof < t8dg_LGL_functionbasis_get_num_dof (problem->functionbasis); idof++) {
    /*get_basisfunction_nodal vertex */
    t8dg_LGL_functionbasis_get_vertex (reference_vertex, problem->functionbasis, idof);

    t8dg_fine_to_coarse_geometry (coarse_vertex, reference_vertex, element_geometry_data);

    /* tree vertices are application data for linear geometry */
    problem->coarse_geometry->geometry (image_vertex, coarse_vertex, tree_vertices);
    T8_ASSERT (element_dof_values != NULL);

    element_dof_values[idof] = problem->initial_condition_fn (image_vertex, problem->t);
  }
}

static void
t8dg_element_set_precalculated_values_1D_linear (t8dg_advect_problem_linear_1D_t * problem, t8_locidx_t ielement, double *tree_vertices)
{
  T8_ASSERT (problem->dim == 1);
  double              coarse_1D_tangential_vector[DIM3];
  double             *transformed_gradient_tangential_vector;
  double             *normal_vector;
  double              element_normal_vector[DIM3];
  double              gram_det, h, norm;
  int                 iquad, iface;
  int                 num_faces;
  int                 num_elem_quad, num_face_quad;
  double             *element_quad_trafo = t8dg_advect_element_get_element_quad_trafo_weights (problem, ielement);
  double             *face_quad_trafo[MAX_FACES];

  num_faces = t8dg_LGL_quadrature_get_num_faces (problem->quadrature);
  num_elem_quad = t8dg_LGL_quadrature_get_num_element_vertices (problem->quadrature);
  /*tangential vector of the coarse element */
  t8_vec_axpyz (tree_vertices, tree_vertices + 3, coarse_1D_tangential_vector, -1);

  t8dg_element_fine_to_coarse_geometry_data_t *geometry_data;
  geometry_data = t8dg_advect_element_get_fine_to_coarse_geometry_data (problem, ielement);

  /*for more complicated geometries these values differ for different quadrature points */
  h = geometry_data->scaling_factor;    /*TODO: implement */
  norm = t8_vec_norm (coarse_1D_tangential_vector);
  gram_det = h * norm;

  t8_vec_axb (coarse_1D_tangential_vector, element_normal_vector, 1. / norm, 0);

  for (iquad = 0; iquad < num_elem_quad; iquad++) {

    transformed_gradient_tangential_vector = t8dg_element_quad_get_transformed_gradient_tangential_vector (problem, ielement, iquad);
    t8_vec_axb (coarse_1D_tangential_vector, transformed_gradient_tangential_vector, 1. / (h * norm * norm), 0);

    element_quad_trafo[iquad] = gram_det * t8dg_LGL_quadrature_get_element_weight (problem->quadrature, iquad);
  }

  for (iface = 0; iface < num_faces; iface++) {
    face_quad_trafo[iface] = t8dg_advect_element_get_face_quad_trafo_weights (problem, ielement, iface);
    num_face_quad = t8dg_LGL_quadrature_get_num_face_vertices (problem->quadrature, iface);
    for (iquad = 0; iquad < num_face_quad; iquad++) {
      normal_vector = t8dg_element_quad_get_normal_vector (problem, ielement, iface, iquad);

      switch (iquad) {
      case 0:
        t8_vec_axb (element_normal_vector, normal_vector, -1, 0);
        break;
      case 1:
        t8_vec_axb (element_normal_vector, normal_vector, 1, 0);        /* vec_copy */
        break;
      default:
        T8DG_ASSERT (0);
        break;
      }
      face_quad_trafo[iface][iquad] = 1;
    }
  }

}

static void
t8dg_advect_evolution (sc_array_t * dudt_array, const sc_array_t * u_array, double t, const void *application_data)
{
  /** TODO: Only for testing purposes!!*/
  T8_ASSERT (dudt_array->elem_count == u_array->elem_count);
  T8_ASSERT (dudt_array->elem_size == u_array->elem_size);
  unsigned            i;
  for (i = 0; i < dudt_array->elem_count * dudt_array->elem_size / sizeof (double); i++) {
    ((double *) dudt_array->array)[i] = (2. / t) * ((double *) u_array->array)[i];
  }
  return;

  /*In reality du/dt = invMassmatrix((A_c)u - Bu + Mg) */
#if 0
  t8dg_advect_problem_linear_1D_t *problem = (t8dg_advect_problem_linear_1D_t *) application_data;
  sc_array_t         *stiffness_u = sc_array_new_count (u_array->elem_size, u_array->elem_count);
  sc_array_t         *boundary_u = sc_array_new_count (u_array->elem_size, u_array->elem_count);
//  sc_array_t * mass_g = sc_array_new_count(u_array->elem_size, u_array->elem_count);
  apply_stiffness (stiffness_u, u_array, t, problem);
  apply_boundary (boundary_u, u_array, t, problem);
//  apply_mass(mass_g,g_array,problem);
  t8dg_sc_array_block_double_axpy (1, boundary_u, stiffness_u);
//  t8dg_sc_array_block_double_axpy(1,mass_g,stiffness_u);
  apply_mass_invers (dudt_array, stiffness_u, problem);
#endif
}

t8dg_advect_problem_linear_1D_t *
t8dg_advect_problem_init_linear_1D (t8_cmesh_t cmesh, t8dg_scalar_function_3d_time_fn u_initial, double flow_velocity,
                                    int level, int number_LGL_points,
                                    double start_time, double end_time, double cfl, int time_order, sc_MPI_Comm comm)
{
  t8dg_advect_problem_linear_1D_t *problem;
  t8_scheme_cxx_t    *default_scheme;
  int                 iface;

  /* allocate problem */
  problem = T8_ALLOC (t8dg_advect_problem_linear_1D_t, 1);

  problem->dim = 1;
  problem->uniform_refinement_level = level;
  problem->comm = comm;

  problem->initial_condition_fn = u_initial;

  problem->time_order = time_order;
  problem->T = end_time;
  problem->t = start_time;
  problem->cfl = cfl;
  problem->delta_t = 0.1;       /* TODO: make dependent on cfl number and element diameter */

  t8_debugf ("start LGL construction\n");
  /* these allocate memory: */
  t8dg_LGL_quadrature_and_functionbasis_new_1D (&problem->quadrature, &problem->functionbasis, number_LGL_points);
  problem->coarse_geometry = t8dg_coarse_geometry_new_1D_linear ();

  problem->numerical_flux_fn = t8dg_upwind_flux_1D;
  problem->evolution_matrix = t8dg_advect_evolution;

  default_scheme = t8_scheme_new_default_cxx ();
  t8_debugf ("create uniform forest\n");

  problem->forest = t8_forest_new_uniform (cmesh, default_scheme, level, 1, comm);

  int                 num_elements = t8_forest_get_num_element (problem->forest);

  t8_debugf ("start creating sc_arrays\n");

  problem->element_fine_to_coarse_geometry_data = sc_array_new_count (sizeof (t8dg_element_fine_to_coarse_geometry_data_t), num_elements);

  problem->element_trafo_quad_weight =
    sc_array_new_count (sizeof (double) * t8dg_LGL_quadrature_get_num_element_vertices (problem->quadrature), num_elements);

  for (iface = 0; iface < t8dg_LGL_quadrature_get_num_faces (problem->quadrature); iface++) {
    problem->face_trafo_quad_weight[iface] =
      sc_array_new_count (sizeof (double) * t8dg_LGL_quadrature_get_num_face_vertices (problem->quadrature, iface), num_elements);
    problem->face_mortar[iface] =
      sc_array_new_count (sizeof (double) * t8dg_LGL_quadrature_get_num_face_vertices (problem->quadrature, iface), num_elements);

    problem->face_normal_vectors[iface] =
      sc_array_new_count (sizeof (double) * 3 * t8dg_LGL_quadrature_get_num_element_vertices (problem->quadrature), num_elements);

  }

  problem->element_transformed_gradient_tangential_vectors =
    sc_array_new_count (sizeof (double) * 3 * t8dg_LGL_quadrature_get_num_element_vertices (problem->quadrature), num_elements);

  /*currently no ghost, since serial, but generally the dof_values need to be ghosted. */
  problem->element_dof_values =
    sc_array_new_count (sizeof (double) * t8dg_LGL_functionbasis_get_num_dof (problem->functionbasis),
                        num_elements + t8_forest_get_num_ghosts (problem->forest));

  problem->element_dof_values_new =
    sc_array_new_count (sizeof (double) * t8dg_LGL_functionbasis_get_num_dof (problem->functionbasis), num_elements);

  t8_debugf ("finished creating sc_arrays\n");

  return problem;
}

void
t8dg_advect_problem_destroy (t8dg_advect_problem_linear_1D_t ** pproblem)
{
  t8dg_advect_problem_linear_1D_t *problem;
  int                 iface;

  T8_ASSERT (pproblem != NULL);
  problem = *pproblem;
  if (problem == NULL) {
    return;
  }

  problem->dim = -1;

  /* Free the arrays */
  sc_array_destroy (problem->element_fine_to_coarse_geometry_data);
  sc_array_destroy (problem->element_dof_values);
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
t8dg_advect_problem_init_elements_linear_1D (t8dg_advect_problem_linear_1D_t * problem)
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

      t8dg_element_set_geometry_data (geometry_data, element, idata, scheme);

      /*precompute values for element idata */

      t8dg_element_set_dofs_initial (problem, idata, tree_vertices);

      t8dg_element_set_precalculated_values_1D_linear (problem, idata, tree_vertices);
    }
  }
  t8_debugf ("End element init \n");
}

void
t8dg_advect_evolve (t8dg_advect_problem_linear_1D_t * problem)
{
//  printf("time %f\n",problem->t);
  if (problem->t + problem->delta_t > problem->T) {
    problem->delta_t = problem->T - problem->t;
  }

  /*TODO: calculate fluxes in mortars */

  t8dg_rungekutta_timestep (problem->time_order, problem->t, problem->delta_t, problem->evolution_matrix,
                            problem->element_dof_values_new, problem->element_dof_values, NULL);
  problem->t += problem->delta_t;
  /*TODO: swap problem.dof_values and problem.dof_new */
  t8dg_sc_array_swap (&problem->element_dof_values, &problem->element_dof_values_new);

}

int
t8dg_advect_problem_endtime_reached (t8dg_advect_problem_linear_1D_t * problem)
{
  return problem->t >= problem->T;
}

void
t8dg_advect_problem_printdof (t8dg_advect_problem_linear_1D_t * problem)
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

  t8dg_advect_problem_linear_1D_t *problem = (t8dg_advect_problem_linear_1D_t *) application_data;
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

  t8dg_advect_problem_linear_1D_t *problem = (t8dg_advect_problem_linear_1D_t *) application_data;
  for (idata = 0; idata < num_data; idata++) {
    element_dof_values = t8dg_sc_array_new_double_element_view (src);
    element_dof_values_new = t8dg_sc_array_new_double_element_view (dest);

    sc_array_destroy (element_dof_values);
    sc_array_destroy (element_dof_values_new);
  }
}

void                t8dg_apply_boundary_fluxes (sc_array_t * dest, const sc_array_t * src, const void *application_data);
#endif
