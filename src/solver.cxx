/*
 * solver.c
 *
 *  Created on: Mar 11, 2020
 *      Author: lukas
 */
#include <t8.h>
#include <example/common/t8_example_common.h>
#include <t8_schemes/t8_default_cxx.hxx>
#include <t8_forest.h>
#include "solver.hxx"
#include "global.h"
#include "sc_array_access.h"





static t8dg_1D_advect_problem_t *
t8dg_1D_advect_problem_init (t8_cmesh_t cmesh, t8_scalar_function_1d_fn u_0, double flow_velocity,
				   int level, int number_LGL_points, sc_MPI_Comm comm)
{
  t8dg_1D_advect_problem_t 		*problem;
  t8_scheme_cxx_t			*default_scheme;
  int					iface;


  /* allocate problem */
  problem = T8_ALLOC (t8dg_1D_advect_problem_t, 1);

  problem->dim = 1;
  problem->flow_velocity = flow_velocity;
  problem->quadrature = t8dg_1D_LGL_quadrature(number_LGL_points);/*allocates*/
  problem->functionbasis = t8dg_1DLGL_functionbasis(number_LGL_points);
  problem->u_0 = u_0;
  problem->comm = comm;

  problem->vandermonde = identity_matrix;
  problem->vandermonde_transpose = identity_matrix;
  problem->face_vandermonde = face_vandermonde_1D_linear_LGL;
  problem->face_vandermonde_transpose = face_vandermonde_transpose_1D_linear_LGL;
//  problem->directional_derivative_matrix = directional_derivative_1D_LGL2_matrix;

  problem->numerical_flux = upwind_flux_1D;


  default_scheme = t8_scheme_new_default_cxx ();
  problem->forest =
    t8_forest_new_uniform (cmesh, default_scheme, level, 1, comm);

  problem->element_values =
    sc_array_new_count (sizeof (t8dg_1D_advect_element_precomputed_values_t),
                        t8_forest_get_num_element (problem->forest));

  problem->advance_element_data =
    sc_array_new_count (sizeof (t8dg_1D_advect_advance_element_data_t),
                        t8_forest_get_num_element (problem->forest));

  problem->jacobian_at_quad =
    sc_array_new_count (sizeof (double) * problem->dim * problem->dim * number_LGL_points,
                        t8_forest_get_num_element (problem->forest));

  problem->element_trafo_quad_weight =
    sc_array_new_count (sizeof (double) * problem->quadrature->number_of_vertices,
                        t8_forest_get_num_element (problem->forest));

  for(iface=0; iface < problem->quadrature->number_of_faces; iface++){
      problem->face_trafo_quad_weight[iface] =
        sc_array_new_count (sizeof (double) * problem->quadrature->number_of_facevertices[iface],
                            t8_forest_get_num_element (problem->forest));
  }

  problem->dof_values =
    sc_array_new_count (sizeof (double) * problem->quadrature->number_of_vertices,
                        t8_forest_get_num_element (problem->forest) +
                        t8_forest_get_num_ghosts (problem->forest));

  return problem;
}

static void
t8dg_1D_advect_problem_destroy (t8dg_1D_advect_problem_t ** pproblem)
{
  t8dg_1D_advect_problem_t *problem;
  int iface;

  T8_ASSERT (pproblem != NULL);
  problem = *pproblem;
  if (problem == NULL) {
    return;
  }
  /* destroy elements */
//  t8_advect_problem_elements_destroy (problem);

  problem->dim = -1;
  problem->flow_velocity = 0;
  problem->u_0 = NULL;

  /* Free the arrays */
  sc_array_destroy (problem->element_values);
  sc_array_destroy (problem->dof_values);
  sc_array_destroy (problem->advance_element_data);
  sc_array_destroy (problem->element_trafo_quad_weight);
  sc_array_destroy (problem->jacobian_at_quad);
  for(iface = 0; iface <problem->quadrature->number_of_faces; iface++){
      sc_array_destroy (problem->face_trafo_quad_weight[iface]);
  }

  t8_dg_quadrature_destroy(&(problem->quadrature));


  /* Unref the forest */
  t8_forest_unref (&problem->forest);/*unrefs coarse mesh as well*/
  /* Free the problem and set pointer to NULL */
  T8_FREE (problem);
  *pproblem = NULL;
}
static void
t8dg_1D_advect_problem_init_elements (t8dg_1D_advect_problem_t * problem)
{
  t8_locidx_t			itree, ielement, idata;
  t8_locidx_t			num_trees, num_elems_in_tree;
  t8_element_t			*element;
//  int				iface, ineigh, isubface;
  t8dg_1D_advect_element_precomputed_values_t	*element_values;
//  t8dg_1D_advect_advance_element_data_t		*element_advance_data;
//what about advance_data?


  t8_eclass_scheme_c 		*scheme;//, *neigh_scheme;
  double			*tree_vertices;
  double			min_delta_t,delta_t;
  double			speed;

  speed = fabs(problem->flow_velocity);
  min_delta_t = problem->T - problem->t;

  num_trees = t8_forest_get_num_local_trees (problem->forest);

  for (itree = 0, idata = 0; itree < num_trees; itree++) {
    scheme =
      t8_forest_get_eclass_scheme (problem->forest,
                                   t8_forest_get_tree_class (problem->forest,
                                                             itree));
    num_elems_in_tree =
      t8_forest_get_tree_num_elements (problem->forest, itree);

    tree_vertices = t8_forest_get_tree_vertices (problem->forest, itree);

    for (ielement = 0; ielement < num_elems_in_tree; ielement++, idata++) {
      element =
        t8_forest_get_element_in_tree (problem->forest, itree, ielement);
      element_values = (t8dg_1D_advect_element_precomputed_values_t *)
            t8_sc_array_index_locidx (problem->element_values, idata);

      element_values->diameter =
        t8_forest_element_diam (problem->forest, itree, element,
                                tree_vertices);
      element_values->level = scheme->t8_element_level(element);
      element_values->num_faces = scheme->t8_element_num_faces(element);
      element_values->scaling_factor = pow(2,-element_values->level);
      t8_forest_element_coordinate(problem->forest,itree,element,tree_vertices,0,element_values->translation_vector);
      if (speed > 0) {
        delta_t = problem->cfl * element_values->diameter / speed;
      }
      min_delta_t = SC_MIN (delta_t, min_delta_t);

      t8dg_element_set_dofs_initial(problem,idata);
      t8dg_element_set_jacobian(problem,idata);
      t8dg_element_set_trafo_weights(problem,idata);

#if 0
      sc_array_t			*jacobian_at_quad;	/* d*d*Q */
       sc_array_t			*element_trafo_quad_weight;
       sc_array_t			*face_trafo_quad_weight[MAX_FACES];
       sc_array_t			*dof_values;		/*those get ghosted*/
#endif

    }
  }
  problem->delta_t = min_delta_t;
}



void t8dg_1D_advect_solve (t8_cmesh_t cmesh, t8_scalar_function_1d_fn u_0, double flow_velocity,
			   int level, int number_LGL_points, sc_MPI_Comm comm)
{
  t8dg_1D_advect_problem_t	*problem;
  problem = t8dg_1D_advect_problem_init (cmesh, u_0, flow_velocity,
  				   level, number_LGL_points, comm);


  t8_debugf("Start Advection Solve\n");
  t8dg_1D_advect_problem_init_elements (problem);


  t8dg_1D_advect_problem_destroy(&problem);
  return;
}


