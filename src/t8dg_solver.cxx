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

#include "t8dg_solver.hxx"
#include "t8dg_global.h"
#include "t8dg_timestepping.h"
#include "t8dg_sc_array.h"
#include "t8dg_geometry.h"

/*Access Functions for the sc_arrays that get partitioned*/
/*TODO: ASSERTS! */
static double *
t8dg_advect_element_get_dof (const t8dg_1D_advect_problem_t * problem,
                           t8_locidx_t ielement)
{
  return ((double *)
           t8_sc_array_index_locidx (problem->dof_values, ielement));
}


static t8dg_1D_advect_element_precomputed_values_t *
t8dg_advect_element_get_element_values (const t8dg_1D_advect_problem_t * problem,
                           t8_locidx_t ielement)
{
  return ((t8dg_1D_advect_element_precomputed_values_t *)
           t8_sc_array_index_locidx (problem->element_values, ielement));
}

static double *
t8dg_advect_element_get_face_quad_trafo_weights (const t8dg_1D_advect_problem_t * problem,
                           t8_locidx_t ielement, int faceindex)
{
  return ((double *)
           t8_sc_array_index_locidx (problem->face_trafo_quad_weight[faceindex], ielement));
}
static double *
t8dg_advect_element_get_element_quad_trafo_weights (const t8dg_1D_advect_problem_t * problem,
                           t8_locidx_t ielement)
{
  return ((double *)
           t8_sc_array_index_locidx (problem->element_trafo_quad_weight, ielement));
}

static double *
t8dg_advect_element_get_element_jacobian_invers_linear_array(const t8dg_1D_advect_problem_t * problem,
                           t8_locidx_t ielement)
{
  return ((double *)
           t8_sc_array_index_locidx (problem->jacobian_invers_at_quad, ielement));
}


/* fill the initial dof_values
 * for each element, iterate over quadrature points, use fine_to_coarse_geo and coarse_geo to find the image vertex;
 * Use problem->u_0 to fill dof array
 */
static void t8dg_element_set_dofs_initial(t8dg_1D_advect_problem *problem,t8_locidx_t ielement, double *tree_vertices){
  int idof;
  double *element_dof_values;
  double reference_vertex[MAX_DIM];
  double coarse_vertex[MAX_DIM];
  double image_vertex[MAX_DIM];

  t8dg_1D_advect_element_precomputed_values_t *element_values = t8dg_advect_element_get_element_values(problem,ielement);

  element_dof_values = t8dg_advect_element_get_dof(problem,ielement);


  for(idof = 0; idof < problem->functionbasis->number_of_dof; idof++){
//    get_functionbasis_vertex(reference_vertex,ielement,idof);
    //get_basisfunction_nodal vertex
    t8dg_refined_to_coarse_geometry(coarse_vertex,reference_vertex,element_values);
    //fine_to_coarse
    problem->coarse_geometry->geometry(image_vertex,coarse_vertex,tree_vertices);
    //coarse

    //u_0(coarse_vertex)

    //TEST!!
    element_dof_values[idof] = ielement * ielement + idof;
  }
}
static void t8dg_flatten_jacobian_matrix(double *flat_array,t8dg_jacobian_matrix_t jacobian_matrix, int dim){
  int ixdim,iydim;
  for(ixdim=0; ixdim < dim; ixdim++){
      for(iydim = 0; iydim < dim ; iydim++){
	  flat_array[ixdim * dim + iydim] = jacobian_matrix[ixdim][iydim];
      }
  }
}

static void t8dg_element_set_jacobian_invers_and_quad_trafo_weights(t8dg_1D_advect_problem *problem,t8_locidx_t ielement,
								    double *tree_vertices)
{
  T8_ASSERT(problem->dim == 1);
  double vertex[MAX_DIM];
  double coarse_vertex[MAX_DIM];
  double det;
  t8dg_jacobian_matrix_t coarse_jacobian_matrix,coarse_jacobian_invers,fine_jacobian_invers;
  int iquad,iface;

  t8dg_1D_advect_element_precomputed_values_t *element_values = t8dg_advect_element_get_element_values(problem,ielement);
  double *jacobian_invers_linear_array = t8dg_advect_element_get_element_jacobian_invers_linear_array(problem,ielement);
  double *element_quad_trafo = t8dg_advect_element_get_element_quad_trafo_weights(problem,ielement);
  double *face_quad_trafo[MAX_FACES];
  for(iface = 0; iface < element_values->num_faces; iface++){
   face_quad_trafo[iface] = t8dg_advect_element_get_face_quad_trafo_weights(problem,ielement,iface);
  }

  for(iquad = 0; iquad < problem->quadrature->number_of_vertices; iquad++){
      /* for each quad point,calculate its position in the coarse element, calculate the coarse jacobian, take its determinant and invert it.
       * save the inverse jacobian in a linear array, use the determinant to calculate the element_quad_trafo_weights
       * In 1D, the facetrafo weights are 1, in 2D and 3D use tangential vektor and gram determinant
       * Additionally in 2/3D we would need normalvectors.
       */

//TODO:      get_quadrature_vertex(vertex,problem,iquad);

      t8dg_refined_to_coarse_geometry(coarse_vertex, vertex,element_values);

      problem->coarse_geometry->jacobian(coarse_jacobian_matrix,coarse_vertex,tree_vertices);
      t8dg_invert_jacobian_matrix(coarse_jacobian_invers, coarse_jacobian_matrix, problem->dim);
      //apply rotation-reflection-matrix R to jacobian invers
//TODO:      apply_rotation_reflection_matrix_to_matrix(fine_jacobian_invers,coarse_jacobian_invers,element_values->idx_rotation_reflection);
      //scale by h
//TODO:      scale_jacobian_matrix(fine_jacobian_invers,element_values->scaling_factor);
      //flatten and save in jacobian_invers_linear array
      t8dg_flatten_jacobian_matrix(jacobian_invers_linear_array, fine_jacobian_invers,problem->dim);

      t8dg_determinant_jacobian_matrix(&det, coarse_jacobian_matrix, problem->dim);
      //multiply det with quadweight and save
/*TODO: get_quadrature_weight */
      element_quad_trafo[iquad] = det * ((double *) problem->quadrature->weights->array)[iquad];
  }
  for(iface = 0; iface < element_values->num_faces; iface ++){
      if(problem->dim == 1){
	  /*only one facequadrature point, with weight 1*/
	  face_quad_trafo[iface][0] = 1;
      }
  }
}

static void t8dg_1D_advect_evolution(sc_array_t *dudt_array, const sc_array_t *u_array, double t, const void *application_data){
  /** TODO: Only for testing purposes!!*/
  T8_ASSERT(dudt_array->elem_count==u_array->elem_count);
  T8_ASSERT(dudt_array->elem_size==u_array->elem_size);
  unsigned i;
  for(i = 0; i < dudt_array->elem_count * dudt_array->elem_size/sizeof(double); i++){
      ((double *)dudt_array->array)[i] =  (2. / t) * ((double *)u_array->array)[i];
  }
  /*In reality du/dt = invMassmatrix(cAu - Bu + Mg)*/

}

static t8dg_1D_advect_problem_t *
t8dg_1D_advect_problem_init (t8_cmesh_t cmesh, t8dg_scalar_function_MAX_DIMd_fn u_0, double flow_velocity,
				   int level, int number_LGL_points, sc_MPI_Comm comm)
{
  t8dg_1D_advect_problem_t 		*problem;
  t8_scheme_cxx_t			*default_scheme;
  int					iface;


  /* allocate problem */
  problem = T8_ALLOC (t8dg_1D_advect_problem_t, 1);


  problem->dim = 1;
  problem->flow_velocity = flow_velocity;
  problem->u_0 = u_0;
  problem->comm = comm;

  /*change to input!*/
  problem->T=2;
  problem->t=1;
  problem->cfl=0.1;


  problem->quadrature = t8dg_1D_LGL_quadrature(number_LGL_points);/*allocates*/
  problem->functionbasis = t8dg_1D_LGL_functionbasis(number_LGL_points);
  problem->coarse_geometry = t8dg_1D_linear_geometry();

  problem->numerical_flux = t8dg_upwind_flux_1D;

  problem->evolution_matrix = t8dg_1D_advect_evolution;

  /*If basisfunctions and quadrature use different vertices, these matrices need to be precomputed once!*/
  problem->vandermonde = t8dg_identity_matrix;
  problem->vandermonde_transpose = t8dg_identity_matrix;
  problem->face_vandermonde = t8dg_face_vandermonde_1D_linear_LGL;
  problem->face_vandermonde_transpose = t8dg_face_vandermonde_transpose_1D_linear_LGL;

  default_scheme = t8_scheme_new_default_cxx ();
  problem->forest =
    t8_forest_new_uniform (cmesh, default_scheme, level, 1, comm);

  int num_elements = t8_forest_get_num_element (problem->forest);


  problem->element_values =
    sc_array_new_count (sizeof (t8dg_1D_advect_element_precomputed_values_t),
			num_elements);

  problem->advance_element_data =
    sc_array_new_count (sizeof (t8dg_1D_advect_advance_element_data_t),
			num_elements);

  problem->jacobian_invers_at_quad =
    sc_array_new_count (sizeof (double) * problem->dim * problem->dim * problem->quadrature->number_of_vertices,
			num_elements);

  problem->element_trafo_quad_weight =
    sc_array_new_count (sizeof (double) * problem->quadrature->number_of_vertices,
			num_elements);

  for(iface=0; iface < problem->quadrature->number_of_faces; iface++){
      problem->face_trafo_quad_weight[iface] =
        sc_array_new_count (sizeof (double) * problem->quadrature->number_of_facevertices[iface],
			    num_elements);
  }/*rest auf NULL setzen ?*/


  /*currently no ghost, since serial, but generally the dof_values need to be ghosted.*/
  problem->dof_values =
    sc_array_new_count (sizeof (double) * problem->functionbasis->number_of_dof,
			num_elements +
                        t8_forest_get_num_ghosts (problem->forest));


  problem->dof_new =
    sc_array_new_count (sizeof (double) * problem->functionbasis->number_of_dof,
			num_elements);


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
  /* destroy advance element data?  */
//  t8dg_advect_problem_elements_destroy (problem);

  problem->dim = -1;
  problem->flow_velocity = 0;
  problem->u_0 = NULL;



  /* Free the arrays */
  sc_array_destroy (problem->element_values);
  sc_array_destroy (problem->dof_values);
  sc_array_destroy (problem->dof_new);
  sc_array_destroy (problem->advance_element_data);
  sc_array_destroy (problem->element_trafo_quad_weight);
  sc_array_destroy (problem->jacobian_invers_at_quad);
  for(iface = 0; iface <problem->quadrature->number_of_faces; iface++){
      sc_array_destroy (problem->face_trafo_quad_weight[iface]);
  }

  t8dg_quadrature_destroy(&(problem->quadrature));
  t8dg_functionbasis_destroy(&(problem->functionbasis));
  t8dg_coarse_geometry_destroy(&(problem->coarse_geometry));


  /* Unref the forest */
  t8_forest_unref (&problem->forest);/*unrefs coarse mesh as well*/
  /* Free the problem and set pointer to NULL */
  T8_FREE (problem);
  *pproblem = NULL;
}
static void
t8dg_1D_advect_problem_init_elements (t8dg_1D_advect_problem_t * problem)
{
  t8_locidx_t			itree, ielement, idata, idim;
  t8_locidx_t			num_trees, num_elems_in_tree;
  t8_element_t			*element;
  t8dg_1D_advect_element_precomputed_values_t	*element_values;


  t8_eclass_scheme_c 		*scheme;
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
      element_values->dim = problem->dim;
      element_values->idx_rotation_reflection = 0;

      /*TODO: is vertex0 really always the translation vector?*/
      double vertex[3];
      t8_forest_element_coordinate(problem->forest,itree,element,tree_vertices,0,vertex);
      for(idim = 0 ; idim < MAX_DIM ; idim++){
	  element_values->translation_vector[idim] = vertex[idim];
      }
      if (speed > 0) {
        delta_t = problem->cfl * element_values->diameter / speed;
      }
      min_delta_t = SC_MIN (delta_t, min_delta_t);

      /*precompute values for element idata*/
      t8dg_element_set_dofs_initial(problem,idata,tree_vertices);
      t8dg_element_set_jacobian_invers_and_quad_trafo_weights(problem,idata,tree_vertices);
    }
  }
  problem->delta_t = min_delta_t;
}



void t8dg_1D_advect_solve (t8_cmesh_t cmesh, t8dg_scalar_function_MAX_DIMd_fn u_0, double flow_velocity,
			   int level, int number_LGL_points, sc_MPI_Comm comm)
{
  int time_order=2; /*TODO: change to input*/
  t8dg_1D_advect_problem_t	*problem;

  t8_debugf("Start Advection Solve\n");

  problem = t8dg_1D_advect_problem_init (cmesh, u_0, flow_velocity,
  				   level, number_LGL_points, comm);

  t8dg_1D_advect_problem_init_elements (problem);

  t8dg_sc_array_block_double_print(problem->dof_values);


  /*Timeloop with Rungekutta timestepping: */
  while(problem->t < problem->T){
      printf("time %f\n",problem->t);
      if(problem->t + problem->delta_t > problem->T){
	  problem->delta_t = problem->T - problem->t;
      }
      t8dg_rungekutta_timestep(time_order,problem->t,problem->delta_t,problem->evolution_matrix,problem->dof_new,problem->dof_values,NULL);
      problem->t += problem->delta_t;
      /*TODO: swap problem.dof_values and problem.dof_new*/
      t8dg_sc_array_swap(&problem->dof_values,&problem->dof_new);
  }

  t8dg_sc_array_block_double_print(problem->dof_values);


  t8dg_1D_advect_problem_destroy(&problem);
  return;
}


