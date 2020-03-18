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
#include "t8dg.h"
#include "t8dg_timestepping.h"
#include "t8dg_sc_array.h"
#include "t8dg_geometry.h"
#include "t8dg_LGL.h"
#include "t8dg_numerical_flux.h"

typedef struct t8dg_advect_problem_linear_1D
{
  int 				dim = 1;

/*TODO: add source term functions*/
  t8dg_scalar_function_3d_fn 	u_0;
  double 			flow_velocity;


  t8_forest_t			forest;

  sc_array_t			*element_fine_to_coarse_geometry_data;		/*those get partitioned*/
  sc_array_t			*element_jacobian_invers_at_quad;	/* d*d*Q */
  sc_array_t			*element_trafo_quad_weight;	/* Q */
  sc_array_t			*face_trafo_quad_weight[MAX_FACES]; /* FQ */

  sc_array_t			*element_dof_values;			/*those get ghosted*/
  sc_array_t			*element_dof_values_new;			/*those are only needed locally to save the result of a rungekutta timestep*/
  sc_array_t			*face_mortar[MAX_FACES];			/*those need to be recalculated for each time step, remain processor local*/

  int 				uniform_refinement_level;				/*uniform refinement level*/

  int				time_order;

  double			delta_t;/*time step*/
  double 			t;/*current time*/
  double 			T;/*end time*/
  double			cfl;

  t8dg_LGL_quadrature_t			*quadrature;
  t8dg_LGL_functionbasis_t		*functionbasis;
  t8dg_coarse_geometry_3D_t		*coarse_geometry;
  t8dg_numerical_flux_1D_fn		numerical_flux;

  t8dg_time_matrix_application		evolution_matrix;

  sc_MPI_Comm				comm;
}t8dg_advect_problem_linear_1D_t;


/*Access Functions for the sc_arrays that get partitioned*/
static double *
t8dg_advect_element_get_element_dof_values (const t8dg_advect_problem_linear_1D_t * problem,
                           t8_locidx_t ielement)
{
  return ((double *)
           t8_sc_array_index_locidx (problem->element_dof_values, ielement));
}

static double *
t8dg_advect_element_get_element_dof_values_new (const t8dg_advect_problem_linear_1D_t * problem,
                           t8_locidx_t ielement)
{
  return ((double *)
           t8_sc_array_index_locidx (problem->element_dof_values_new, ielement));
}


static double *
t8dg_advect_element_get_face_quad_trafo_weights (const t8dg_advect_problem_linear_1D_t * problem,
                           t8_locidx_t ielement, int faceindex)
{
  return ((double *)
           t8_sc_array_index_locidx (problem->face_trafo_quad_weight[faceindex], ielement));
}
static double *
t8dg_advect_element_get_element_quad_trafo_weights (const t8dg_advect_problem_linear_1D_t * problem,
                           t8_locidx_t ielement)
{
  return ((double *)
           t8_sc_array_index_locidx (problem->element_trafo_quad_weight, ielement));
}

static double *
t8dg_advect_element_get_element_jacobian_invers_linear_array(const t8dg_advect_problem_linear_1D_t * problem,
                           t8_locidx_t ielement)
{
  return ((double *)
           t8_sc_array_index_locidx (problem->element_jacobian_invers_at_quad, ielement));
}

static t8dg_element_fine_to_coarse_geometry_data *
t8dg_advect_element_get_fine_to_coarse_geometry_data(const t8dg_advect_problem_linear_1D_t * problem,
				                           t8_locidx_t ielement){
  return ((t8dg_element_fine_to_coarse_geometry_data *)
      t8_sc_array_index_locidx(problem->element_fine_to_coarse_geometry_data, ielement));
}


/* fill the initial dof_values
 * for each element, iterate over quadrature points, use fine_to_coarse_geo and coarse_geo to find the image vertex;
 * Use problem->u_0 to fill dof array
 */
static void t8dg_element_set_dofs_initial(t8dg_advect_problem_linear_1D_t *problem,t8_locidx_t ielement, double *tree_vertices){
  int idof;
  double *element_dof_values;
  double reference_vertex[DIM3];
  double coarse_vertex[DIM3];
  double image_vertex[DIM3];

//  t8dg_element_fine_to_coarse_geometry_data_t *element_data;

//  element_data = t8dg_advect_element_get_fine_to_coarse_geometry_data(problem,ielement);
  element_dof_values = t8dg_advect_element_get_element_dof_values(problem,ielement);


  for(idof = 0; idof < problem->functionbasis->number_of_dof; idof++){
//    get_functionbasis_vertex(reference_vertex,ielement,idof);
    //get_basisfunction_nodal vertex
//    t8dg_fine_to_coarse_geometry(coarse_vertex,reference_vertex,element_data);
    //fine_to_coarse
//    problem->coarse_geometry->geometry(image_vertex,coarse_vertex,tree_vertices);
    //coarse
    element_dof_values[idof] = ielement * ielement + idof; //problem->u_0(image_vertex);
  }
}
static void t8dg_flatten_jacobian_matrix(double *flat_array,t8dg_square_3D_matrix_t jacobian_matrix, int dim){
  int ixdim,iydim;
  for(ixdim=0; ixdim < dim; ixdim++){
    for(iydim = 0; iydim < dim ; iydim++){
      flat_array[ixdim * dim + iydim] = jacobian_matrix[ixdim][iydim];
    }
  }
}

static void t8dg_element_set_jacobian_invers_and_quad_trafo_weights(t8dg_advect_problem_linear_1D_t *problem,t8_locidx_t ielement,
								    double *tree_vertices)
{
  return;



  T8_ASSERT(problem->dim == 1);
  double vertex[DIM3];
  double coarse_vertex[DIM3];
  double det;
  t8dg_square_3D_matrix_t coarse_jacobian_matrix,coarse_jacobian_invers,fine_jacobian_invers;
  int iquad,iface;

  t8dg_element_fine_to_coarse_geometry_data_t *element_data;

//  element_data = t8dg_advect_element_get_fine_to_coarse_geometry_data(problem,ielement);



  double *jacobian_invers_linear_array = t8dg_advect_element_get_element_jacobian_invers_linear_array(problem,ielement);
  double *element_quad_trafo = t8dg_advect_element_get_element_quad_trafo_weights(problem,ielement);
  double *face_quad_trafo[MAX_FACES];
  for(iface = 0; iface < problem->quadrature->vertices->number_of_faces; iface++){
   face_quad_trafo[iface] = t8dg_advect_element_get_face_quad_trafo_weights(problem,ielement,iface);
  }

  for(iquad = 0; iquad < problem->quadrature->vertices->number_of_vertices; iquad++){
      /* for each quad point,calculate its position in the coarse element, calculate the coarse jacobian, take its determinant and invert it.
       * save the inverse jacobian in a linear array, use the determinant to calculate the element_quad_trafo_weights
       * In 1D, the facetrafo weights are 1, in 2D and 3D use tangential vektor and gram determinant
       * Additionally in 2/3D we would need normalvectors.
       */

//TODO:      get_quadrature_vertex(vertex,problem,iquad);

      t8dg_fine_to_coarse_geometry(coarse_vertex, vertex,element_data);

      problem->coarse_geometry->jacobian(coarse_jacobian_matrix,coarse_vertex,tree_vertices);
//      t8dg_invert_jacobian_matrix(coarse_jacobian_invers, coarse_jacobian_matrix, problem->dim);
      //apply rotation-reflection-matrix R to jacobian invers
//TODO:      apply_rotation_reflection_matrix_to_matrix(fine_jacobian_invers,coarse_jacobian_invers,element_values->idx_rotation_reflection);
      //scale by h
//TODO:      scale_jacobian_matrix(fine_jacobian_invers,element_values->scaling_factor);
      //flatten and save in jacobian_invers_linear array
      t8dg_flatten_jacobian_matrix(jacobian_invers_linear_array, fine_jacobian_invers,problem->dim);

//      t8dg_determinant_jacobian_matrix(&det, coarse_jacobian_matrix, problem->dim);
      //multiply det with quadweight and save
/*TODO: get_quadrature_weight */
      element_quad_trafo[iquad] = det * ((double *) problem->quadrature->weights->array)[iquad];
  }
  if(problem->dim == 1){
      /*only one facequadrature point, with weight 1*/
    for(iface = 0; iface < problem->quadrature->vertices->number_of_faces; iface ++){
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

static t8dg_advect_problem_linear_1D_t *
t8dg_1D_advect_problem_init (t8_cmesh_t cmesh, t8dg_scalar_function_3d_fn u_0, double flow_velocity,
				   int level, int number_LGL_points,
				   double start_time, double end_time, double cfl, int time_order,
				   sc_MPI_Comm comm)
{
  t8dg_advect_problem_linear_1D_t 	*problem;
  t8_scheme_cxx_t			*default_scheme;
  int					iface;


  /* allocate problem */
  problem = T8_ALLOC (t8dg_advect_problem_linear_1D_t, 1);


  problem->dim = 1;
  problem->uniform_refinement_level = level;
  problem->flow_velocity = flow_velocity;
  problem->u_0 = u_0;
  problem->comm = comm;

  problem->time_order = time_order;
  problem->T=end_time;
  problem->t=start_time;
  problem->cfl=cfl;
  //problem->delta_t is dependent on elements

  t8_debugf("start LGL construction\n");
  /* these allocate memory: */
  t8dg_LGL_quadrature_and_functionbasis_new_1D(&problem->quadrature,&problem->functionbasis,number_LGL_points);
  problem->coarse_geometry = t8dg_coarse_geometry_new_1D_linear();

  problem->numerical_flux = t8dg_upwind_flux_1D;
  problem->evolution_matrix = t8dg_1D_advect_evolution;

  default_scheme = t8_scheme_new_default_cxx ();
  t8_debugf("create uniform forest\n");

  problem->forest =
    t8_forest_new_uniform (cmesh, default_scheme, level, 1, comm);

  int num_elements = t8_forest_get_num_element (problem->forest);

  t8_debugf("start creating sc_arrays\n");

  problem->element_fine_to_coarse_geometry_data =
    sc_array_new_count (sizeof (t8dg_element_fine_to_coarse_geometry_data_t),
			num_elements);

  problem->element_jacobian_invers_at_quad =
    sc_array_new_count (sizeof (double) * problem->dim * problem->dim * problem->quadrature->number_of_quadrature_points,
			num_elements);

  problem->element_trafo_quad_weight =
    sc_array_new_count (sizeof (double) * problem->quadrature->number_of_quadrature_points,
			num_elements);

  for(iface=0; iface < problem->quadrature->vertices->number_of_faces; iface++){
      problem->face_trafo_quad_weight[iface] =
        sc_array_new_count (sizeof (double) * problem->quadrature->vertices->number_of_facevertices[iface],
			    num_elements);
      problem->face_mortar[iface] =
        sc_array_new_count (sizeof (double) * problem->quadrature->vertices->number_of_facevertices[iface],
			    num_elements);
  }/*rest auf NULL setzen ?*/


  /*currently no ghost, since serial, but generally the dof_values need to be ghosted.*/
  problem->element_dof_values =
    sc_array_new_count (sizeof (double) * problem->functionbasis->number_of_dof,
			num_elements +
                        t8_forest_get_num_ghosts (problem->forest));



  problem->element_dof_values_new =
    sc_array_new_count (sizeof (double) * problem->functionbasis->number_of_dof,
			num_elements);

  t8_debugf("finished creating sc_arrays\n");


  return problem;
}

static void
t8dg_1D_advect_problem_destroy (t8dg_advect_problem_linear_1D_t ** pproblem)
{
  t8dg_advect_problem_linear_1D_t *problem;
  int iface;

  T8_ASSERT (pproblem != NULL);
  problem = *pproblem;
  if (problem == NULL) {
    return;
  }

  problem->dim = -1;
  problem->flow_velocity = 0;
  problem->u_0 = NULL;



  /* Free the arrays */
  sc_array_destroy (problem->element_fine_to_coarse_geometry_data);
  sc_array_destroy (problem->element_dof_values);
  sc_array_destroy (problem->element_dof_values_new);
  sc_array_destroy (problem->element_trafo_quad_weight);
  sc_array_destroy (problem->element_jacobian_invers_at_quad);
  for(iface = 0; iface <problem->quadrature->vertices->number_of_faces; iface++){
      sc_array_destroy (problem->face_trafo_quad_weight[iface]);
      sc_array_destroy (problem->face_mortar[iface]);
  }

  t8dg_LGL_quadrature_and_functionbasis_destroy(&problem->quadrature,&problem->functionbasis);
  t8dg_coarse_geometry_destroy(&(problem->coarse_geometry));

  /* Unref the forest */
  t8_forest_unref (&problem->forest);/*unrefs coarse mesh as well*/
  /* Free the problem and set pointer to NULL */
  T8_FREE (problem);
  *pproblem = NULL;
}
static void
t8dg_advect_problem_linear_1D_init_elements (t8dg_advect_problem_linear_1D_t * problem)
{
  t8_locidx_t			itree, ielement, idata, idim;
  t8_locidx_t			num_trees, num_elems_in_tree;
  t8_element_t			*element;


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
#if 0
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
      for(idim = 0 ; idim < DIM3 ; idim++){
	  element_values->translation_vector[idim] = vertex[idim];
      }
#endif
      if (speed > 0) {
        delta_t = 0.1 ;//problem->cfl * element_values->diameter / speed;
      }
      min_delta_t = SC_MIN (delta_t, min_delta_t);

      /*precompute values for element idata*/
      t8dg_element_set_dofs_initial(problem,idata,tree_vertices);
      t8dg_element_set_jacobian_invers_and_quad_trafo_weights(problem,idata,tree_vertices);
    }
  }
  problem->delta_t = min_delta_t;
}

void t8dg_advect_evolve(t8dg_advect_problem_linear_1D_t *problem){
  t8dg_rungekutta_timestep(problem->time_order,problem->t,problem->delta_t,problem->evolution_matrix,
			   problem->element_dof_values_new,problem->element_dof_values,NULL);
  problem->t += problem->delta_t;
  /*TODO: swap problem.dof_values and problem.dof_new*/
  t8dg_sc_array_swap(&problem->element_dof_values,&problem->element_dof_values_new);

}

void t8dg_1D_advect_solve (t8_cmesh_t cmesh, t8dg_scalar_function_3d_fn u_0, double flow_velocity,
			   int level, int number_LGL_points,
			   double start_time, double end_time, double cfl, int time_order,
			   sc_MPI_Comm comm)
{
  t8dg_advect_problem_linear_1D_t	*problem;

  t8_debugf("Start Advection Solve\n");

  problem = t8dg_1D_advect_problem_init (cmesh, u_0, flow_velocity,
  				   level, number_LGL_points,
				   start_time, end_time, cfl, time_order,
				   comm);
  t8dg_advect_problem_linear_1D_init_elements (problem);

  t8dg_sc_array_block_double_print(problem->element_dof_values);


  /*Timeloop with Rungekutta timestepping: */
  while(problem->t < problem->T){
      printf("time %f\n",problem->t);
      if(problem->t + problem->delta_t > problem->T){
	  problem->delta_t = problem->T - problem->t;
      }
      t8dg_advect_evolve(problem);
  }

  t8dg_sc_array_block_double_print(problem->element_dof_values);

  t8dg_1D_advect_problem_destroy(&problem);
  return;
}


