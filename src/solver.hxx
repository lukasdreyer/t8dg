/*
 * solver.h
 *
 *  Created on: Mar 11, 2020
 *      Author: lukas
 */

#ifndef SRC_SOLVER_H_
#define SRC_SOLVER_H_

#include <t8_cmesh.h>
#include <t8_forest.h>
#include <example/common/t8_example_common.h>
#include <sc_dmatrix.h>
#include <sc_containers.h>
#include "global.h"
#include "t8dg_geometry.h"


typedef struct t8dg_1D_advect_problem t8dg_1D_advect_problem_t;
typedef struct t8dg_1D_advect_element_precomputed_values t8dg_1D_advect_element_precomputed_values_t;
typedef struct t8dg_1D_advect_advance_element_data t8dg_1D_advect_advance_element_data_t;

struct t8dg_1D_advect_problem
{
  t8_scalar_function_1d_fn 	u_0;
  double flow_velocity;

  t8_forest_t			forest;

  sc_array_t			*element_values;		/*those get partitioned*/
  sc_array_t			*jacobian_invers_at_quad;	/* d*d*Q */
  sc_array_t			*element_trafo_quad_weight;	/* Q */
  sc_array_t			*face_trafo_quad_weight[MAX_FACES]; /* FQ */

  sc_array_t			*dof_values;			/*those get ghosted*/
  sc_array_t			*dof_new;			/*those are only needed locally to save the result of a rungekutta timestep*/
  sc_array_t			*advance_element_data;		/*those need to be recalculated for each time step, remain processor local*/

  int 				level;				/*uniform refinement level*/
  int 				dim;

  double			delta_t;/*time step*/
  double 			t;/*current time*/
  double 			T;/*end time*/
  double			cfl;

  t8dg_quadrature_t			*quadrature;
  t8dg_functionbasis_t			*functionbasis;
  t8dg_coarse_geometry_t		*coarse_geometry;
  t8dg_matrix_application		vandermonde, vandermonde_transpose;
  t8dg_matrix_application		face_vandermonde, face_vandermonde_transpose;
  t8dg_numerical_flux			numerical_flux;

  t8dg_time_matrix_application		evolution_matrix;

  sc_MPI_Comm				comm;
};


/*maybe obsolete, since t8code gives all these values*/
struct t8dg_1D_advect_element_precomputed_values
{
  int level;
  int num_faces;
  double diameter;
  double scaling_factor;
  double translation_vector[3];
  /* Jacobian, quadrature/trafo weights as seperate sc_array because of size dependency on number of quadrature points*/
};



struct t8dg_1D_advect_advance_element_data
{
  sc_array_t			*dof_new;/*the new values, maybe move to single outside sc_array */
  sc_array_t 			*fluxes[MAX_FACES];
  /*sc_array_t			*subfluxes[MAX_FACES][MAX_SUBFACES]; */
  t8_locidx_t			same_level_neighbour[MAX_FACES]; /**< Indices of the neighbor elements */
};



void t8dg_1D_advect_solve (t8_cmesh_t cmesh, t8_scalar_function_1d_fn u_0, double flow_velocity,
			   int level, int number_LGL_points, sc_MPI_Comm comm);


#endif /* SRC_SOLVER_H_ */



