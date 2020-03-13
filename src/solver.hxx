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


#define MAX_FACES 2

typedef struct t8dg_quadrature t8dg_quadrature_t;
typedef struct t8dg_functionbasis t8dg_functionbasis_t;
typedef struct t8dg_coarse_geometry t8dg_coarse_geometry_t;
typedef struct t8dg_1D_advect_problem t8dg_1D_advect_problem_t;
typedef struct t8dg_1D_advect_element_precomputed_values t8dg_1D_advect_element_precomputed_values_t;
typedef struct t8dg_1D_advect_advance_element_data t8dg_1D_advect_advance_element_data_t;

typedef void (*t8dg_matrix_application) (sc_array_t *dest, const sc_array_t *src, const void *application_data);
typedef double (*t8dg_numerical_flux)(const double u_minus,const double u_plus, const void *application_data);
typedef sc_dmatrix_t* (*jacobian_fn)(const double vertex[3]);
typedef sc_array_t* (*geometry_fn)(const double vertex[3]);/*TODO: maybe change interface? */


struct t8dg_coarse_geometry
{
  jacobian_fn		jacobian;
  geometry_fn		geometry;
  void			*geometry_data;
  t8dg_1D_advect_problem_t	*problem;

};

struct t8dg_quadrature
{
  int 			tensorflag;
  t8dg_quadrature_t 	*tensor1;
  t8dg_quadrature_t 	*tensor2;
  int			number_of_vertices;
  int			number_of_faces;
  int			number_of_facevertices[MAX_FACES];
  sc_array_t		*vertices;
  sc_array_t		*weights;
  /* how do i get to face-vertices/weights */
};

struct t8dg_functionbasis
{
  int number_of_dof;
  t8dg_matrix_application		directional_derivative_matrix;
};

struct t8dg_1D_advect_problem
{
  t8_scalar_function_1d_fn 	u_0;
  double flow_velocity;

  t8_forest_t			forest;
  sc_array_t			*element_values;	/*those get partitioned*/
  sc_array_t			*jacobian_at_quad;	/* d*d*Q */
  sc_array_t			*element_trafo_quad_weight;
  sc_array_t			*face_trafo_quad_weight[MAX_FACES];
  sc_array_t			*dof_values;		/*those get ghosted*/
  sc_array_t			*advance_element_data;/*those need to be recalculated for each time step, remain processor local*/

  int 				level;	/*uniform refinement level*/
  int 				dim;
//  int				number_LGL_points;

  double			delta_t,t,T;/*current time, end time*/
  double			cfl;

  t8dg_quadrature_t			*quadrature;
  t8dg_functionbasis_t			*functionbasis;
  t8dg_matrix_application		vandermonde, vandermonde_transpose;
  t8dg_matrix_application		face_vandermonde, face_vandermonde_transpose;
  t8dg_numerical_flux			numerical_flux;

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
#if 0
  int orientation;/*defines orthonormal matrix R*/
#endif
  /* Jacobian, quadrature/trafo weights as seperate sc_array because of variable size*/
};



struct t8dg_1D_advect_advance_element_data
{
  sc_array_t			dof_new;/* */
  sc_array_t 			*fluxes[MAX_FACES];
  /*sc_array_t			*subfluxes[MAX_FACES][MAX_SUBFACES]; */
  t8_locidx_t			*neighs[MAX_FACES]; /**< Indices of the neighbor elements */
};

/*
t8dg_1D_advect_problem_t *
t8dg_1D_advect_problem_init (t8_cmesh_t cmesh, t8_scalar_function_1d_fn u_0, double flow_velocity,
				   int level, int number_LGL_points, sc_MPI_Comm comm);
*/

void t8dg_1D_advect_solve (t8_cmesh_t cmesh, t8_scalar_function_1d_fn u_0, double flow_velocity,
			   int level, int number_LGL_points, sc_MPI_Comm comm);

t8dg_quadrature_t * t8dg_1D_LGL_quadrature(int number_of_LGL);

#endif /* SRC_SOLVER_H_ */



