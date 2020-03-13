/*
 * solver.h
 *
 *  Created on: Mar 11, 2020
 *      Author: lukas
 */

#ifndef SRC_SOLVER_H_
#define SRC_SOLVER_H_

#include <t8_cmesh.h>
#include <example/common/t8_example_common.h>

#define MAX_FACES 2

typedef struct
{
  t8_scalar_function_1d_fn 	u_0;
  double flow_velocity;

  t8_forest_t			forest;
  sc_array_t			*element_values;	/*those get partitioned*/
  sc_array_t			*dof_values;		/*those get ghosted*/
  sc_array_t			*advance_element_data;/*those need to be recalculated for each time step, remain processor local*/

  int 				level;
  int 				dim;
  int				number_LGL_points;
}t8dg_1D_advect_problem_t;

typedef struct
{
  int level;
  int num_faces;
  /* Jacobian, quadrature/trafo weights*/
}t8dg_1D_advect_element_precomputed_values_t;

typedef struct
{
  sc_array_t			u_new;
  double 			*fluxes[MAX_FACES];
  /*double			*subfluxes[MAX_FACES][MAX_SUBFACES]; */

  t8_locidx_t			*neighs[MAX_FACES]; /**< Indices of the neighbor elements */
}t8dg_1D_advect_advance_element_data_t;


t8dg_1D_advect_problem_t *
t8dg_1D_advect_problem_init (t8_cmesh_t cmesh, t8_scalar_function_1d_fn u_0, double flow_velocity,
				   int level, int number_LGL_points, sc_MPI_Comm comm);

void t8dg_1D_advect_solve (t8_cmesh_t cmesh, t8_scalar_function_1d_fn u_0, double flow_velocity,
			   int level, int number_LGL_points, sc_MPI_Comm comm);


#endif /* SRC_SOLVER_H_ */



