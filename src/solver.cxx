/*
 * solver.c
 *
 *  Created on: Mar 11, 2020
 *      Author: lukas
 */
#include <t8.h>
#include <example/common/t8_example_common.h>
#include <t8_schemes/t8_default_cxx.hxx>
#include "solver.hxx"

static void
t8dg_1D_advect_problem_destroy (t8dg_1D_advect_problem_t ** pproblem)
{
  t8dg_1D_advect_problem_t *problem;

  T8_ASSERT (pproblem != NULL);
  problem = *pproblem;
  if (problem == NULL) {
    return;
  }
  /* destroy elements */
//  t8_advect_problem_elements_destroy (problem);

  problem->dim = -1;
  problem->flow_velocity = 0;
  problem->number_LGL_points = -1;
  problem->u_0 = NULL;
  /* Free the element array */
  sc_array_destroy (problem->element_values);
  sc_array_destroy (problem->dof_values);
  sc_array_destroy (problem->advance_element_data);
  /* Unref the forest */
  t8_forest_unref (&problem->forest);//coarse mesh?
  /* Free the problem and set pointer to NULL */
  T8_FREE (problem);
  *pproblem = NULL;
}

void t8dg_1D_advect_solve (t8_cmesh_t cmesh, t8_scalar_function_1d_fn u_0, double flow_velocity,
			   int level, int number_LGL_points, sc_MPI_Comm comm)
{
  t8dg_1D_advect_problem_t	*problem;
  problem = t8dg_1D_advect_problem_init (cmesh, u_0, flow_velocity,
  				   level, number_LGL_points, comm);


  t8_debugf("Start Advection Solve\n");


  t8dg_1D_advect_problem_destroy(&problem);
  return;
}


t8dg_1D_advect_problem_t *
t8dg_1D_advect_problem_init (t8_cmesh_t cmesh, t8_scalar_function_1d_fn u_0, double flow_velocity,
				   int level, int number_LGL_points, sc_MPI_Comm comm)
{
  t8dg_1D_advect_problem_t 		*problem;
  t8_scheme_cxx_t			*default_scheme;



  /* allocate problem */
  problem = T8_ALLOC (t8dg_1D_advect_problem_t, 1);

  problem->dim = 1;
  problem->flow_velocity = flow_velocity;
  problem->number_LGL_points = number_LGL_points;
  problem->u_0 = u_0;

  default_scheme = t8_scheme_new_default_cxx ();
  problem->forest =
    t8_forest_new_uniform (cmesh, default_scheme, level, 1, comm);

  problem->element_values =
    sc_array_new_count (sizeof (t8dg_1D_advect_element_precomputed_values_t),
                        t8_forest_get_num_element (problem->forest));

  problem->advance_element_data =
    sc_array_new_count (sizeof (t8dg_1D_advect_advance_element_data_t),
                        t8_forest_get_num_element (problem->forest));

  problem->dof_values =
    sc_array_new_count (sizeof (double) * problem->number_LGL_points,
                        t8_forest_get_num_element (problem->forest) +
                        t8_forest_get_num_ghosts (problem->forest));


  return problem;
}
