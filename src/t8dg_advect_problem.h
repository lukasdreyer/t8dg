/*
 * t8dg_advect.h
 *
 *  Created on: Mar 22, 2020
 *      Author: lukas
 */
/** @file t8dg_advect_problem.h */

#ifndef SRC_T8DG_ADVECT_H_
#define SRC_T8DG_ADVECT_H_

#include <t8_cmesh.h>
#include <sc.h>

#include "t8dg_flux.h"
#include "t8dg_timestepping.h"
#include "t8dg_coarse_geometry.h"

T8DG_EXTERN_C_BEGIN ();

typedef struct t8dg_linear_advection_problem t8dg_linear_advection_problem_t;

typedef struct t8dg_linear_advection_problem_description
{
  t8dg_scalar_function_3d_time_fn initial_condition_fn;             /**< Initial condition function */

  t8dg_linear_flux3D_fn velocity_field;
  void               *flux_data;

  t8dg_scalar_function_3d_time_fn source_sink_fn;
  t8dg_scalar_function_3d_time_fn analytical_sol_fn;             /**< Analytical solution function */
} t8dg_linear_advection_problem_description_t;

t8dg_linear_advection_problem_t *t8dg_advect_problem_init (t8_cmesh_t cmesh,
                                                           t8dg_coarse_geometry_t * coarse_geometry,
                                                           int dim,
                                                           t8dg_linear_advection_problem_description_t * description,
                                                           int uniform_level, int max_level,
                                                           int number_LGL_points, t8dg_timestepping_data_t * time_data,
                                                           t8_forest_adapt_t adapt_fn, sc_MPI_Comm comm);

void                t8dg_advect_problem_destroy (t8dg_linear_advection_problem_t ** pproblem);

/*do the important stuff*/
void                t8dg_advect_problem_advance_timestep (t8dg_linear_advection_problem_t * problem);

void                t8dg_advect_problem_adapt (t8dg_linear_advection_problem_t * problem, int measure_time);

void                t8dg_advect_problem_partition (t8dg_linear_advection_problem_t * problem, int measure_time);

/*stats*/
void                t8dg_advect_problem_compute_and_print_stats (t8dg_linear_advection_problem_t * problem);

/*error*/
double              t8dg_advect_problem_l_infty_rel (t8dg_linear_advection_problem_t * problem);

double              t8dg_advect_problem_l2_rel (t8dg_linear_advection_problem_t * problem);

/*output*/
void                t8dg_advect_problem_printdof (t8dg_linear_advection_problem_t * problem);

void                t8dg_advect_problem_write_vtk (t8dg_linear_advection_problem_t * problem);

/*getter*/
int                 t8dg_advect_problem_endtime_reached (t8dg_linear_advection_problem_t * problem);

int                 t8dg_advect_problem_get_stepnumber (t8dg_linear_advection_problem_t * problem);

T8DG_EXTERN_C_END ();

#endif /* SRC_T8DG_ADVECT_H_ */
