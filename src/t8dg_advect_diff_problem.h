/*
 *  This file is part of t8dg.
 *  t8dg is a discontinuous galerkin solver on adaptive meshes
 *  that uses the adaptive meshing library t8code to distribute
 *  the workload equally among parallel processes to achieve good scaling.
 *
 *  Copyright (C) 2020 the developers
 *
 *  t8dg is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  t8dg is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with t8dg.  If not, see <http://www.gnu.org/licenses/>.
 */

/** @file t8dg_advect_diff_problem.h */

#ifndef SRC_T8DG_ADVECT_DIFF_H_
#define SRC_T8DG_ADVECT_DIFF_H_

#include <t8_cmesh.h>
#include <sc.h>

#include "t8dg_flux.h"
#include "t8dg_timestepping.h"
#include "t8dg_coarse_geometry.h"
#include "t8dg_values.h"
#include "t8dg_adapt.h"
#include "t8dg_output.h"

T8DG_EXTERN_C_BEGIN ();

typedef struct t8dg_linear_advection_diffusion_problem t8dg_linear_advection_diffusion_problem_t;

typedef struct t8dg_linear_advection_diffusion_problem_description
{
  int                 dim;
  t8dg_scalar_function_3d_time_fn initial_condition_fn;             /**< Initial condition function */
  void               *initial_condition_data;
  t8dg_scalar_function_3d_time_fn source_sink_fn;
  void               *source_sink_data;
  t8dg_scalar_function_3d_time_fn analytical_sol_fn;             /**< Analytical solution function */
  void               *analytical_sol_data;

  t8dg_linear_flux3D_fn velocity_field;
  void               *flux_data;

  double              diffusion_coefficient;
  t8dg_numerical_linear_flux3D_fn numerical_flux_advection;
  void               *numerical_flux_advection_data;
  t8dg_numerical_flux1D_fn numerical_flux_diffusion_gradient;
  void               *numerical_flux_diffusion_gradient_data;
  t8dg_numerical_flux1D_fn numerical_flux_diffusion_concentration;
  void               *numerical_flux_diffusion_concentration_data;

} t8dg_linear_advection_diffusion_problem_description_t;

t8dg_linear_advection_diffusion_problem_t *t8dg_advect_diff_problem_init_arguments (int icmesh,
                                                                                    int initial_level,
                                                                                    int number_LGL_points,
                                                                                    int initial_cond_arg,
                                                                                    double flow_speed,
                                                                                    double diffusion_coefficient,
                                                                                    double start_time,
                                                                                    double end_time,
                                                                                    double cfl,
                                                                                    double delta_t,
                                                                                    int time_steps,
                                                                                    int time_order,
                                                                                    int min_level,
                                                                                    int max_level,
                                                                                    int adapt_arg,
                                                                                    int adapt_freq,
                                                                                    const char *prefix,
                                                                                    int vtk_freq, int numerical_flux_arg,
                                                                                    int source_sink_arg, int refine_error,
                                                                                    sc_MPI_Comm comm);

t8dg_linear_advection_diffusion_problem_t *t8dg_advect_diff_problem_init (t8_forest_t forest,
                                                                          t8dg_linear_advection_diffusion_problem_description_t *
                                                                          description, t8dg_values_t * dg_values,
                                                                          t8dg_timestepping_data_t * time_data,
                                                                          t8dg_adapt_data_t * adapt_data, t8dg_vtk_data_t * vtk_data,
                                                                          double init_time, int refine_error, sc_MPI_Comm comm);

void                t8dg_advect_diff_problem_destroy (t8dg_linear_advection_diffusion_problem_t ** pproblem);

/*do the important stuff*/
void                t8dg_advect_diff_solve (t8dg_linear_advection_diffusion_problem_t * problem);

void                t8dg_advect_diff_problem_advance_timestep (t8dg_linear_advection_diffusion_problem_t * problem);

void                t8dg_advect_diff_problem_adapt (t8dg_linear_advection_diffusion_problem_t * problem, int measure_time);

void                t8dg_advect_diff_problem_adapt_uniform (t8dg_linear_advection_diffusion_problem_t * problem, int uniform_level);

void                t8dg_advect_diff_problem_partition (t8dg_linear_advection_diffusion_problem_t * problem, int measure_time);

/*stats*/
void                t8dg_advect_diff_problem_compute_and_print_stats (t8dg_linear_advection_diffusion_problem_t * problem);

/*error*/
double              t8dg_advect_diff_problem_l_infty_rel (t8dg_linear_advection_diffusion_problem_t * problem);

double              t8dg_advect_diff_problem_l2_rel (t8dg_linear_advection_diffusion_problem_t * problem);

/*output*/
void                t8dg_advect_diff_problem_printdof (t8dg_linear_advection_diffusion_problem_t * problem);

void                t8dg_advect_diff_problem_write_vtk (t8dg_linear_advection_diffusion_problem_t * problem);

/*getter*/
int                 t8dg_advect_diff_problem_endtime_reached (t8dg_linear_advection_diffusion_problem_t * problem);

int                 t8dg_advect_diff_problem_get_stepnumber (t8dg_linear_advection_diffusion_problem_t * problem);

void                t8dg_advect_diff_problem_set_time_step (t8dg_linear_advection_diffusion_problem_t * problem);

int                 t8dg_advect_diff_problem_get_apx_total_steps (t8dg_linear_advection_diffusion_problem_t * problem);

T8DG_EXTERN_C_END ();

#endif /* SRC_T8DG_ADVECT_DIFF_H_ */
