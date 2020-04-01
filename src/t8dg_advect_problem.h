/*
 * t8dg_advect.h
 *
 *  Created on: Mar 22, 2020
 *      Author: lukas
 */
/** @file t8dg_advect_problem.h */

#ifndef SRC_T8DG_ADVECT_H_
#define SRC_T8DG_ADVECT_H_
#include "t8dg_numerical_flux.h"
#include "t8dg_geometry.h"

typedef struct t8dg_linear_advection_problem t8dg_linear_advection_problem_t;

t8dg_linear_advection_problem_t *t8dg_advect_problem_init_linear_geometry_1D (t8_cmesh_t cmesh,
                                                                              t8dg_scalar_function_3d_time_fn u_0,
                                                                              double flow_velocity, int level,
                                                                              int number_LGL_points, double start_time,
                                                                              double end_time, double cfl, int time_order,
                                                                              sc_MPI_Comm comm);

t8dg_linear_advection_problem_t *t8dg_advect_problem_init (t8_cmesh_t cmesh,
                                                           t8dg_coarse_geometry_3D_t * coarse_geometry,
                                                           int dim,
                                                           t8dg_scalar_function_3d_time_fn u_initial,
                                                           double flow_speed,
                                                           int level,
                                                           int number_LGL_points,
                                                           double start_time,
                                                           double end_time, double cfl, int time_order, sc_MPI_Comm comm);

void                t8dg_advect_problem_init_elements (t8dg_linear_advection_problem_t * problem);

void                t8dg_advect_problem_destroy (t8dg_linear_advection_problem_t ** pproblem);

void                t8dg_advect_runge_kutta_step (t8dg_linear_advection_problem_t * problem);

int                 t8dg_advect_problem_endtime_reached (t8dg_linear_advection_problem_t * problem);

void                t8dg_advect_problem_printdof (t8dg_linear_advection_problem_t * problem);

void                t8dg_advect_write_vtk (t8dg_linear_advection_problem_t * problem);

#if 0
t8dg_mortar_t      *t8dg_advect_element_get_face_mortar (const t8dg_linear_advection_problem_t * problem, t8_locidx_t idata, int iface);

void
 
 
 
 
 
 
 
 t8dg_advect_element_set_face_mortar (const t8dg_linear_advection_problem_t * problem, t8_locidx_t idata, int iface,
                                      t8dg_mortar_t * mortar);

t8_forest_t         t8dg_advect_problem_get_forest (t8dg_linear_advection_problem_t * problem);
#endif

#endif /* SRC_T8DG_ADVECT_H_ */
