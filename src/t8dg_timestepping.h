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



/** @file t8dg_timestepping.h */
#ifndef SRC_TIMESTEPPING_H_
#define SRC_TIMESTEPPING_H_
#include <sc_containers.h>
#include "t8dg.h"
#include "t8dg_dof.h"

T8DG_EXTERN_C_BEGIN ();

typedef struct t8dg_timestepping_data t8dg_timestepping_data_t;

typedef void        (*t8dg_time_matrix_application) (t8dg_dof_values_t * src_dof, t8dg_dof_values_t * dest_dof, const double t,
                                                     const void *application_data);

/** implements a single step of runge kutta with a-values in the butcher-tableau only on the first minor diagonal
 * The time derivative application is given by time_derivative
 */
void
 
 
 
 
 
 
 
 t8dg_timestepping_runge_kutta_step (t8dg_time_matrix_application time_derivative,
                                     t8dg_timestepping_data_t * time_data, t8dg_dof_values_t ** pdof_array, void *user_data);

void                t8dg_runge_kutta_fill_coefficients (int time_order, double **prk_a, double **prk_b, double **prk_c);

t8dg_timestepping_data_t *t8dg_timestepping_data_new_cfl (int time_order, double start_time, double end_time, double cfl);

t8dg_timestepping_data_t *t8dg_timestepping_data_new_constant_timestep (int time_order, double start_time, double end_time, double delta_t);

void                t8dg_timestepping_data_destroy (t8dg_timestepping_data_t ** ptime_data);

double              t8dg_timestepping_data_get_current_time (t8dg_timestepping_data_t * time_data);

void                t8dg_timestepping_data_set_current_time (t8dg_timestepping_data_t * time_data, double current_time);

double              t8dg_timestepping_data_get_end_time (t8dg_timestepping_data_t * time_data);

double              t8dg_timestepping_data_get_cfl (t8dg_timestepping_data_t * time_data);

int                 t8dg_timestepping_data_get_time_order (t8dg_timestepping_data_t * time_data);

double              t8dg_timestepping_data_get_time_step (t8dg_timestepping_data_t * time_data);

void                t8dg_timestepping_data_set_time_step (t8dg_timestepping_data_t * time_data, double delta_t);

int                 t8dg_timestepping_data_is_endtime_reached (t8dg_timestepping_data_t * time_data);

int                 t8dg_timestepping_data_get_step_number (t8dg_timestepping_data_t * time_data);

void                t8dg_timestepping_data_increase_step_number (t8dg_timestepping_data_t * time_data);

double              t8dg_timestepping_data_get_time_left (t8dg_timestepping_data_t * time_data);

T8DG_EXTERN_C_END ();

#endif /* SRC_TIMESTEPPING_H_ */
