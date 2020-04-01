/*
 * timestepping.h
 *
 *  Created on: Mar 14, 2020
 *      Author: lukas
 */

/** @file t8dg_timestepping.h */
#ifndef SRC_TIMESTEPPING_H_
#define SRC_TIMESTEPPING_H_
#include <sc_containers.h>
#include "t8dg.h"

typedef void        (*t8dg_time_matrix_application) (sc_array_t * dest, const sc_array_t * src, double t, const void *application_data);

#if 0
/** implements a single step of runge kutta with a-values in the butcher-tableau only on the first minor diagonal
 * The time derivative application is given by f_matrix
 */
void                t8dg_rungekutta_timestep (int order, const double t, const double delta_t,
                                              const t8dg_time_matrix_application f_matrix, sc_array_t * dest,
                                              sc_array_t * src, const void *application_data);

#endif

void                t8dg_runge_kutta_fill_coefficients (int time_order, double **prk_a, double **prk_b, double **prk_c);

#endif /* SRC_TIMESTEPPING_H_ */
