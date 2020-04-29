/*
 * t8dg_numerical_flux.h
 *
 *  Created on: Mar 17, 2020
 *      Author: lukas
 */
/** @file t8dg_numerical_flux.h */

#ifndef SRC_T8DG_NUMERICAL_FLUX_H_
#define SRC_T8DG_NUMERICAL_FLUX_H_

#include <t8.h>
#include <t8_forest.h>

#include "t8dg.h"

typedef struct t8dg_flux t8dg_flux_t;

t8dg_flux_t        *t8dg_flux_new_linear_constant_flux (const double flow_direction[3], const double flow_velocity);

void                t8dg_flux_destroy (t8dg_flux_t ** pflux);

void                t8dg_flux_calulate_flux (const t8dg_flux_t * flux, const double x_vec[3], double flux_vec[3], const double t);

double              t8dg_flux_calculate_numerical_flux_value (const t8dg_flux_t * flux, const double u_minus, const double u_plus,
                                                              const double flow_vector[3], const double normal_vector[3]);
#endif /* SRC_T8DG_NUMERICAL_FLUX_H_ */
