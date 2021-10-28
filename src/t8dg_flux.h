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

#include <t8dg.h>
#include <t8dg_flux_data_base.h>

T8DG_EXTERN_C_BEGIN ();

typedef double      (*t8dg_numerical_linear_flux3D_fn) (const double u_minus, const double u_plus,
                                                        const double flow_vector[3], const double normal_vector[3],
                                                        const void *numerical_flux_data);

typedef double      (*t8dg_numerical_linear_flux1D_fn) (const double u_minus, const double u_plus,
                                                        const double flow_constant, const double normal_component,
                                                        const void *numerical_flux_data);

typedef double      (*t8dg_numerical_flux1D_fn) (const double u_minus, const double u_plus,
                                                 const double normal_component, const void *numerical_flux_data, t8_locidx_t itree, t8_locidx_t ielement);

typedef void        (*t8dg_linear_flux1D_fn) (const double x_vec[3], double *flux_velocity, const double t, const void *flux_data, t8_locidx_t itree, t8_locidx_t ielement);

typedef void        (*t8dg_linear_flux3D_fn) (double x_vec[3], double flux_velocity[3], double t, const t8dg_flux_data_base *flux_data, t8_locidx_t itree, t8_locidx_t ielement);

T8DG_EXTERN_C_END ();

#endif /* SRC_T8DG_NUMERICAL_FLUX_H_ */
