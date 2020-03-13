/*
 * global.h
 *
 *  Created on: Mar 13, 2020
 *      Author: lukas
 */

#include <sc_containers.h>
#include "solver.hxx"

#ifndef SRC_GLOBAL_H_
#define SRC_GLOBAL_H_
void identity_matrix (sc_array_t *dest, const sc_array_t *src, const void *application_data);
void face_vandermonde_1D_linear_LGL (sc_array_t *dest, const sc_array_t *src, const void *application_data);
void face_vandermonde_transpose_1D_linear_LGL (sc_array_t *dest, const sc_array_t *src, const void *application_data);
void directional_derivative_1D_LGL2_matrix(sc_array_t *dest, const sc_array_t *src, const void *application_data);
double upwind_flux_1D(const double u_minus,const double u_plus, const void *application_data);

t8dg_quadrature_t * t8dg_1D_LGL_quadrature(int number_of_LGL);
void t8_dg_quadrature_destroy(t8dg_quadrature_t **pquadrature);


t8dg_functionbasis_t * t8dg_1D_LGL_functionbasis(int number_of_LGL);
void t8dg_functionbasis_destroy(t8dg_functionbasis_t **pfunctionbasis);



#endif /* SRC_GLOBAL_H_ */
