/*
 * t8dg_numerical_flux.h
 *
 *  Created on: Mar 17, 2020
 *      Author: lukas
 */

#ifndef SRC_T8DG_NUMERICAL_FLUX_H_
#define SRC_T8DG_NUMERICAL_FLUX_H_

#include <sc_containers.h>
#include "t8dg.h"

typedef struct t8dg_mortar t8dg_mortar_t;
typedef double (*t8dg_linear_numerical_flux_1D_fn)(const double u_minus, const double u_plus, const double flow_vector[3], const double normal_vector[3]);
typedef double (*t8dg_linear_flux_velocity_3D_time_fn)(const double flux_velocity[3], const double x_vec[3], double t);



struct t8dg_mortar
{
  int				number_face_quadrature_points;
  t8dg_locidx_t			elem_idx_minus, elem_idx_plus;

  /*one value for each quadrature point*/
  sc_array_t			*u_minus;
  sc_array_t			*u_plus;

  sc_array_t 			*fluxes;


  /*sc_array_t			*subfluxes[MAX_FACES][MAX_SUBFACES]; */
  /*subface_indices*/
};


double t8dg_upwind_flux_1D(const double u_minus,const double u_plus, const double flow_vector[3], const double normal_vector[3]);


#endif /* SRC_T8DG_NUMERICAL_FLUX_H_ */
