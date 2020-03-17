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


struct t8dg_advect_flux_mortar
{
  int				number_face_quadrature_points;
  t8dg_locidx_t			elem_idx_minus, elem_idx_plus;
  sc_array_t			*u_minus;
  sc_array_t 			*fluxes;
  sc_array_t			*u_plus;
  sc_array_t			*normal_vectors;/*from e- to e+ */
  t8dg_numerical_flux_1D_fn	numerical_flux;


  /*sc_array_t			*subfluxes[MAX_FACES][MAX_SUBFACES]; */
  /*subface_indices*/
};


double t8dg_upwind_flux_1D(const double u_minus,const double u_plus, const void *application_data);


#endif /* SRC_T8DG_NUMERICAL_FLUX_H_ */
