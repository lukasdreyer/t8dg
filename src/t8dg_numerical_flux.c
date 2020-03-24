/*
 * t8dg_numerical_flux.c
 *
 *  Created on: Mar 17, 2020
 *      Author: lukas
 */

#include "t8dg.h"
#include "t8dg_numerical_flux.h"
#include <t8_vec.h>

double t8dg_upwind_flux_1D(const double u_minus,const double u_plus, const double flow_vector[3], const double normal_vector[3]){
  /*TODO: ASSERT linear dependence of flow_vector and normal_vector*/
  if(t8_vec_dot(flow_vector,normal_vector)>0) return u_minus * t8_vec_dot(flow_vector,normal_vector);
  else return u_plus * t8_vec_dot(flow_vector,normal_vector);
}
