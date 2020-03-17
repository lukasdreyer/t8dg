/*
 * t8dg_numerical_flux.c
 *
 *  Created on: Mar 17, 2020
 *      Author: lukas
 */

#include "t8dg.h"
#include "t8dg_numerical_flux.h"

double t8dg_upwind_flux_1D(const double u_minus,const double u_plus, const void *application_data){
  T8DG_ASSERT(application_data!=NULL);
  double *c = (double*) application_data;
  if(*c > 0)return u_minus;
  return u_plus;
}
