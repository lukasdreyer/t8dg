/*
 * timestepping.h
 *
 *  Created on: Mar 14, 2020
 *      Author: lukas
 */

#ifndef SRC_TIMESTEPPING_H_
#define SRC_TIMESTEPPING_H_
#include <sc_containers.h>
#include "t8dg_global.h"


/*src nicht const wegen sc_array_copy*/
/* implements a single step of runge kutta with a-values in the butcher-tableau only on the first minor diagonal
 * The time derivative application is given by f_matrix
 */
void t8dg_rungekutta_timestep(int order,const double t,const double delta_t,const t8dg_time_matrix_application f_matrix ,
			 sc_array_t *dest, sc_array_t *src, const void *application_data);


#endif /* SRC_TIMESTEPPING_H_ */
