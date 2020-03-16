/*
 * timestepping.h
 *
 *  Created on: Mar 14, 2020
 *      Author: lukas
 */

#ifndef SRC_TIMESTEPPING_H_
#define SRC_TIMESTEPPING_H_
#include <sc_containers.h>
#include "global.h"

void rungekutta_timestep(int order,const double t,const double delta_t,const t8dg_time_matrix_application f_matrix ,
			 sc_array_t *dest, sc_array_t *src, const void *application_data);


#endif /* SRC_TIMESTEPPING_H_ */
