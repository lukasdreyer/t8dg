/*
 * t8dg_sc_array.h
 *
 *  Created on: Mar 14, 2020
 *      Author: lukas
 */

#ifndef SRC_T8DG_SC_ARRAY_H_
#define SRC_T8DG_SC_ARRAY_H_

#include<sc_containers.h>

/* All array need to be allocated before*/
/* y = ax+y, where both array have elem_size = multiple of sizeof(double)*/
void                t8dg_sc_array_block_double_axpy (double a, const sc_array_t * x, sc_array_t * y);

/*z = ax+y, where both array have elem_size = multiple of sizeof(double)*/
void                t8dg_sc_array_block_double_zaxpy (sc_array_t * z, double a, const sc_array_t * x, const sc_array_t * y);

void                t8dg_sc_array_block_double_print (sc_array_t * array);

void                t8dg_sc_array_swap (sc_array_t ** parray1, sc_array_t ** parray2);

#endif /* SRC_T8DG_SC_ARRAY_H_ */
