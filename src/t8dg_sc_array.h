/*
 * t8dg_sc_array.h
 *
 *  Created on: Mar 14, 2020
 *      Author: lukas
 */
/** @file t8dg_sc_array.h */

#ifndef SRC_T8DG_SC_ARRAY_H_
#define SRC_T8DG_SC_ARRAY_H_

#include <t8.h>
#include <sc_containers.h>

/** Calculates Y = a * X + Y for sc_arrays, where each element is assumed to be a block of double
 *
 * \param [in] a            		scalar
 * \param [in] x            		allocated sc_array with elem_size % sizeof(double) == 0
 * \param [in,out] y        		allocated sc_array with elem_size % sizeof(double) == 0
 */
void                t8dg_sc_array_block_double_axpy (double a, const sc_array_t * x, sc_array_t * y);

/** Calculates Z = a * X + Y for sc_arrays, where each element is assumed to be a block of double
 *
 * \param [in] a            		scalar
 * \param [in] x           		allocated sc_array with elem_size % sizeof(double) == 0
 * \param [in] y        		allocated sc_array with elem_size % sizeof(double) == 0
 * \param [out] z        		allocated sc_array with elem_size % sizeof(double) == 0
 */
void                t8dg_sc_array_block_double_axpyz (double a, const sc_array_t * x, const sc_array_t * y, sc_array_t * z);

/** Prints the sc_array array, where each block is printed in a line
 *
 * \param [in] array        		allocated sc_array with elem_size % sizeof(double) == 0
*/
void                t8dg_sc_array_block_double_debug_print (sc_array_t * array);
void                t8dg_sc_array_block_double_print (sc_array_t * array);

/** Swaps 2 sc_arrays
 *
 * \param [in, out] parray1       	Pointer to the first sc_array
 * \param [in, out] parray2        	Pointer to the second sc_array
*/

void                t8dg_sc_array_swap (sc_array_t ** parray1, sc_array_t ** parray2);

sc_array_t         *t8dg_sc_array_clone (sc_array_t * src);

sc_array_t         *t8dg_sc_array_duplicate (sc_array_t * src);

sc_array_t         *t8dg_sc_array_new_block_double_view (sc_array_t * src, t8_locidx_t idata);

void                t8dg_sc_array_block_double_copy_subarray_into_array (const sc_array_t * src, sc_array_t * dest);
void                t8dg_sc_array_block_double_copy_array_into_subarray (const sc_array_t * src, sc_array_t * dest);
void                t8dg_sc_array_block_double_set_zero (sc_array_t * array);

#endif /* SRC_T8DG_SC_ARRAY_H_ */
