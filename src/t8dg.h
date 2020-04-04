/*
 * t8dg.h
 *
 *  Created on: Mar 17, 2020
 *      Author: lukas
 */
/** @file t8dg.h */

#ifndef SRC_T8DG_H_
#define SRC_T8DG_H_

#include <t8.h>
#include <sc.h>
#include <t8_forest.h>          /*maybe replace by typedef struct t8_forest *t8_forest_t */

/** This macro opens the extern "C" brace needed in C headers.
 * It needs to be followed by a semicolon to look like a statement. */
#define T8DG_EXTERN_C_BEGIN() T8_EXTERN_C_BEGIN()

/** This macro closes the extern "C" brace needed in C headers.
 * It needs to be followed by a semicolon to look like a statement. */
#define T8DG_EXTERN_C_END() T8_EXTERN_C_END()

/* call this after including all headers */
T8_EXTERN_C_BEGIN ();

#define T8DG_ASSERT T8_ASSERT
#define T8DG_ALLOC T8_ALLOC
#define T8DG_FREE T8_FREE

#define T8DG_CHECK_ABORT SC_CHECK_ABORT
#define T8DG_ABORT SC_ABORT

#define DIM3 3
#define MAX_FACES 2
#define MAX_SUBFACES 1
#define MAX_SUBELEMENTS 2

/**A timedependent scalar function f:R^3 x R^+ -> R*/
typedef double      (*t8dg_scalar_function_3d_time_fn) (const double x[DIM3], const double t);

t8_locidx_t         t8dg_itree_ielement_to_idata (t8_forest_t forest, t8_locidx_t itree, t8_locidx_t ielement);
void                t8dg_vec_print (double x[3]);

/* call this at the end of a header file to match T8_EXTERN_C_BEGIN (). */
T8DG_EXTERN_C_END ();

#endif /* SRC_T8DG_H_ */
