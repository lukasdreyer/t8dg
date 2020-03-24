/*
 * t8dg.h
 *
 *  Created on: Mar 17, 2020
 *      Author: lukas
 */

#ifndef SRC_T8DG_H_
#define SRC_T8DG_H_

#include <t8.h>

/** This macro opens the extern "C" brace needed in C headers.
 * It needs to be followed by a semicolon to look like a statement. */
#define T8DG_EXTERN_C_BEGIN() T8_EXTERN_C_BEGIN()

/** This macro closes the extern "C" brace needed in C headers.
 * It needs to be followed by a semicolon to look like a statement. */
#define T8DG_EXTERN_C_END() T8_EXTERN_C_END()

/* call this after including all headers */
T8DG_EXTERN_C_BEGIN ();

#define T8DG_ASSERT T8_ASSERT
#define T8DG_ALLOC T8_ALLOC
#define T8DG_FREE T8_FREE

#define DIM3 3
#define MAX_FACES 2

typedef t8_locidx_t t8dg_locidx_t;

typedef double      (*t8dg_scalar_function_3d_fn) (const double x[DIM3]);
typedef double      (*t8dg_scalar_function_3d_time_fn) (const double x[DIM3],const double t);

typedef void (*t8dg_matrix_application) (sc_array_t *dest, const sc_array_t *src, const void *application_data);
typedef void (*t8dg_time_matrix_application) (sc_array_t *dest, const sc_array_t *src, double t, const void *application_data);

/* call this at the end of a header file to match T8_EXTERN_C_BEGIN (). */
T8DG_EXTERN_C_END ();


#endif /* SRC_T8DG_H_ */
