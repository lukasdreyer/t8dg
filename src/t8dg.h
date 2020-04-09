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
#define T8DG_ALLOC_ZERO T8_ALLOC_ZERO
#define T8DG_FREE T8_FREE

#define T8DG_CHECK_ABORT SC_CHECK_ABORT
#define T8DG_ABORT SC_ABORT

#define DIM3 3
#define MAX_FACES 2
#define MAX_SUBFACES 1
#define MAX_SUBELEMENTS 2

typedef t8_locidx_t t8dg_locidx_t;

/**A timedependent scalar function f:R^3 x R^+ -> R*/
typedef double      (*t8dg_scalar_function_3d_time_fn) (const double x[DIM3], const double t);
typedef double      (*t8dg_scalar_function_3d_fn) (const double x[DIM3]);

void               *t8dg_sc_array_index_locidx (const sc_array_t * array, t8dg_locidx_t it);

t8_locidx_t         t8dg_itree_ielement_to_idata (t8_forest_t forest, t8_locidx_t itree, t8_locidx_t ielement);
void                t8dg_vec_print (double x[3]);

/** Query the package identity as registered in libsc.
 * \return          This is -1 before \ref t8dg_init has been called
 *                  and a proper package identifier afterwards.
 */
int                 t8dg_get_package_id (void);

/** Logging function parametrized by local/global category and priority.
 * \param [in] category     Either SC_LC_NORMAL for outputting on every rank
 *                          or SC_LC_GLOBAL for outputting on the root rank.
 * \param [in] priority     Please see sc.h for legal log priorities.
 * \param [in] fmt          Printf-style format string.
 * \param [in] ap           Argument list; see stdarg.h.
 */
void                t8dg_logv (int category, int priority, const char *fmt, va_list ap);

/* *INDENT-OFF* */
/** Logging function parametrized by local/global category and priority.
 * \param [in] category     Either SC_LC_NORMAL for outputting on every rank
 *                          or SC_LC_GLOBAL for outputting on the root rank.
 * \param [in] priority     Please see sc.h for legal log priorities.
 * \param [in] fmt          Printf-style format string.
 */
void                t8dg_logf (int category, int priority, const char *fmt, ...)
#ifndef T8_DOXYGEN
  __attribute__ ((format (printf, 3, 4)))
#endif
  ;
#if 0
/** Add one space to the start of t8's default log format. */
void                t8dg_log_indent_push (void);

/** Remove one space from the start of a t8's default log format. */
void                t8dg_log_indent_pop (void);
#endif

/** Log a message on the root rank with priority SC_LP_ERROR.
 * \param [in] fmt          Printf-style format string.
 */
void                t8dg_global_errorf (const char *fmt, ...)
#ifndef T8_DOXYGEN
  __attribute__ ((format (printf, 1, 2)))
#endif
  ;

/** Log a message on the root rank with priority SC_LP_ESSENTIAL.
 * \param [in] fmt          Printf-style format string.
 */
void                t8dg_global_essentialf (const char *fmt, ...)
#ifndef T8_DOXYGEN
  __attribute__ ((format (printf, 1, 2)))
#endif
  ;

/** Log a message on the root rank with priority SC_LP_PRODUCTION.
 * \param [in] fmt          Printf-style format string.
 */
void                t8dg_global_productionf (const char *fmt, ...)
#ifndef T8_DOXYGEN
  __attribute__ ((format (printf, 1, 2)))
#endif
  ;

/** Log a message on the root rank with priority SC_LP_INFO.
 * \param [in] fmt          Printf-style format string.
 */
void                t8dg_global_infof (const char *fmt, ...)
#ifndef T8_DOXYGEN
  __attribute__ ((format (printf, 1, 2)))
#endif
  ;

/** Log a message, no matter what rank, with priority SC_LP_INFO.
 * \param [in] fmt          Printf-style format string.
 */
void                t8dg_infof (const char *fmt, ...)
#ifndef T8_DOXYGEN
  __attribute__ ((format (printf, 1, 2)))
#endif
  ;

/** Log a message, no matter what rank, with priority SC_LP_DEBUG.
 * \param [in] fmt          Printf-style format string.
 * \note This function does not print anything unless t8code was compiled
 * in debug mode (--enable-debug, T8_ENABLE_DEBUG was defined).
 */
void                t8dg_debugf (const char *fmt, ...)
#ifndef T8_DOXYGEN
  __attribute__ ((format (printf, 1, 2)))
#endif
  ;

/** Log a message, no matter what rank, with priority SC_LP_ERROR.
 * \param [in] fmt          Printf-style format string.
 */
void                t8dg_errorf (const char *fmt, ...)
#ifndef T8_DOXYGEN
  __attribute__ ((format (printf, 1, 2)))
#endif
  ;

/* *INDENT-ON* */

/** Register t8code with libsc and print version and variable information.
 * \param [in] log_threshold Declared in sc.h.  SC_LP_DEFAULT is fine.
 *                           You can also choose from log levels SC_LP_*.
 */
void                t8dg_init (int log_threshold);

/* call this at the end of a header file to match T8_EXTERN_C_BEGIN (). */
T8DG_EXTERN_C_END ();

#endif /* SRC_T8DG_H_ */
