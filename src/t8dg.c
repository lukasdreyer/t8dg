/*
 * t8dg.c
 *
 *  Created on: Mar 11, 2020
 *      Author: lukas
 */

#include <t8_forest.h>
#include <t8.h>
#include "t8dg.h"

static int          t8dg_package_id = -1;

int
t8dg_get_package_id (void)
{
  return t8dg_package_id;
}

void
t8dg_logv (int category, int priority, const char *fmt, va_list ap)
{
  sc_logv ("unknown", -1, t8dg_package_id, category, priority, fmt, ap);
}

void
t8dg_logf (int category, int priority, const char *fmt, ...)
{
  va_list             ap;

  va_start (ap, fmt);
  sc_logv ("unknown", -1, t8dg_package_id, category, priority, fmt, ap);
  va_end (ap);
}

void
t8dg_log_indent_push (void)
{
  sc_log_indent_push_count (t8dg_get_package_id (), 1);
}

void
t8dg_log_indent_pop (void)
{
  sc_log_indent_pop_count (t8dg_get_package_id (), 1);
}

void
t8dg_global_errorf (const char *fmt, ...)
{
  va_list             ap;

  va_start (ap, fmt);
  t8dg_logv (SC_LC_GLOBAL, SC_LP_ERROR, fmt, ap);
  va_end (ap);
}

void
t8dg_global_essentialf (const char *fmt, ...)
{
  va_list             ap;

  va_start (ap, fmt);
  t8dg_logv (SC_LC_GLOBAL, SC_LP_ESSENTIAL, fmt, ap);
  va_end (ap);
}

void
t8dg_global_productionf (const char *fmt, ...)
{
  va_list             ap;

  va_start (ap, fmt);
  t8dg_logv (SC_LC_GLOBAL, SC_LP_PRODUCTION, fmt, ap);
  va_end (ap);
}

void
t8dg_global_infof (const char *fmt, ...)
{
  va_list             ap;

  va_start (ap, fmt);
  t8dg_logv (SC_LC_GLOBAL, SC_LP_INFO, fmt, ap);
  va_end (ap);
}

void
t8dg_infof (const char *fmt, ...)
{
  va_list             ap;

  va_start (ap, fmt);
  t8dg_logv (SC_LC_NORMAL, SC_LP_INFO, fmt, ap);
  va_end (ap);
}

void
t8dg_debugf (const char *fmt, ...)
{
#ifdef T8_ENABLE_DEBUG
  va_list             ap;

  va_start (ap, fmt);
  t8dg_logv (SC_LC_NORMAL, SC_LP_DEBUG, fmt, ap);
  va_end (ap);
#endif
}

void
t8dg_errorf (const char *fmt, ...)
{
  va_list             ap;

  va_start (ap, fmt);
  t8dg_logv (SC_LC_NORMAL, SC_LP_ERROR, fmt, ap);
  va_end (ap);
}

void
t8dg_init (int log_threshold)
{
  t8dg_package_id = sc_package_register (NULL, log_threshold, "t8dg", "t8DG solver");
}

t8_locidx_t
t8dg_itree_ielement_to_idata (t8_forest_t forest, t8_locidx_t itree, t8_locidx_t ielement)
{
  return t8_forest_get_tree_element_offset (forest, itree) + ielement;
}

void
t8dg_vec_print (double x[3])
{
  t8dg_debugf ("%f %f %f\n", x[0], x[1], x[2]);
}

void               *
t8dg_sc_array_index_locidx (const sc_array_t * array, t8dg_locidx_t it)
{
  T8DG_ASSERT (it >= 0 && (size_t) it < array->elem_count);
  return array->array + array->elem_size * (size_t) it;
}

void
t8dg_transform_3tensoridx (int idx, const int tensordims[DIM3], int tensoridx[DIM3])
{
  /*TODO: Check */
  int                 numtensor = 0, itensor;
  while (numtensor < 3 && tensordims[numtensor] > 0) {
    numtensor++;
  }
  for (itensor = 0; itensor < numtensor; itensor++) {
    tensoridx[itensor] = idx % tensordims[itensor];
    idx /= tensordims[itensor];
  }

}
