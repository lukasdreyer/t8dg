/*
 * t8dg.c
 *
 *  Created on: Mar 11, 2020
 *      Author: lukas
 */

#include <t8_forest.h>
#include <t8_forest/t8_forest_ghost.h>
#include <t8_element_cxx.hxx>
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

t8_locidx_t
t8dg_ighosttree_ielement_to_idata (t8_forest_t forest, t8_locidx_t ighosttree, t8_locidx_t ielement)
{
  return t8_forest_get_local_num_elements (forest) + t8_forest_ghost_get_tree_element_offset (forest, ighosttree) + ielement;
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

int
t8dg_almost_equal (double x1, double x2)
{
  double              eps = 1e-10;
  return fabs (x1 - x2) < eps;
}

t8_eclass_t
t8dg_forest_get_eclass (t8_forest_t forest, t8_locidx_t itree, t8_locidx_t ielement)
{
  if (itree == -1) {
    /*TODO: is there a way to find out the eclass of a ghost element */
    return t8_forest_get_eclass (forest, 0);
  }
  return t8_forest_get_eclass (forest, itree);
}

t8_eclass_t
t8dg_eclass_from_gloidx_element (const t8_forest_t forest, const t8_gloidx_t iglobaltree, const t8_element_t * element)
{
  t8_eclass_t         tree_eclass;
  t8_eclass_scheme_c *scheme;
  t8_locidx_t         ilocaltree;
  t8_locidx_t         ighosttree;
  ilocaltree = t8_forest_get_local_id (forest, iglobaltree);
  if (ilocaltree >= 0) {
    tree_eclass = t8_forest_get_eclass (forest, ilocaltree);
  }
  else {
    ighosttree = t8_forest_ghost_get_ghost_treeid (forest, iglobaltree);
    T8DG_ASSERT (ighosttree >= 0);
    tree_eclass = t8_forest_ghost_get_tree_class (forest, ighosttree);
  }
  scheme = t8_forest_get_eclass_scheme (forest, tree_eclass);
  return scheme->t8_element_shape (element);
}

double             *
t8dg_forest_get_tree_vertices_gloidx (t8_forest_t forest, t8_gloidx_t iglobaltree)
{
  t8_locidx_t         localcmeshtreeid;
  localcmeshtreeid = t8_cmesh_get_local_id (forest->cmesh, iglobaltree);
  return t8_cmesh_get_tree_vertices (forest->cmesh, localcmeshtreeid);
}
