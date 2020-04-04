/*
 * t8dg.c
 *
 *  Created on: Mar 11, 2020
 *      Author: lukas
 */

#include <t8_forest.h>
#include <t8.h>
#include "t8dg.h"

t8_locidx_t
t8dg_itree_ielement_to_idata (t8_forest_t forest, t8_locidx_t itree, t8_locidx_t ielement)
{
  return t8_forest_get_tree_element_offset (forest, itree) + ielement;
}

void
t8dg_vec_print (double x[3])
{
  t8_debugf ("%f %f %f\n", x[0], x[1], x[2]);
}
