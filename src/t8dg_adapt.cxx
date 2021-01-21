#include <t8_forest.h>
#include "t8dg_dof.h"
#include "t8dg_values.h"
#include "t8dg_adapt.h"
#include <t8_element_cxx.hxx>

t8_forest_adapt_t
t8dg_adapt_fn_arg (int adapt_arg)
{
  switch (adapt_arg) {
  case 0:
    return t8dg_adapt_mass;
    break;
  case 1:
    return t8dg_adapt_rel_min_max;
    break;

  default:
    T8DG_ABORT ("Wrong adapt fn arg");
    break;
  }
}

t8dg_adapt_data_t  *
t8dg_adapt_data_new (t8dg_values_t * dg_values, int initial_level, int min_level, int max_level, int adapt_fn_arg, int adapt_freq)
{
  t8dg_adapt_data_t  *adapt_data;
  adapt_data = T8DG_ALLOC_ZERO (t8dg_adapt_data_t, 1);
  adapt_data->dg_values = dg_values;
  adapt_data->minimum_refinement_level = min_level;
  adapt_data->initial_refinement_level = initial_level;
  adapt_data->maximum_refinement_level = max_level;
  adapt_data->adapt_fn = t8dg_adapt_fn_arg (adapt_fn_arg);
  adapt_data->adapt_freq = adapt_freq;
  return adapt_data;
}

void
t8dg_adapt_data_destroy (t8dg_adapt_data_t ** p_adapt_data)
{
  t8dg_adapt_data_t  *adapt_data;
  adapt_data = *p_adapt_data;
  T8DG_FREE (adapt_data);
  p_adapt_data = NULL;
}

int
t8dg_adapt_mass (t8_forest_t forest,
                 t8_forest_t forest_from,
                 t8_locidx_t itree, t8_locidx_t ielement, t8_eclass_scheme_c * ts, int num_elements, t8_element_t * elements[])
{
  t8dg_adapt_data_t  *adapt_data;
  sc_array_t         *element_dof;
  int                 level;
  double              norm, area;
  int                 ifamilyelement;

  adapt_data = (t8dg_adapt_data_t *) t8_forest_get_user_data (forest);
  level = ts->t8_element_level (elements[0]);

  if (level == adapt_data->maximum_refinement_level && num_elements == 1) {
    /* It is not possible to refine this level */
    return 0;
  }

  element_dof = t8dg_dof_values_new_element_dof_values_view (adapt_data->dof_values, itree, ielement);

  norm = t8dg_values_element_norm_l2_squared (adapt_data->dg_values, element_dof, itree, ielement);
  area = t8dg_values_element_area (adapt_data->dg_values, itree, ielement);

  sc_array_destroy (element_dof);

  if (norm / area > 0.2) {
    return level < adapt_data->maximum_refinement_level;
  }
  if (num_elements > 1) {
    if (level == adapt_data->minimum_refinement_level) {
      return 0;                 /* It is not possible to coarsen this element. If this is wanted, balance is needed outside */
    }

    for (ifamilyelement = 0; ifamilyelement < num_elements; ifamilyelement++) {
      element_dof = t8dg_dof_values_new_element_dof_values_view (adapt_data->dof_values, itree, ielement + ifamilyelement);
      norm = t8dg_values_element_norm_l2_squared (adapt_data->dg_values, element_dof, itree, ielement + ifamilyelement);
      area = t8dg_values_element_area (adapt_data->dg_values, itree, ielement + ifamilyelement);
      sc_array_destroy (element_dof);

      if (norm / area > 0.1) {
        return 0;
      }
    }
    return -1;
  }
  return 0;
}

int
t8dg_adapt_rel_min_max (t8_forest_t forest,
                        t8_forest_t forest_from,
                        t8_locidx_t itree, t8_locidx_t lelement_id, t8_eclass_scheme_c * ts, int num_elements, t8_element_t * elements[])
{
  t8dg_adapt_data_t  *adapt_data;
  int                 level;

  double              max_value;
  double              min_value;
  double              rel_gradient;
  /*Move thresholds into adapt_data? */
  double              gradient_threshold_refine = 0.6;
  double              gradient_threshold_coarsen = 0.5;
  double              min_value_threshold = 0.05;

  adapt_data = (t8dg_adapt_data_t *) t8_forest_get_user_data (forest);

  level = ts->t8_element_level (elements[0]);

  if (level < adapt_data->maximum_refinement_level) {
    max_value = t8dg_dof_values_get_max_value (adapt_data->dof_values, itree, lelement_id);
    min_value = t8dg_dof_values_get_min_value (adapt_data->dof_values, itree, lelement_id);
    if (min_value <= min_value_threshold) {
      if (max_value > min_value_threshold)
        return 1;
    }
    else {
      rel_gradient = (max_value - min_value) / min_value;
      if (rel_gradient > gradient_threshold_refine) {
        return 1;
      }
    }
  }
  //not refine, check coarsen
  if (num_elements == 1 || level <= adapt_data->minimum_refinement_level) {
    return 0;
  }
  else {
    int                 ichild;
    for (ichild = 1; ichild < num_elements; ichild++) {
      max_value = t8dg_dof_values_get_max_value (adapt_data->dof_values, itree, lelement_id + ichild);
      min_value = t8dg_dof_values_get_min_value (adapt_data->dof_values, itree, lelement_id + ichild);
      if (min_value <= min_value_threshold && max_value > min_value_threshold)
        return 0;
      rel_gradient = (max_value - min_value) / min_value;
      if (rel_gradient > gradient_threshold_coarsen)
        return 0;
    }
    return -1;
  }
}

void
t8dg_adapt_replace (t8_forest_t forest_old,
                    t8_forest_t forest_new,
                    t8_locidx_t itree,
                    t8_eclass_scheme_c * ts, int num_elems_old, t8_locidx_t first_ielem_old, int num_elems_new, t8_locidx_t first_ielem_new)
{
  t8dg_adapt_data_t  *adapt_data;
  t8_locidx_t         first_idata_old, first_idata_new;
  int                 ichild, num_children;

  adapt_data = (t8dg_adapt_data_t *) t8_forest_get_user_data (forest_new);

  first_idata_old = t8dg_itree_ielement_to_idata (forest_old, itree, first_ielem_old);
  first_idata_new = t8dg_itree_ielement_to_idata (forest_new, itree, first_ielem_new);

  if (num_elems_old == num_elems_new && num_elems_old == 1) {
    t8dg_values_copy_element_values (adapt_data->dg_values, first_idata_old, first_idata_new);
    t8dg_dof_values_copy_from_index_to_index (adapt_data->dof_values, first_idata_old, adapt_data->dof_values_adapt, first_idata_new);
  }
  else if (num_elems_old == 1) {
    num_children = t8dg_global_values_get_num_children (t8dg_values_get_global_values (adapt_data->dg_values, itree, first_ielem_old));
    T8DG_ASSERT (num_children == num_elems_new);
    for (ichild = 0; ichild < num_children; ichild++) {
      t8dg_values_set_element_adapt (adapt_data->dg_values, itree, first_ielem_new + ichild);
      t8dg_values_transform_parent_dof_to_child_dof (adapt_data->dg_values, adapt_data->dof_values, adapt_data->dof_values_adapt, itree,
                                                     first_ielem_old, first_ielem_new + ichild, ichild);
    }
  }
  else {
    num_children =
      t8dg_global_values_get_num_children (t8dg_values_get_global_values_adapt (adapt_data->dg_values, itree, first_ielem_new));
    T8DG_ASSERT (num_children == num_elems_old && num_elems_new == 1);
    /* Needed before! we calculate the dof_values */
    t8dg_values_set_element_adapt (adapt_data->dg_values, itree, first_ielem_new);
    t8dg_values_transform_child_dof_to_parent_dof (adapt_data->dg_values, adapt_data->dof_values, adapt_data->dof_values_adapt, itree,
                                                   num_children, first_ielem_old, first_ielem_new);
  }
}
