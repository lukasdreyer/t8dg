#include <t8_forest.h>
#include "t8dg_dof.h"
#include "t8dg_values.h"
#include "t8dg_adapt.h"
#include <t8_element_cxx.hxx>

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
    if (level == adapt_data->uniform_refinement_level) {
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
t8dg_adapt_gradient (t8_forest_t forest,
                     t8_forest_t forest_from,
                     t8_locidx_t itree, t8_locidx_t lelement_id, t8_eclass_scheme_c * ts, int num_elements, t8_element_t * elements[])
{
  t8dg_adapt_data_t  *adapt_data;
  t8_locidx_t         first_idata;
  double             *dof_values;
  int                 level;
  double              diam;
  double             *tree_vertices;

  double              gradient_threshold_refine = 0.6;
  double              gradient_threshold_coarsen = 0.2;

  int                 num_dof;

  first_idata = t8dg_itree_ielement_to_idata (forest_from, itree, lelement_id);
  adapt_data = (t8dg_adapt_data_t *) t8_forest_get_user_data (forest);

  T8DG_CHECK_ABORT (t8_eclass_to_dimension[t8_forest_get_eclass (forest, itree)] == 1, "Not yet implemented");

  num_dof = t8dg_global_values_get_num_dof (t8dg_values_get_global_values (adapt_data->dg_values, itree, lelement_id));

  level = ts->t8_element_level (elements[0]);

  if (level == adapt_data->maximum_refinement_level && num_elements == 1) {
    /* It is not possible to refine this level */
    return 0;
  }
  tree_vertices = t8_forest_get_tree_vertices (forest_from, itree);

  if (num_elements == 1) {
    dof_values = t8dg_dof_values_get_double_pointer (adapt_data->dof_values, first_idata);
    diam = t8_forest_element_diam (forest_from, itree, elements[0], tree_vertices);

    double              gradient = fabs (dof_values[0] - dof_values[num_dof - 1]) / diam;
    return gradient > gradient_threshold_refine;
  }
  else {
    dof_values = t8dg_dof_values_get_double_pointer (adapt_data->dof_values, first_idata);
    diam = t8_forest_element_diam (forest_from, itree, elements[0], tree_vertices);

    double              gradient_left = fabs (dof_values[0] - dof_values[num_dof - 1]) / diam;

    dof_values = t8dg_dof_values_get_double_pointer (adapt_data->dof_values, first_idata);
    diam = t8_forest_element_diam (forest_from, itree, elements[1], tree_vertices);

    double              gradient_right = fabs (dof_values[0] - dof_values[num_dof - 1]) / diam;
    if (gradient_left > gradient_threshold_refine && level < adapt_data->maximum_refinement_level)
      return 1;
    return -(gradient_left < gradient_threshold_coarsen && gradient_right < gradient_threshold_coarsen
             && level > adapt_data->uniform_refinement_level);
  }
  return 0;
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
