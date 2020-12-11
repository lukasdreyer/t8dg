#include "t8dg_dof.h"
#include "t8dg_values.h"
#include <t8_forest.h>
#include <t8_forest_vtk.h>
#include <t8_vtk.h>

void
t8dg_output_write_vtk (t8dg_dof_values_t * dof_values, t8_forest_t forest, char prefix[BUFSIZ], int *vtk_count)
{
  double             *dof_array;
  t8_locidx_t         num_local_elements, itree, ielement, idata, num_trees, num_elems_in_tree;
  t8_vtk_data_field_t vtk_data;
  double              average;
  int                 idof, number_of_dof;
  t8dg_element_dof_values_t *element_dof_values;
  char                fileprefix[BUFSIZ];

  num_local_elements = t8_forest_get_num_element (forest);
  dof_array = T8_ALLOC_ZERO (double, num_local_elements);

  num_trees = t8_forest_get_num_local_trees (forest);
  for (itree = 0, idata = 0; itree < num_trees; itree++) {
    num_elems_in_tree = t8_forest_get_tree_num_elements (forest, itree);
    for (ielement = 0; ielement < num_elems_in_tree; ielement++, idata++) {
      element_dof_values = t8dg_dof_values_new_element_dof_values_view (dof_values, itree, ielement);
      average = 0;
      number_of_dof = t8dg_element_dof_values_get_num_dof (element_dof_values);
      for (idof = 0; idof < number_of_dof; idof++) {
        average += t8dg_element_dof_values_get_value (element_dof_values, idof);
      }
      t8dg_element_dof_values_destroy (&element_dof_values);
      average /= number_of_dof;
      dof_array[idata] = average;
    }
  }

  /* Write meta data for vtk */
  snprintf (vtk_data.description, BUFSIZ, "Num. Solution");
  vtk_data.type = T8_VTK_SCALAR;
  vtk_data.data = dof_array;
  /* Write filename */
  snprintf (fileprefix, BUFSIZ, "%s_%03i", prefix, *vtk_count);
  /* Write vtk files */
  if (t8_forest_vtk_write_file (forest, fileprefix, 1, 1, 1, 1, 0, 1, &vtk_data)) {
    t8_debugf ("[Advect] Wrote pvtu to files %s\n", fileprefix);
  }
  else {
    t8_errorf ("[Advect] Error writing to files %s\n", fileprefix);
  }
  /*increase vtk_count */
  (*vtk_count)++;
  /* clean-up */
  T8_FREE (dof_array);
}
