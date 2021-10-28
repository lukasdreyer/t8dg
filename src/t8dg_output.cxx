#include "t8dg_dof.h"
#include "t8dg_values.h"
#include <t8_forest.h>
#include <t8_forest_vtk.h>
#include <t8_vtk.h>
#include <t8dg_advect_diff_problem.h>
#include "t8dg_output.h"

t8dg_vtk_data_t    *
t8dg_output_vtk_data_new (const char *prefix, int vtk_freq)
{
  t8dg_vtk_data_t    *vtk_data;
  vtk_data = T8DG_ALLOC_ZERO (t8dg_vtk_data_t, 1);
  vtk_data->prefix = prefix;
  vtk_data->vtk_count = 0;
  vtk_data->vtk_freq = vtk_freq;
  return vtk_data;
}

void
t8dg_output_vtk_data_destroy (t8dg_vtk_data_t ** p_vtk_data)
{
  t8dg_vtk_data_t    *vtk_data;
  vtk_data = *p_vtk_data;
  T8DG_FREE (vtk_data);
  vtk_data = NULL;
}

void
t8dg_output_write_vtk (const t8dg_dof_values_t * dof_values, t8dg_vtk_data_t * output_data, int write_flow, t8dg_linear_flux3D_fn flow_field, const double time, t8dg_flux_data_base *flux_data)
{
  t8_forest_t         forest;
  double             *dof_array;
  double             *flow_array;
  t8_locidx_t         num_local_elements, itree, ielement, idata, num_trees, num_elems_in_tree;
  t8_vtk_data_field_t vtk_data[2];
  double              average;
  int                 idof, number_of_dof;
  t8dg_element_dof_values_t *element_dof_values;
  char                fileprefix[BUFSIZ];
  const int           num_data_fields = write_flow ? 2 : 1;

  forest = t8dg_dof_values_get_forest (dof_values);

  num_local_elements = t8_forest_get_local_num_elements (forest);
  dof_array = T8_ALLOC_ZERO (double, num_local_elements);
  if (write_flow) {
    flow_array = T8_ALLOC (double, 3 * num_local_elements);
  }

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
      if (write_flow) {
        /* Calculate flow at element midpoint. */
        double midpoint[3];
        const t8_element_t *element = t8_forest_get_element_in_tree (forest, itree, ielement);
        t8_forest_element_centroid (forest, itree, element, midpoint);
        flow_field (midpoint, flow_array + 3 * idata, time, flux_data, itree, ielement);
      }
    }
  }

  /* Write meta data for vtk */
  snprintf (vtk_data[0].description, BUFSIZ, "Num. Solution");
  vtk_data[0].type = T8_VTK_SCALAR;
  vtk_data[0].data = dof_array;
  if (write_flow) {
    snprintf (vtk_data[1].description, BUFSIZ, "Flow field");
    vtk_data[1].type = T8_VTK_VECTOR;
    vtk_data[1].data = flow_array;
  }
  /* Write filename */
  snprintf (fileprefix, BUFSIZ, "%s_%03i", output_data->prefix, output_data->vtk_count);
  /* Write vtk files */
  if (t8_forest_vtk_write_file (forest, fileprefix, 1, 1, 1, 1, 0, num_data_fields, vtk_data)) {
    t8_debugf ("[Advect] Wrote pvtu to files %s\n", fileprefix);
  }
  else {
    t8_errorf ("[Advect] Error writing to files %s\n", fileprefix);
  }
  /*increase vtk_count */
  (output_data->vtk_count)++;
  /* clean-up */
  T8_FREE (dof_array);
  if (write_flow) {
    T8_FREE (flow_array);
  }
}
