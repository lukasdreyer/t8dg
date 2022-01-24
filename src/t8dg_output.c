/*
 *  This file is part of t8dg.
 *  t8dg is a discontinuous galerkin solver on adaptive meshes
 *  that uses the adaptive meshing library t8code to distribute
 *  the workload equally among parallel processes to achieve good scaling.
 *
 *  Copyright (C) 2020 the developers
 *
 *  t8dg is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  t8dg is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with t8dg.  If not, see <http://www.gnu.org/licenses/>.
 */

#include "t8dg_dof.h"
#include "t8dg_values.h"
#include <t8_forest.h>
#include <t8_forest_vtk.h>
#include <t8_vtk.h>
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
t8dg_output_write_vtk (t8dg_dof_values_t * dof_values, t8dg_vtk_data_t * output_data)
{
  t8_forest_t         forest;
  double             *dof_array;
  t8_locidx_t         num_local_elements, itree, ielement, idata, num_trees, num_elems_in_tree;
  t8_vtk_data_field_t vtk_data;
  double              average;
  int                 idof, number_of_dof;
  t8dg_element_dof_values_t *element_dof_values;
  char                fileprefix[BUFSIZ];

  forest = t8dg_dof_values_get_forest (dof_values);

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
  snprintf (fileprefix, BUFSIZ, "%s_%03i", output_data->prefix, output_data->vtk_count);
  /* Write vtk files */
  if (t8_forest_vtk_write_file (forest, fileprefix, 1, 1, 1, 1, 0, 1, &vtk_data)) {
    t8_debugf ("[Advect] Wrote pvtu to files %s\n", fileprefix);
  }
  else {
    t8_errorf ("[Advect] Error writing to files %s\n", fileprefix);
  }
  /*increase vtk_count */
  (output_data->vtk_count)++;
  /* clean-up */
  T8_FREE (dof_array);
}
