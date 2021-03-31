#include <gtest/gtest.h>
#include "../src/t8dg.h"
#include "../src/t8dg_mortar.h"
#include "../src/t8dg_flux_implementation.h"
#include <sc_containers.h>

#include <t8_cmesh.h>
#include <t8_schemes/t8_default_cxx.hxx>
#include <t8_forest.h>

static double const x_unit_vec[3] = { 1, 0, 0 };

TEST (mortar_arrray1D, creation)
{
  t8dg_mortar_array_t *mortar_array;

  t8_forest_t         forest;
  t8_cmesh_t          cmesh;
  t8_scheme_cxx_t    *default_scheme;

  int                 num_lgl_vertices = 3;
  double              dof_array[6] = { 0, 2, 4, 7, 4, 1 };
  t8dg_dof_values_t  *dof_values;

  t8dg_local_values_t *local_values;
  t8dg_global_values_t *global_values;
  t8dg_global_values_t *global_values_array[T8_ECLASS_COUNT] = { 0 };

  t8dg_coarse_geometry_t *coarse_geometry;

  t8dg_linear_flux3D_fn linear_flux = t8dg_linear_flux3D_constant_flux_fn;
  t8dg_linear_flux3D_constant_flux_data_t flux_data;
  flux_data.flow_direction[0] = 1;
  flux_data.flow_direction[1] = 0;
  flux_data.flow_direction[2] = 0;
  flux_data.flow_velocity = 1;

  t8dg_numerical_linear_flux3D_fn numerical_flux = t8dg_linear_numerical_flux3D_lax_friedrich_fn;
  double              velocity_bound = 1;

  cmesh = t8_cmesh_new_hypercube (T8_ECLASS_LINE, sc_MPI_COMM_WORLD, 0, 0, 1);
  default_scheme = t8_scheme_new_default_cxx ();
  forest = t8_forest_new_uniform (cmesh, default_scheme, 1, 1, sc_MPI_COMM_WORLD);
  coarse_geometry = t8dg_coarse_geometry_new_1D_linear ();

  global_values = t8dg_global_values_new_1D_LGL (num_lgl_vertices);
  global_values_array[T8_ECLASS_LINE] = global_values;

  local_values = t8dg_local_values_new (forest, global_values_array, coarse_geometry);
  t8dg_local_values_set_all_elements (local_values);

  mortar_array = t8dg_mortar_array_new_empty (forest, local_values);
  dof_values = t8dg_dof_values_new_data_local (forest, global_values_array, dof_array, 6);

  t8dg_mortar_array_calculate_linear_flux3D (mortar_array, dof_values, linear_flux, &flux_data, numerical_flux, &velocity_bound, 0);

  /*Check locally */
  double              result[6] = { -1, 4, 0, -4, 1, 0 };
  double              flux_value;
  double              result_value;
  t8dg_dof_values_t  *result_dof;

  result_dof = t8dg_dof_values_new_data_local (forest, global_values_array, result, 2 * t8_forest_get_num_element (forest));
  t8dg_element_dof_values_t *element_result_dof;

  t8_locidx_t         num_elems_in_tree, ielement;
  int                 iface;
  num_elems_in_tree = t8_forest_get_tree_num_elements (forest, 0);
  for (ielement = 0; ielement < num_elems_in_tree; ielement++) {
    for (iface = 0; iface < 2; iface++) {
      flux_value = t8dg_face_dof_values_get_value (t8dg_mortar_array_get_oriented_flux (mortar_array, ielement, iface), 0);
      element_result_dof = t8dg_dof_values_new_element_dof_values_view (result_dof, 0, ielement);
      result_value = t8dg_element_dof_values_get_value (element_result_dof, iface);
#ifdef T8DG_ENABLE_MPI
      EXPECT_DOUBLE_EQ_MPI (flux_value, result_value);
#else
      EXPECT_DOUBLE_EQ (flux_value, result_value);
#endif
      t8dg_element_dof_values_destroy (&element_result_dof);
    }
  }
  t8dg_dof_values_destroy (&result_dof);
  t8dg_mortar_array_destroy (&mortar_array);
  t8dg_local_values_destroy (&local_values);
  t8dg_global_values_destroy (&global_values);
  t8dg_coarse_geometry_destroy (&coarse_geometry);
  t8dg_dof_values_destroy (&dof_values);
  t8_forest_unref (&forest);
}
