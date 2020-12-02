#include <gtest/gtest.h>
#include "../src/t8dg.h"
#include "../src/t8dg_mortar.h"
#include <sc_containers.h>

#include <t8_cmesh.h>
#include <t8_schemes/t8_default_cxx.hxx>
#include <t8_forest.h>

static double const x_unit_vec[3] = { 1, 0, 0 };

TEST (mortar_arrray1D, creation)
{
/*
  t8dg_mortar_array_t *mortar_array;

  t8_forest_t         forest;
  t8_cmesh_t          cmesh;
  t8_scheme_cxx_t    *default_scheme;

  double              dof_values[6] = { 0, 2, 4, 7, 4, 1 };
  sc_array_t         *dof_array;

  t8dg_global_precomputed_values_t *global_values;
  t8dg_local_precomputed_values_t *local_values;

  t8dg_coarse_geometry_t *coarse_geometry;

  t8dg_flux_t        *flux;

  cmesh = t8_cmesh_new_hypercube (T8_ECLASS_LINE, sc_MPI_COMM_WORLD, 0, 0, 1);
  default_scheme = t8_scheme_new_default_cxx ();
  forest = t8_forest_new_uniform (cmesh, default_scheme, 1, 1, sc_MPI_COMM_WORLD);
  coarse_geometry = t8dg_coarse_geometry_new_1D_linear ();
  flux = t8dg_flux_new_linear_constant_flux (x_unit_vec, 1);

  mortar_array = t8dg_mortar_array_new_empty (forest, 2);

  dof_array = sc_array_new_data (dof_values, sizeof (double) * 3, 2);

  global_values = t8dg_global_precomputed_values_new_1D_LGL (3);
  local_values = t8dg_local_precomputed_values_new (forest, global_values);

  t8dg_geometry_transformation_data_t geometry_data = { coarse_geometry, forest, 0, 0 };

  t8dg_local_precomputed_values_set_element (local_values, &geometry_data, global_values);
  geometry_data.ielement = 1;
  t8dg_local_precomputed_values_set_element (local_values, &geometry_data, global_values);

  t8dg_mortar_fill_data_t mortar_fill_data = { global_values, local_values, &geometry_data, flux, dof_array, 0 };

  t8dg_mortar_array_fill (mortar_array, &mortar_fill_data);

  EXPECT_DOUBLE_EQ (*(double *) sc_array_index (t8dg_mortar_array_get_oriented_flux (mortar_array, 0, 0), 0), -1);
  EXPECT_DOUBLE_EQ (*(double *) sc_array_index (t8dg_mortar_array_get_oriented_flux (mortar_array, 0, 1), 0), 4);
  EXPECT_DOUBLE_EQ (*(double *) sc_array_index (t8dg_mortar_array_get_oriented_flux (mortar_array, 1, 0), 0), -4);
  EXPECT_DOUBLE_EQ (*(double *) sc_array_index (t8dg_mortar_array_get_oriented_flux (mortar_array, 1, 1), 0), 1);

  t8dg_mortar_array_destroy (&mortar_array);
  t8dg_flux_destroy (&flux);
  t8dg_local_precomputed_values_destroy (&local_values);
  t8dg_global_precomputed_values_destroy (&global_values);
  t8dg_coarse_geometry_destroy (&coarse_geometry);
  sc_array_destroy_null (&dof_array);
  t8_forest_unref (&forest);
*/
}
