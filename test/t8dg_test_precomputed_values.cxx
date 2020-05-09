#include <gtest/gtest.h>
#include "../src/t8dg.h"
#include "../src/t8dg_global_precomputed_values.h"
#include "../src/t8dg_local_precomputed_values.h"
#include "../src/t8dg_precomputed_values.h"
#include "../src/t8dg_sc_array.h"
#include "../src/t8dg_coarse_geometry.h"
#include <example/common/t8_example_common.h>

#include <t8_eclass.h>
#include <t8_cmesh.h>
#include <t8_forest.h>
#include <t8_schemes/t8_default_cxx.hxx>
#include <mpi.h>

/* *INDENT-OFF* */
class PrecomputedValuesProjectionTest1D : public ::testing::TestWithParam<int>
/* *INDENT-ON* */
{
protected:
  /* *INDENT-OFF* */
  void SetUp () override
  /* *INDENT-ON* */
  {
    int                 num_lgl;
    t8_cmesh_t          cmesh;
    t8_scheme_cxx_t    *default_scheme;
    t8dg_quadrature_t  *quadrature;

    cmesh = t8_cmesh_new_hypercube (T8_ECLASS_LINE, sc_MPI_COMM_WORLD, 0, 0, 1);
    default_scheme = t8_scheme_new_default_cxx ();
    t8_scheme_cxx_ref (default_scheme);
    t8_cmesh_ref (cmesh);
    forest_adapt = t8_forest_new_uniform (cmesh, default_scheme, 1, 1, sc_MPI_COMM_WORLD);
    forest = t8_forest_new_uniform (cmesh, default_scheme, 0, 1, sc_MPI_COMM_WORLD);

    coarse_geometry = t8dg_coarse_geometry_new_1D_linear ();
    t8dg_geometry_transformation_data_t geometry_data = { coarse_geometry, forest, 0, 0 };
    t8dg_geometry_transformation_data_t geometry_data_adapt = { coarse_geometry, forest_adapt, 0, 0 };

    num_lgl = GetParam ();

    global_values = t8dg_global_precomputed_values_new_1D_LGL (num_lgl);
    quadrature = t8dg_global_precomputed_values_get_quadrature (global_values);
    local_values = t8dg_local_precomputed_values_new (quadrature, 1);
    local_values_adapt = t8dg_local_precomputed_values_new (quadrature, 2);
    t8dg_local_precomputed_values_set_element (local_values, &geometry_data, quadrature);
    t8dg_local_precomputed_values_set_element (local_values_adapt, &geometry_data_adapt, quadrature);
    geometry_data_adapt.ielement = 1;
    t8dg_local_precomputed_values_set_element (local_values_adapt, &geometry_data_adapt, quadrature);

    dof_values_adapt = sc_array_new_count (sizeof (double) * num_lgl, 2);
    dof_values = sc_array_new_count (sizeof (double) * num_lgl, 1);

  }

  /* *INDENT-OFF* */
  void TearDown () override
  /* *INDENT-ON* */
  {
    t8dg_global_precomputed_values_destroy (&global_values);
    t8dg_local_precomputed_values_destroy (&local_values);
    t8dg_local_precomputed_values_destroy (&local_values_adapt);
    t8dg_coarse_geometry_destroy (&coarse_geometry);
    sc_array_destroy (dof_values_adapt);
    sc_array_destroy (dof_values);
    t8_forest_unref (&forest_adapt);
    t8_forest_unref (&forest);
  }

  t8dg_global_precomputed_values_t *global_values;
  t8dg_local_precomputed_values_t *local_values;
  t8dg_local_precomputed_values_t *local_values_adapt;
  t8dg_coarse_geometry_t *coarse_geometry;

  sc_array_t         *dof_values_adapt;
  sc_array_t         *dof_values;
  t8_forest_t         forest, forest_adapt;
};

INSTANTIATE_TEST_SUITE_P (lglRange, PrecomputedValuesProjectionTest1D, testing::Range (2, MAX_LGL_NUMBER + 1));

TEST_P (PrecomputedValuesProjectionTest1D, const_one)
{
  int                 num_lgl, idof, ichild;
  sc_array_t         *element_dof_values_child[2];
  sc_array_t         *element_dof_values_parent;

  num_lgl = GetParam ();
  for (ichild = 0; ichild < 2; ichild++) {
    element_dof_values_child[ichild] = t8dg_sc_array_block_double_new_view (dof_values_adapt, ichild);
    for (idof = 0; idof < num_lgl; idof++) {
      *(double *) sc_array_index (element_dof_values_child[ichild], idof) = 1;  //idof + ichild * (i_lgl + 1);
    }
    sc_array_destroy (element_dof_values_child[ichild]);
  }

  for (ichild = 0; ichild < 2; ichild++) {
    element_dof_values_child[ichild] = t8dg_sc_array_block_double_new_view (dof_values_adapt, ichild);
  }
  element_dof_values_parent = t8dg_sc_array_block_double_new_view (dof_values, 0);
  t8dg_precomputed_values_transform_child_dof_to_parent_dof (global_values, element_dof_values_child, element_dof_values_parent,
                                                             2, local_values_adapt, local_values, 0, 0);

  for (idof = 0; idof < num_lgl; idof++) {
    EXPECT_NEAR (*(double *) sc_array_index_int (element_dof_values_parent, idof), 1, 1e-10);
  }
  for (ichild = 0; ichild < 2; ichild++) {
    sc_array_destroy (element_dof_values_child[ichild]);
  }
  sc_array_destroy (element_dof_values_parent);
}

/* *INDENT-OFF* */
class PrecomputedValuesL2norm1D : public ::testing::TestWithParam<std::tuple<int,std::tuple<t8dg_scalar_function_3d_time_fn,double>>>
/* *INDENT-ON* */

{
protected:
  /* *INDENT-OFF* */
  void SetUp () override
  /* *INDENT-ON* */
  {
    int                 num_lgl;
    int                 level = 4;
    int                 itree, ielement;
    t8_cmesh_t          cmesh;
    t8_scheme_cxx_t    *default_scheme;

    t8dg_quadrature_t  *quadrature;

    cmesh = t8_cmesh_new_periodic_line_more_trees (sc_MPI_COMM_WORLD);
    default_scheme = t8_scheme_new_default_cxx ();
    forest = t8_forest_new_uniform (cmesh, default_scheme, level, 1, sc_MPI_COMM_WORLD);
    coarse_geometry = t8dg_coarse_geometry_new_1D_linear ();

    t8dg_geometry_transformation_data_t geometry_data = { coarse_geometry, forest, 0, 0 };

    num_lgl = std::get < 0 > (GetParam ());

    global_values = t8dg_global_precomputed_values_new_1D_LGL (num_lgl);
    quadrature = t8dg_global_precomputed_values_get_quadrature (global_values);
    local_values = t8dg_local_precomputed_values_new (quadrature, t8_forest_get_num_element (forest));
    for (itree = 0; itree < t8_forest_get_num_local_trees (forest); itree++) {
      geometry_data.itree = itree;
      for (ielement = 0; ielement < t8_forest_get_tree_num_elements (forest, itree); ielement++) {
        geometry_data.ielement = ielement;
        t8dg_local_precomputed_values_set_element (local_values, &geometry_data, quadrature);
      }
    }

    dof_values = sc_array_new_count (sizeof (double) * num_lgl, t8_forest_get_num_element (forest));
  }

  /* *INDENT-OFF* */
  void TearDown () override
  /* *INDENT-ON* */
  {
    t8dg_global_precomputed_values_destroy (&global_values);
    t8dg_local_precomputed_values_destroy (&local_values);
    t8dg_coarse_geometry_destroy (&coarse_geometry);
    sc_array_destroy (dof_values);
    t8_forest_unref (&forest);
  }

  t8dg_global_precomputed_values_t *global_values;
  t8dg_local_precomputed_values_t *local_values;
  t8dg_coarse_geometry_t *coarse_geometry;

  sc_array_t         *dof_values;
  t8_forest_t         forest;
};

static double
t8dg_test_const_one (const double x[3], const double t)
{
  return 1;
}

static double
t8dg_test_sinx (const double x[3], const double t)
{
  return sin (2 * M_PI * x[0]);
}

static double
t8dg_test_expx (const double x[3], const double t)
{
  return exp (x[0]);
}

static const
  std::tuple <
t8dg_scalar_function_3d_time_fn, double >
function1 (t8dg_test_const_one, 1);
static const
  std::tuple <
t8dg_scalar_function_3d_time_fn, double >function2 (t8dg_test_sinx, sqrt (0.5));
static const
  std::tuple <
t8dg_scalar_function_3d_time_fn, double >function3 (t8dg_test_expx, sqrt ((exp (2) - 1) / 2));

/* *INDENT-OFF* */
INSTANTIATE_TEST_SUITE_P (lgl_and_functions, PrecomputedValuesL2norm1D,
    ::testing::Combine (
	::testing::Range (2, 8),
	::testing::Values (function1, function2, function3)));
/* *INDENT-ON* */

TEST_P (PrecomputedValuesL2norm1D, test_functions)
{
  t8dg_functionbasis_t *functionbasis;
  t8_locidx_t         itree = 0, ielement, idata;
  double              time = 0;
  double              norm_squared = 0, norm;
  sc_array_t         *element_dof_values;

  t8dg_geometry_transformation_data_t geometry_data = { coarse_geometry, forest, 0, 0 };
  t8dg_precomputed_values_fn_evaluation_data_t evaluation_data = { &geometry_data, std::get < 0 > (std::get < 1 > (GetParam ())), time };

  functionbasis = t8dg_global_precomputed_values_get_functionbasis (global_values);

  for (itree = 0, idata = 0; itree < t8_forest_get_num_local_trees (forest); itree++) {
    geometry_data.itree = itree;
    for (ielement = 0; ielement < t8_forest_get_tree_num_elements (forest, itree); ielement++, idata++) {
      geometry_data.ielement = ielement;

      element_dof_values = t8dg_sc_array_block_double_new_view (dof_values, idata);
      t8dg_functionbasis_interpolate_scalar_fn (functionbasis, t8dg_precomputed_values_transform_reference_vertex_and_evaluate,
                                                &evaluation_data, element_dof_values);
      norm_squared += t8dg_precomputed_values_element_norm_l2_squared (element_dof_values, global_values, local_values, idata);
      sc_array_destroy (element_dof_values);
    }
  }
  norm = sqrt (norm_squared);
  EXPECT_NEAR (norm, std::get < 1 > (std::get < 1 > (GetParam ())), exp (-4 * std::get < 0 > (GetParam ())));
}
