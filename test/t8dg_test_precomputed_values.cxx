#include <gtest/gtest.h>
#include "../src/t8dg.h"
#include "../src/t8dg_global_precomputed_values.h"
#include "../src/t8dg_local_precomputed_values.h"
#include "../src/t8dg_precomputed_values.h"
#include "../src/t8dg_sc_array.h"

#include <t8_cmesh.h>
#include <t8_forest.h>
#include <t8_schemes/t8_default_cxx.hxx>
#include <mpi.h>

#define TEST_MAX_LGL_NUMBER 4
/* *INDENT-OFF* */
class PrecomputedValuesProjectionTest : public ::testing::Test
/* *INDENT-ON* */
{
protected:
  /* *INDENT-OFF* */
  void SetUp () override
  /* *INDENT-ON* */
  {
    int                 i_lgl;
    t8_cmesh_t          cmesh;
    t8_scheme_cxx_t    *default_scheme;

    t8_eclass_scheme_c *scheme;
    t8dg_quadrature_t  *quadrature;

    cmesh = t8_cmesh_new_hypercube (T8_ECLASS_LINE, sc_MPI_COMM_WORLD, 0, 0, 1);
    default_scheme = t8_scheme_new_default_cxx ();
    t8_scheme_cxx_ref (default_scheme);
    t8_cmesh_ref (cmesh);
    forest_adapt = t8_forest_new_uniform (cmesh, default_scheme, 1, 1, sc_MPI_COMM_WORLD);
    forest = t8_forest_new_uniform (cmesh, default_scheme, 0, 1, sc_MPI_COMM_WORLD);

    scheme = t8_forest_get_eclass_scheme (forest, T8_ECLASS_LINE);

    for (i_lgl = 0; i_lgl < TEST_MAX_LGL_NUMBER; i_lgl++) {
      global_values[i_lgl] = t8dg_global_precomputed_values_new_1D_LGL (i_lgl + 1);
      quadrature = t8dg_global_precomputed_values_get_quadrature (global_values[i_lgl]);
      local_values[i_lgl] = t8dg_local_precomputed_values_new (quadrature, 1);
      local_values_adapt[i_lgl] = t8dg_local_precomputed_values_new (quadrature, 2);
      t8dg_local_precomputed_values_set_element (local_values[i_lgl], forest, 0, scheme, 0, quadrature);
      t8dg_local_precomputed_values_set_element (local_values_adapt[i_lgl], forest_adapt, 0, scheme, 0, quadrature);
      t8dg_local_precomputed_values_set_element (local_values_adapt[i_lgl], forest_adapt, 0, scheme, 1, quadrature);

      dof_values_adapt[i_lgl] = sc_array_new_count (sizeof (double) * (i_lgl + 1), 2);
      dof_values[i_lgl] = sc_array_new_count (sizeof (double) * (i_lgl + 1), 1);
    }

  }

  /* *INDENT-OFF* */
  void TearDown () override
  /* *INDENT-ON* */
  {
    int                 i_lgl;
    for (i_lgl = 0; i_lgl < TEST_MAX_LGL_NUMBER; i_lgl++) {
      t8dg_global_precomputed_values_destroy (global_values + i_lgl);
      t8dg_local_precomputed_values_destroy (local_values + i_lgl);
      t8dg_local_precomputed_values_destroy (local_values_adapt + i_lgl);
      sc_array_destroy (dof_values_adapt[i_lgl]);
      sc_array_destroy (dof_values[i_lgl]);
    }
    t8_forest_unref (&forest_adapt);
    t8_forest_unref (&forest);
  }

  t8dg_global_precomputed_values_t *global_values[TEST_MAX_LGL_NUMBER];
  t8dg_local_precomputed_values_t *local_values[TEST_MAX_LGL_NUMBER];
  t8dg_local_precomputed_values_t *local_values_adapt[TEST_MAX_LGL_NUMBER];

  sc_array_t         *dof_values_adapt[TEST_MAX_LGL_NUMBER];
  sc_array_t         *dof_values[TEST_MAX_LGL_NUMBER];
  t8_forest_t         forest, forest_adapt;
};

TEST_F (PrecomputedValuesProjectionTest, const_one)
{
  int                 i_lgl, idof, ichild;
  sc_array_t         *element_dof_values_child[2];
  sc_array_t         *element_dof_values_parent;

  for (i_lgl = 0; i_lgl < TEST_MAX_LGL_NUMBER; i_lgl++) {
    for (ichild = 0; ichild < 2; ichild++) {
      element_dof_values_child[ichild] = t8dg_sc_array_block_double_new_view (dof_values_adapt[i_lgl], ichild);
      for (idof = 0; idof <= i_lgl; idof++) {
        *(double *) sc_array_index (element_dof_values_child[ichild], idof) = idof + ichild * (i_lgl + 1);
      }
      sc_array_destroy (element_dof_values_child[ichild]);
    }
  }

  for (i_lgl = 0; i_lgl < TEST_MAX_LGL_NUMBER; i_lgl++) {
    for (ichild = 0; ichild < 2; ichild++) {
      element_dof_values_child[ichild] = t8dg_sc_array_block_double_new_view (dof_values_adapt[i_lgl], ichild);
    }
    element_dof_values_parent = t8dg_sc_array_block_double_new_view (dof_values[i_lgl], 0);
    t8dg_precomputed_values_transform_child_dof_to_parent_dof (global_values[i_lgl], element_dof_values_child, element_dof_values_parent,
                                                               2, local_values_adapt[i_lgl], local_values[i_lgl], 0, 0);

    for (ichild = 0; ichild < 2; ichild++) {
      sc_array_destroy (element_dof_values_child[ichild]);
    }
    t8dg_sc_array_block_double_debug_print (element_dof_values_parent);
    sc_array_destroy (element_dof_values_parent);
  }
}
