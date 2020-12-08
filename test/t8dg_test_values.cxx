#include <gtest/gtest.h>
#include "../src/t8dg.h"
#include "../src/t8dg_global_values.h"
#include "../src/t8dg_values.h"
#include "../src/t8dg_local_values.h"
#include "../src/t8dg_values.h"
#include "../src/t8dg_sc_array.h"
#include "../src/t8dg_coarse_geometry.h"
//#include <example/common/t8_example_common.h>

#include <t8_eclass.h>
#include <t8_cmesh.h>
#include <t8_forest.h>
#include <t8_forest/t8_forest_iterate.h>
#include <t8_schemes/t8_default_cxx.hxx>
#include <mpi.h>
#include <tuple>



/* *INDENT-OFF* */
class ValuesChildInterpolationTest2D : public ::testing::Test
/* *INDENT-ON* */
{
protected:
  /* *INDENT-OFF* */
  void SetUp () override
  /* *INDENT-ON* */
  {
    t8_cmesh_t          cmesh;
    t8_scheme_cxx_t    *default_scheme;

    cmesh = t8_cmesh_new_hypercube (T8_ECLASS_QUAD, sc_MPI_COMM_WORLD, 0, 0, 1);
    default_scheme = t8_scheme_new_default_cxx ();

    forest = t8_forest_new_uniform (cmesh, default_scheme, 0, 1, sc_MPI_COMM_WORLD);

    coarse_geometry = t8dg_coarse_geometry_new_2D_linear ();

    values = t8dg_values_new_LGL_hypercube (2, 2, coarse_geometry, forest);

  }

  /* *INDENT-OFF* */
  void TearDown () override
  /* *INDENT-ON* */
  {
    t8dg_values_destroy (&values);
    t8dg_coarse_geometry_destroy (&coarse_geometry);
    t8_forest_unref (&forest);
  }

  t8dg_values_t      *values;
  t8dg_coarse_geometry_t *coarse_geometry;
  t8_forest_t         forest;
};

typedef struct t8dg_test_values_problem
{
  t8dg_values_t      *values;
  t8dg_dof_values_t  *dof_values;
  t8dg_dof_values_t  *dof_values_adapt;
} t8dg_test_values_problem_t;

static void
t8dg_test_values_replace_all (t8_forest_t forest_old,
                              t8_forest_t forest_new,
                              t8_locidx_t itree,
                              t8_eclass_scheme_c * ts,
                              int num_elems_old, t8_locidx_t first_ielem_old, int num_elems_new, t8_locidx_t first_ielem_new)
{
  t8dg_test_values_problem_t *problem;
  int                 ichild, num_children;

  problem = (t8dg_test_values_problem_t *) t8_forest_get_user_data (forest_new);

  t8_element_t       *element = t8_forest_get_element_in_tree (forest_old, itree, first_ielem_old);
  num_children = ts->t8_element_num_children (element);

  T8DG_ASSERT (num_children == num_elems_new && num_elems_old == 1);

  for (ichild = 0; ichild < num_children; ichild++) {
    t8dg_values_set_element_adapt (problem->values, itree, first_ielem_new + ichild);
    t8dg_values_transform_parent_dof_to_child_dof (problem->values, problem->dof_values, problem->dof_values_adapt, itree,
                                                   first_ielem_old, first_ielem_new + ichild, ichild);
  }
}

int
t8dg_test_values_refine_all (t8_forest_t forest,
                             t8_forest_t forest_from,
                             t8_locidx_t itree,
                             t8_locidx_t lelement_id, t8_eclass_scheme_c * ts, int num_elements, t8_element_t * elements[])
{
  return 1;
}

TEST_F (ValuesChildInterpolationTest2D, lgl2_element)
{
  t8dg_global_values_t *global_values = t8dg_values_get_global_values_array (values)[T8_ECLASS_QUAD];

  ASSERT_EQ (t8dg_global_values_get_num_faces (global_values), 4);

  t8dg_dof_values_t  *dof_values_adapt;
  t8dg_dof_values_t  *dof_children_expected;
  t8dg_dof_values_t  *dof_values;

  double              parent_dof[4] = { 1, 2, 3, 5 };
  double              child_dof[16] = { 1, 1.5, 2, 2.75,
    1.5, 2, 2.75, 3.5,
    2, 2.75, 3, 4,
    2.75, 3.5, 4, 5
  };

  dof_values = t8dg_dof_values_new_data_local (forest, t8dg_values_get_global_values_array (values), parent_dof, 4);

  t8_forest_t         forest_adapt;

  t8dg_test_values_problem_t problem = { values, dof_values, NULL };

  t8_forest_init (&forest_adapt);
  t8_forest_set_user_data (forest_adapt, &problem);
  t8_forest_set_adapt (forest_adapt, forest, t8dg_test_values_refine_all, 0);
  t8_forest_commit (forest_adapt);

  t8dg_values_allocate_adapt (values, forest_adapt);

  dof_values_adapt = t8dg_dof_values_new (forest_adapt, t8dg_values_get_global_values_array (values));
  problem.dof_values_adapt = dof_values_adapt;
  t8_forest_iterate_replace (forest_adapt, forest, t8dg_test_values_replace_all);

  t8dg_values_cleanup_adapt (values);

  forest = forest_adapt;
  forest_adapt = NULL;

  dof_children_expected = t8dg_dof_values_new_data_local (forest, t8dg_values_get_global_values_array (values), child_dof, 16);
  EXPECT_TRUE (t8dg_dof_values_equal (dof_values_adapt, dof_children_expected));

  t8dg_dof_values_destroy (&dof_values);
  t8dg_dof_values_destroy (&dof_values_adapt);
  t8dg_dof_values_destroy (&dof_children_expected);
  EXPECT_TRUE (1);
}

/*TODO: do a test for Face child interpolation*/
/*TOOD: do a test for projection*/


/* *INDENT-OFF* */
class PrecomputedValuesProjectionTest1D : public ::testing::TestWithParam<int>
/* *INDENT-ON* */
{
protected:
  /* *INDENT-OFF* */
  void SetUp () override
  /* *INDENT-ON* */
  {
    t8_cmesh_t          cmesh;
    t8_scheme_cxx_t    *default_scheme;

    cmesh = t8_cmesh_new_hypercube (T8_ECLASS_LINE, sc_MPI_COMM_WORLD, 0, 0, 1);
    default_scheme = t8_scheme_new_default_cxx ();

    forest = t8_forest_new_uniform (cmesh, default_scheme, 0, 1, sc_MPI_COMM_WORLD);
    coarse_geometry = t8dg_coarse_geometry_new_1D_linear ();

    int                 num_lgl = GetParam ();
    values = t8dg_values_new_LGL_hypercube (1, num_lgl, coarse_geometry, forest);

    dof_values = t8dg_dof_values_new (forest, t8dg_values_get_global_values_array (values));
  }

  /* *INDENT-OFF* */
  void TearDown () override
  /* *INDENT-ON* */
  {
    t8dg_values_destroy (&values);
    t8dg_coarse_geometry_destroy (&coarse_geometry);
    t8_forest_unref (&forest);
    t8dg_dof_values_destroy (&dof_values);
  }
  t8dg_values_t      *values;
  t8dg_coarse_geometry_t *coarse_geometry;
  t8_forest_t         forest, forest_adapt;
  t8dg_dof_values_t  *dof_values;
};

INSTANTIATE_TEST_SUITE_P (lglRange, PrecomputedValuesProjectionTest1D, testing::Range (2, MAX_LGL_NUMBER + 1));

TEST_P (PrecomputedValuesProjectionTest1D, const_one)
{
  int                 num_lgl, idof, ichild;
  sc_array_t         *element_dof_values_child[2];
  sc_array_t         *element_dof_values_parent;

  t8_forest_t         forest_adapt;
  t8dg_dof_values_t  *dof_values_adapt;

  t8dg_test_values_problem_t problem = { values, dof_values, NULL };

  t8_forest_init (&forest_adapt);
  t8_forest_set_user_data (forest_adapt, &problem);
  t8_forest_set_adapt (forest_adapt, forest, t8dg_test_values_refine_all, 0);
  t8_forest_commit (forest_adapt);

  t8dg_values_allocate_adapt (values, forest_adapt);

  dof_values_adapt = t8dg_dof_values_new (forest_adapt, t8dg_values_get_global_values_array (values));
  problem.dof_values_adapt = dof_values_adapt;
  t8_forest_iterate_replace (forest_adapt, forest, t8dg_test_values_replace_all);

  t8dg_values_cleanup_adapt (values);

  forest = forest_adapt;
  forest_adapt = NULL;

  num_lgl = GetParam ();
#if 0
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
#endif
  t8dg_dof_values_destroy (&dof_values_adapt);
  EXPECT_TRUE (1);
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
    int                 level = 4;      //dependent on parameter?
    int                 itree, ielement;
    t8_cmesh_t          cmesh;
    t8_scheme_cxx_t    *default_scheme;

    cmesh = t8_cmesh_new_periodic_line_more_trees (sc_MPI_COMM_WORLD);
    default_scheme = t8_scheme_new_default_cxx ();
    forest = t8_forest_new_uniform (cmesh, default_scheme, level, 1, sc_MPI_COMM_WORLD);
    coarse_geometry = t8dg_coarse_geometry_new_1D_linear ();

    num_lgl = std::get < 0 > (GetParam ());

    values = t8dg_values_new_LGL_hypercube (1, num_lgl, coarse_geometry, forest);
    dof_values = t8dg_dof_values_new (forest, t8dg_values_get_global_values_array (values));
  }

  /* *INDENT-OFF* */
  void TearDown () override
  /* *INDENT-ON* */
  {
    t8dg_values_destroy (&values);
    t8dg_coarse_geometry_destroy (&coarse_geometry);
    t8dg_dof_values_destroy (&dof_values);
    t8_forest_unref (&forest);
  }

  t8dg_values_t      *values;
  t8dg_coarse_geometry_t *coarse_geometry;

  t8dg_dof_values_t  *dof_values;
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
  t8_locidx_t         itree = 0, ielement;
  double              time = 0;
  double              norm;
  t8dg_element_dof_values_t *element_dof_values;

  t8dg_values_interpolate_scalar_function_3d_time (values, std::get < 0 > (std::get < 1 > (GetParam ())), time, dof_values);

  norm = t8dg_values_norm_l2 (values, dof_values, sc_MPI_COMM_WORLD);
  EXPECT_NEAR (norm, std::get < 1 > (std::get < 1 > (GetParam ())), exp (-4 * std::get < 0 > (GetParam ())));
}
