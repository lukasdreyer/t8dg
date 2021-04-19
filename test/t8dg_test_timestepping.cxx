#include <gtest/gtest.h>
#include <sc_containers.h>
#include "../src/t8dg.h"
#include "../src/t8dg_timestepping.h"
#include <t8_schemes/t8_default_cxx.hxx>
#include "../src/t8dg_global_values.h"

void
square_derivative (t8dg_dof_values_t * src_dof, t8dg_dof_values_t * dest_dof, const double t, const void *application_data)
{
  t8_forest_t         forest = (t8_forest_t) application_data;
  double              derivative_value;
  t8_locidx_t         itree, ielement, num_elems_in_tree;
  t8_locidx_t         num_trees = t8_forest_get_num_local_trees (forest);
  t8dg_element_dof_values_t *src_element_dof_view;
  t8dg_element_dof_values_t *dest_element_dof_view;
  t8dg_dofidx_t       idof;

  for (itree = 0; itree < num_trees; itree++) {
    num_elems_in_tree = t8_forest_get_tree_num_elements (forest, itree);
    for (ielement = 0; ielement < num_elems_in_tree; ielement++) {
      src_element_dof_view = t8dg_dof_values_new_element_dof_values_view (src_dof, itree, ielement);
      dest_element_dof_view = t8dg_dof_values_new_element_dof_values_view (dest_dof, itree, ielement);
      for (idof = 0; idof < t8dg_element_dof_values_get_num_dof (src_element_dof_view); idof++) {
        derivative_value = t8dg_element_dof_values_get_value (src_element_dof_view, idof) * 2 / t;
        t8dg_element_dof_values_set_value (dest_element_dof_view, idof, derivative_value);
      }
      t8dg_element_dof_values_destroy (&src_element_dof_view);
      t8dg_element_dof_values_destroy (&dest_element_dof_view);
    }
  }
}

TEST (timestepping, ode_x0_times_tsquared)
{
  int                 order;
  int                 N = 1000;
  int                 step;
  int                 uniform_level = 0;
  double              tolerance = 0.02;
  t8_cmesh_t          cmesh = t8_cmesh_new_hypercube (T8_ECLASS_LINE, sc_MPI_COMM_WORLD, 0, 0, 1);
  t8dg_timestepping_data_t *time_data;
  t8_scheme_cxx_t    *default_scheme;
  default_scheme = t8_scheme_new_default_cxx ();

  t8_forest_t         forest = t8_forest_new_uniform (cmesh, default_scheme, uniform_level, 1, sc_MPI_COMM_WORLD);
  t8dg_global_values_t *global_values_array[T8_ECLASS_COUNT] = { NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL };
  global_values_array[T8_ECLASS_LINE] = t8dg_global_values_new_hypercube_LGL (1, 2);
  t8dg_dof_values_t  *y = t8dg_dof_values_new (forest, global_values_array);

  t8_locidx_t         itree, ielement, num_elems_in_tree;
  t8_locidx_t         num_trees = t8_forest_get_num_local_trees (forest);
  t8dg_element_dof_values_t *element_dof_view;
  t8dg_dofidx_t       idof;

  for (order = 1; order <= 4; order++) {

    t8dg_dof_values_set_all_values (y, 1);

    for (itree = 0; itree < num_trees; itree++) {
      num_elems_in_tree = t8_forest_get_tree_num_elements (forest, itree);
      for (ielement = 0; ielement < num_elems_in_tree; ielement++) {
        element_dof_view = t8dg_dof_values_new_element_dof_values_view (y, itree, ielement);
        for (idof = 0; idof < t8dg_element_dof_values_get_num_dof (element_dof_view); idof++) {
          EXPECT_EQ (t8dg_element_dof_values_get_value (element_dof_view, idof), 1);
        }
        t8dg_element_dof_values_destroy (&element_dof_view);
      }
    }

    time_data = t8dg_timestepping_data_new_constant_timestep (order, 1, 2, 1. / N);

    for (step = 0; step < N; step++) {
      t8dg_timestepping_runge_kutta_step (square_derivative, time_data, &y, forest);
    }

    for (itree = 0; itree < num_trees; itree++) {
      num_elems_in_tree = t8_forest_get_tree_num_elements (forest, itree);
      for (ielement = 0; ielement < num_elems_in_tree; ielement++) {
        element_dof_view = t8dg_dof_values_new_element_dof_values_view (y, itree, ielement);
        for (idof = 0; idof < t8dg_element_dof_values_get_num_dof (element_dof_view); idof++) {
          EXPECT_NEAR (t8dg_element_dof_values_get_value (element_dof_view, idof), 4, tolerance);
        }
        t8dg_element_dof_values_destroy (&element_dof_view);
      }
    }
    tolerance /= N;
    t8dg_timestepping_data_destroy (&time_data);
  }
  t8dg_dof_values_destroy (&y);
  t8dg_global_values_destroy (&global_values_array[T8_ECLASS_LINE]);
  t8_forest_unref (&forest);
}
