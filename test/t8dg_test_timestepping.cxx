#include <gtest/gtest.h>
#include <sc_containers.h>
#include "../src/t8dg.h"
#include "../src/t8dg_timestepping.h"

void
square_derivative (const sc_array_t * src, sc_array_t * dest, const double t, const void *application_data)
{
  *(double *) t8dg_sc_array_index_locidx (dest, 0) = *(double *) t8dg_sc_array_index_locidx (src, 0) * 2 / t;
}

TEST (timestepping, ode_x0_times_tsquared)
{
  int                 order;
  int                 N = 1000;
  int                 step;
  double              tolerance = 0.02;
  t8dg_timestepping_data_t *time_data;
  sc_array_t         *y;
  y = sc_array_new_count (sizeof (double), 1);

  for (order = 1; order <= 4; order++) {
    *(double *) sc_array_index (y, 0) = 1;
    time_data = t8dg_timestepping_data_new (order, 1, 2, -1);
    t8dg_timestepping_data_set_time_step (time_data, 1. / N);
    for (step = 0; step < N; step++) {
      t8dg_timestepping_runge_kutta_step (square_derivative, time_data, &y, NULL);
    }
    EXPECT_NEAR (*(double *) sc_array_index (y, 0), 4, tolerance);
    tolerance /= N;
    t8dg_timestepping_data_destroy (&time_data);
  }
  sc_array_destroy (y);
}
