/*
 * timestepping.c
 *
 *  Created on: Mar 14, 2020
 *      Author: lukas
 */

#include<sc_containers.h>
#include "t8dg_sc_array.h"
#include "t8dg.h"
#include "t8dg_timestepping.h"

struct t8dg_timestepping_data
{
  int                 time_order;/**< time order of the Runge-kutta timestepping*/
  double              delta_t;  /**< time step */
  double              t;        /**< current time */
  double              T;        /**< end time */
  double              cfl;      /**< cfl number*/
  int                 step_number;
};

/*pre computed butcher tableau values for rk with values only on first diagonal*/
double              rk1_b[1] = { 1 };

double              rk2_a[2] = { 1 };
double              rk2_b[2] = { 0.5, 0.5 };

double              rk3_a[3] = { 1. / 3, 2. / 3 };
double              rk3_b[3] = { 1. / 4, 0, 3. / 4 };

double              rk4_a[4] = { 0.5, 0.5, 1 };
double              rk4_b[4] = { 1. / 6, 1. / 3, 1. / 3, 1. / 6 };

void
t8dg_runge_kutta_fill_coefficients (int time_order, double **prk_a, double **prk_b, double **prk_c)
{
  switch (time_order) {
  case 1:
    *prk_a = NULL;
    *prk_b = rk1_b;
    break;
  case 2:
    *prk_a = rk2_a;
    *prk_b = rk2_b;
    break;
  case 3:
    *prk_a = rk3_a;
    *prk_b = rk3_b;
    break;
  case 4:
    *prk_a = rk4_a;
    *prk_b = rk4_b;
    break;
  default:
    break;
  }
  *prk_c = *prk_a;
}

void
t8dg_timestepping_runge_kutta_step (t8dg_time_matrix_application time_derivative,
                                    t8dg_timestepping_data_t * time_data, sc_array_t ** pdof_array, void *user_data)
{
  int                 istep;
  double             *rk_a, *rk_b, *rk_c;
  sc_array_t         *element_dof_beginning;
  sc_array_t         *element_dof_change;
  sc_array_t         *element_dof_new;
  sc_array_t         *element_dof_step;
  double              time_beginning, time_current, time_step;
  int                 time_order = t8dg_timestepping_data_get_time_order (time_data);

  time_beginning = t8dg_timestepping_data_get_current_time (time_data);
  time_current = time_beginning;
  time_step = t8dg_timestepping_data_get_time_step (time_data);

  t8dg_runge_kutta_fill_coefficients (time_order, &rk_a, &rk_b, &rk_c);

  element_dof_beginning = t8dg_sc_array_clone (*pdof_array);
  element_dof_new = t8dg_sc_array_duplicate (element_dof_beginning);
  element_dof_step = t8dg_sc_array_duplicate (element_dof_beginning);
  element_dof_change = t8dg_sc_array_duplicate (element_dof_beginning);

  time_derivative (*pdof_array, element_dof_change, time_current, user_data);
  t8dg_sc_array_block_double_axpyz (rk_b[0] * time_step, element_dof_change, element_dof_beginning, element_dof_new);

  for (istep = 0; istep < time_order - 1; istep++) {
    /*calculate the y-value for which the derivative needs to be evaluated
     * since a has only values on the first minor diagonal, only the k from the step before and the original y is needed*/

    t8dg_sc_array_block_double_axpyz (rk_a[istep] * time_step, element_dof_change, element_dof_beginning, element_dof_step);
    t8dg_sc_array_swap (pdof_array, &element_dof_step);
    /* calculate the derivative at the step time and y value */

    time_current = time_beginning + rk_c[istep] * time_step;
    t8dg_timestepping_data_set_current_time (time_data, time_current);
    time_derivative (*pdof_array, element_dof_change, time_current, user_data);
    /*add weighted summand to result */
    t8dg_sc_array_block_double_axpy (rk_b[istep + 1] * time_step, element_dof_change, element_dof_new);
  }

  t8dg_timestepping_data_set_current_time (time_data, time_beginning + time_step);
  t8dg_sc_array_swap (pdof_array, &element_dof_new);

  sc_array_destroy (element_dof_beginning);
  sc_array_destroy (element_dof_change);
  sc_array_destroy (element_dof_new);
  sc_array_destroy (element_dof_step);

}

t8dg_timestepping_data_t *
t8dg_timestepping_data_new (int time_order, double start_time, double end_time, double cfl)
{
  T8DG_ASSERT (time_order > 0);
  t8dg_timestepping_data_t *time_data = T8DG_ALLOC (t8dg_timestepping_data_t, 1);
  time_data->time_order = time_order;
  time_data->t = start_time;
  time_data->T = end_time;
  time_data->cfl = cfl;
  time_data->step_number = 0;
  time_data->delta_t = -1;
  return time_data;
}

void
t8dg_timestepping_data_destroy (t8dg_timestepping_data_t ** ptime_data)
{
  t8dg_timestepping_data_t *time_data = *ptime_data;
  time_data->time_order = -1;
  time_data->t = -1;
  time_data->T = -1;
  time_data->cfl = -1;
  time_data->delta_t = -1;
  T8_FREE (time_data);
  *ptime_data = NULL;
}

double
t8dg_timestepping_data_get_current_time (t8dg_timestepping_data_t * time_data)
{
  return time_data->t;
}

void
t8dg_timestepping_data_set_current_time (t8dg_timestepping_data_t * time_data, double current_time)
{
  time_data->t = current_time;
}

double
t8dg_timestepping_data_get_end_time (t8dg_timestepping_data_t * time_data)
{
  return time_data->T;
}

double
t8dg_timestepping_data_get_cfl (t8dg_timestepping_data_t * time_data)
{
  return time_data->cfl;
}

int
t8dg_timestepping_data_get_time_order (t8dg_timestepping_data_t * time_data)
{
  return time_data->time_order;
}

double
t8dg_timestepping_data_get_time_step (t8dg_timestepping_data_t * time_data)
{
  return time_data->delta_t;
}

void
t8dg_timestepping_data_set_time_step (t8dg_timestepping_data_t * time_data, double delta_t)
{
  if (delta_t < (time_data->T - time_data->t)) {
    time_data->delta_t = delta_t;
  }
  else {
    time_data->delta_t = time_data->T - time_data->t;
  }
}

int
t8dg_timestepping_data_is_endtime_reached (t8dg_timestepping_data_t * time_data)
{
  return !(time_data->t < time_data->T);
}

int
t8dg_timestepping_data_get_step_number (t8dg_timestepping_data_t * time_data)
{
  return time_data->step_number;
}

void
t8dg_timestepping_data_increase_step_number (t8dg_timestepping_data_t * time_data)
{
  time_data->step_number++;
}

double
t8dg_timestepping_data_get_time_left (t8dg_timestepping_data_t * time_data)
{
  return time_data->T - time_data->t;
}
