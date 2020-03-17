/*
 * timestepping.c
 *
 *  Created on: Mar 14, 2020
 *      Author: lukas
 */

#include<sc_containers.h>
#include <t8.h>
#include "t8dg_sc_array.h"
#include "t8dg_global.h"


/*pre computed butcher tableau values for rk with values only on first diagonal*/
static const double rk1_b[1] = {1};

static const double rk2_a[2] = {1};
static const double rk2_b[2] = {0.5,0.5};

static const double rk3_a[3] = {1./3,2./3};
static const double rk3_b[3] = {1./4,0,3./4};

static const double rk4_a[4] = {0.5,0.5,1};
static const double rk4_b[4] = {1./6,1./3,1./3,1./6};

void t8dg_rungekutta_timestep(int order,const double t,const double delta_t,const t8dg_time_matrix_application f_matrix ,
			 sc_array_t *dest, sc_array_t *src, const void *application_data){
  T8_ASSERT(order>0&&order<=4);
  T8_ASSERT(f_matrix!=NULL);
  T8_ASSERT(dest!=NULL&&src!=NULL);
  T8_ASSERT(dest->elem_count==src->elem_count);
  T8_ASSERT(dest->elem_size==src->elem_size);

  const double *rk_a,*rk_b,*rk_c;
  size_t count=src->elem_count;
  size_t size = src->elem_size;
  int i;

  switch(order){
    case 1:
      rk_a=NULL;
      rk_b=rk1_b;
      break;
    case 2:
      rk_a=rk2_a;
      rk_b=rk2_b;
      break;
    case 3:
      rk_a=rk3_a;
      rk_b=rk3_b;
      break;
    case 4:
      rk_a=rk4_a;
      rk_b=rk4_b;
      break;
    default:
      break;
  }
  rk_c = rk_a;

  sc_array_t *y_start,*y_step,*y_res,*k_step;
  y_start = src;
  y_res = dest;
  k_step = sc_array_new_count(size,count);
  y_step = sc_array_new_count(size,count);


  /* k_0 = f(y_n,t)*/
  f_matrix(k_step,y_start,t,application_data);

  /*y_res has partial sums by weighting the intermediate values k*/
  t8dg_sc_array_block_double_zaxpy(y_res,rk_b[0]*delta_t,k_step,y_start);

  for(i=0;i<order-1;i++){
    /*calculate the y-value for which the derivative needs to be evaluated
     * since a has only values on the first minor diagonal, only the k from the step before and the original y is needed*/
    t8dg_sc_array_block_double_zaxpy(y_step,rk_a[i]*delta_t,k_step,y_start);
    /* calculate the derivative at the step time and y value */
    f_matrix(k_step, y_step, t + rk_c[i] *delta_t , application_data);
      /*add weighted summand to result*/
    t8dg_sc_array_block_double_axpy(rk_b[i+1]*delta_t,k_step,y_res);
  }

  sc_array_destroy(k_step);
  sc_array_destroy(y_step);
}
