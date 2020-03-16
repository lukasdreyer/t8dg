/*
 * t8dg_sc_array.c
 *
 *  Created on: Mar 14, 2020
 *      Author: lukas
 */
#include<t8.h>

void t8dg_sc_array_block_double_axpy(const double a, const sc_array_t *x, sc_array_t *y){
  T8_ASSERT(x!=NULL&&y!=NULL);
  T8_ASSERT(x->array!=NULL&&y->array!=NULL);
  T8_ASSERT(x->elem_size%sizeof(double)==0&&y->elem_size%sizeof(double)==0);
  T8_ASSERT(x->elem_count*x->elem_size == y->elem_count * y->elem_size);
  double *x_double, *y_double;
  int double_count,i;

  /*View array as double array*/
  x_double = (double*) x->array;
  y_double = (double*) y->array;
  double_count = x->elem_size/sizeof(double)*x->elem_count;/*total number of doubles*/

  for(i=0;i<double_count;i++){
      y_double[i] = a*x_double[i] + y_double[i];
  }
}

void t8dg_sc_array_block_double_zaxpy(sc_array_t *z, double a, const sc_array_t *x, const sc_array_t *y){
  T8_ASSERT(z!=NULL&&x!=NULL&&y!=NULL);
  T8_ASSERT(z->array!=NULL&&x->array!=NULL&&y->array!=NULL);
  T8_ASSERT(x->elem_size%sizeof(double)==0&&y->elem_size%sizeof(double)==0&&z->elem_size%sizeof(double)==0);
  T8_ASSERT(x->elem_count*x->elem_size == y->elem_count * y->elem_size);
  T8_ASSERT(x->elem_count*x->elem_size == z->elem_count * z->elem_size);

  /*View array as double array*/
  double *x_double, *y_double, *z_double;
  int double_count,i;

  x_double = (double*) (x->array);
  y_double = (double*) (y->array);
  z_double = (double*) (z->array);

  double_count = (x->elem_size/sizeof(double))*x->elem_count; /*total number of doubles*/

  for(i=0;i<double_count;i++){
      z_double[i] = a*x_double[i]+y_double[i];
  }
}

void t8dg_sc_array_swap(sc_array_t ** parray1,sc_array_t ** parray2){
  sc_array_t *temp;
  temp = *parray1;
  *parray1 = *parray2;
  *parray2 = temp;
}
