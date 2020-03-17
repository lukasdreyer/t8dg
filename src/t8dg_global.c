/*
 * global.c
 *
 *  Created on: Mar 13, 2020
 *      Author: lukas
 */
#include "t8dg_global.h"
#include <t8.h>
#include <sc_containers.h>

void t8dg_identity_matrix (sc_array_t *dest, const sc_array_t *src, const void *application_data){
  T8_ASSERT (dest->elem_size == src->elem_size);
  T8_ASSERT (dest->elem_count == src->elem_count);
  /* Copy the elem_count*elem_size many bits from src to dest*/
  memcpy (dest->array, src->array, src->elem_count * src->elem_size);
}

/*application_data is faceindex integer*/
void t8dg_face_vandermonde_1D_linear_LGL (sc_array_t *dest, const sc_array_t *src, const void *application_data){
  T8_ASSERT (dest->elem_size == src->elem_size);
  T8_ASSERT (dest->elem_count == 1);
  T8_ASSERT (src->elem_count == 2);
  int faceindex = *((int*)  application_data);
  T8_ASSERT (faceindex >=0 && faceindex <=1);
  double *double_dest = (double*) dest->array;
  const double *double_src  = (double*) src->array;
  double_dest[0]=double_src[faceindex];
}

/*application_data: faceindex integer*/
void t8dg_face_vandermonde_transpose_1D_linear_LGL (sc_array_t *dest, const sc_array_t *src, const void *application_data){
  T8_ASSERT (dest->elem_size == src->elem_size);
  T8_ASSERT (dest->elem_count == 2);
  T8_ASSERT (src->elem_count == 1);
  int faceindex = *((int*)  application_data);
  T8_ASSERT (faceindex >=0 && faceindex <=1);

  double *double_dest = (double*) dest->array;
  const double *double_src = (double*) src->array;
  double_dest[faceindex] = double_src[0];
  double_dest[1-faceindex] = 0;
}

/* f_0' =1/2*(f_0+f_1) for basisfunctions f_0(x) = x, f_1(x) = 1-x */
void t8dg_directional_derivative_1D_LGL2_matrix(sc_array_t *dest, const sc_array_t *src, const void *application_data){
  T8_ASSERT (dest->elem_size == src->elem_size);
  T8_ASSERT (dest->elem_count == 2);
  T8_ASSERT (src->elem_count == 2);

  double *double_dest = (double*) dest->array;
  const double *double_src = (double*) src->array;
  double_dest[0] = double_src[0]+double_src[1];
  double_dest[1] = -double_dest[0];
}

/*TODO: Include normal vector or guarantee that u_minus is always left and u_plus always right.
 * Better: u_minus own value, u_plus neighbour value, normal vector and flow determine which value to chose!
 * */
/*application data is flow_velocity*/
double t8dg_upwind_flux_1D(const double u_minus,const double u_plus, const void *application_data){
  T8_ASSERT(application_data!=NULL);
  double *c = (double*) application_data;
  if(*c > 0)return u_minus;
  return u_plus;
}

static t8dg_LGL_vertex_set_t *
t8dg_1D_LGL_vertex_set(int number_of_LGL_vertices){
  T8_ASSERT(number_of_LGL_vertices >=1 && number_of_LGL_vertices<=4);
  t8dg_LGL_vertex_set_t *vertices = T8_ALLOC(t8dg_LGL_vertex_set_t, 1);
  vertices->dim = 1;


}
static t8dg_LGL_quadrature_t *
t8dg_1D_LGL_quadrature(t8dg_LGL_vertex_set_t *vertex_set)
{
  T8_ASSERT(vertex_set->dim == 1);
  t8dg_LGL_quadrature_t *rquad = T8_ALLOC(t8dg_LGL_quadrature_t,1);
  rquad->vertices = vertex_set;

  rquad->weights = sc_array_new_count(sizeof(double),vertex_set->number_of_vertices);

  double  *weights;
  weights = (double *)rquad->weights->array;

  switch(vertex_set->number_of_vertices){
    case(1):
	weights[0]= 1;
	break;
    case(2):
	weights[0] = 0.5;//Auf Referenzelement [0,1]
	weights[1] = 0.5;
	break;
    case(3):
	weights[0] = 1./3;
	weights[1] = 4./3;
	weights[2] = 1./3;
	break;
    case(4):
	weights[0] = 1./6;
	weights[1] = 5./6;
	weights[2] = 5./6;
	weights[3] = 1./6;
	break;
    default:
      printf("Not yet implemented!\n");
      T8_ASSERT(0);
  }
  return rquad;
}
void t8dg_quadrature_destroy(t8dg_LGL_quadrature_t **pquadrature){
  sc_array_destroy((*pquadrature)->vertices);
  sc_array_destroy((*pquadrature)->weights);
  T8_FREE(*pquadrature);
  *pquadrature = NULL;
}

t8dg_LGL_functionbasis_t * t8dg_1D_LGL_functionbasis(int number_of_LGL){
  T8_ASSERT(number_of_LGL==2);
  t8dg_LGL_functionbasis_t *functionbasis = T8_ALLOC(t8dg_LGL_functionbasis_t,1);
  functionbasis->number_of_dof = number_of_LGL;
  switch(number_of_LGL){
    case(2):
	functionbasis->directional_derivative_matrix = t8dg_directional_derivative_1D_LGL2_matrix;
	break;
  }
  return functionbasis;
}
void t8dg_functionbasis_destroy(t8dg_LGL_functionbasis_t **pfunctionbasis){
  T8_FREE(*pfunctionbasis);
  *pfunctionbasis = NULL;
}


