/*
 * global.c
 *
 *  Created on: Mar 13, 2020
 *      Author: lukas
 */
#include "t8dg.h"
#include "t8dg_LGL.h"
#include <sc_containers.h>

/*those can be generalized with vertex set!*/

/*application_data is faceindex integer*/
void t8dg_face_vandermonde_1D_linear_LGL (sc_array_t *dest, const sc_array_t *src, const void *application_data){
  T8DG_ASSERT (dest->elem_size == src->elem_size);
  T8DG_ASSERT (dest->elem_count == 1);
  T8DG_ASSERT (src->elem_count == 2);
  int faceindex = *((int*)  application_data);
  T8DG_ASSERT (faceindex >=0 && faceindex <=1);
  double *double_dest = (double*) dest->array;
  const double *double_src  = (double*) src->array;
  double_dest[0]=double_src[faceindex];
}

/*application_data: faceindex integer*/
void t8dg_face_vandermonde_transpose_1D_linear_LGL (sc_array_t *dest, const sc_array_t *src, const void *application_data){
  T8DG_ASSERT (dest->elem_size == src->elem_size);
  T8DG_ASSERT (dest->elem_count == 2);
  T8DG_ASSERT (src->elem_count == 1);
  int faceindex = *((int*)  application_data);
  T8DG_ASSERT (faceindex >=0 && faceindex <=1);

  double *double_dest = (double*) dest->array;
  const double *double_src = (double*) src->array;
  double_dest[faceindex] = double_src[0];
  double_dest[1-faceindex] = 0;
}

/* f_0' =1/2*(f_0+f_1) for basisfunctions f_0(x) = x, f_1(x) = 1-x */
void t8dg_directional_derivative_1D_LGL2_matrix(sc_array_t *dest, const sc_array_t *src, const void *application_data){
  T8DG_ASSERT (dest->elem_size == src->elem_size);
  T8DG_ASSERT (dest->elem_count == 2);
  T8DG_ASSERT (src->elem_count == 2);

  double *double_dest = (double*) dest->array;
  const double *double_src = (double*) src->array;
  double_dest[0] = double_src[0]+double_src[1];
  double_dest[1] = -double_dest[0];
}

static t8dg_LGL_vertex_set_t *
t8dg_1D_LGL_vertex_set(int number_of_LGL_vertices){
  T8DG_ASSERT(number_of_LGL_vertices >=1 && number_of_LGL_vertices<=4);
  t8dg_LGL_vertex_set_t *vertices = T8DG_ALLOC(t8dg_LGL_vertex_set_t, 1);
  vertices->dim = 1;
  return vertices;
}
static t8dg_LGL_quadrature_t *
t8dg_1D_LGL_quadrature(t8dg_LGL_vertex_set_t *vertex_set)
{
  T8DG_ASSERT(vertex_set->dim == 1);
  t8dg_LGL_quadrature_t *rquad = T8DG_ALLOC(t8dg_LGL_quadrature_t,1);
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
      T8DG_ASSERT(0);
  }
  return rquad;
}
static t8dg_LGL_functionbasis_t *
t8dg_1D_LGL_functionbasis(t8dg_LGL_vertex_set_t *vertex_set){
  t8dg_LGL_functionbasis_t *rfunctionbasis = T8DG_ALLOC(t8dg_LGL_functionbasis_t,1);
  rfunctionbasis->vertices = vertex_set;
  rfunctionbasis->directional_derivative_matrix = t8dg_directional_derivative_1D_LGL2_matrix;
  return rfunctionbasis;
}


void t8dg_1D_LGL_new_quadrature_and_functionbasis(t8dg_LGL_quadrature_t **pquadrature,t8dg_LGL_functionbasis_t **pfunctionbasis,int number_of_LGL_vertices){
  t8dg_LGL_vertex_set_t *vertices;
  t8dg_LGL_quadrature_t *quadrature;
  t8dg_LGL_functionbasis_t *functionbasis;
  vertices = t8dg_1D_LGL_vertex_set(number_of_LGL_vertices);
  quadrature = t8dg_1D_LGL_quadrature(vertices);
  functionbasis = t8dg_1D_LGL_functionbasis(vertices);
  *pquadrature = quadrature;
  *pfunctionbasis = functionbasis;
}

void t8dg_LGL_quadrature_and_functionbasis_destroy(t8dg_LGL_quadrature_t **pquadrature,t8dg_LGL_functionbasis_t **pfunctionbasis){

}


