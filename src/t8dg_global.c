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

t8dg_quadrature_t * t8dg_1D_LGL_quadrature(int number_of_LGL)
{
  T8_ASSERT(number_of_LGL ==2);
  t8dg_quadrature_t *rquad = T8_ALLOC(t8dg_quadrature_t,1);
  rquad->tensorflag =0;
  rquad->tensor1=NULL;
  rquad->tensor2=NULL;
  rquad->number_of_vertices = number_of_LGL;
  rquad->number_of_faces = 2;
  rquad->number_of_facevertices[0]=1;
  rquad->number_of_facevertices[1]=1;
  rquad->vertices = sc_array_new_count(sizeof(double),rquad->number_of_vertices);
  rquad->weights = sc_array_new_count(sizeof(double),rquad->number_of_vertices);

  double *vertices, *weights;
  vertices = (double *)rquad->vertices->array;
  weights = (double *)rquad->weights->array;

  switch(number_of_LGL){
    case(2):
	vertices[0] = 0;
	vertices[1] = 1;
	weights[0] = 0.5;
	weights[1] = 0.5;
	break;
  }
  return rquad;
}
void t8dg_quadrature_destroy(t8dg_quadrature_t **pquadrature){
  sc_array_destroy((*pquadrature)->vertices);
  sc_array_destroy((*pquadrature)->weights);
  T8_FREE(*pquadrature);
  *pquadrature = NULL;
}

t8dg_functionbasis_t * t8dg_1D_LGL_functionbasis(int number_of_LGL){
  T8_ASSERT(number_of_LGL==2);
  t8dg_functionbasis_t *functionbasis = T8_ALLOC(t8dg_functionbasis_t,1);
  functionbasis->number_of_dof = number_of_LGL;
  switch(number_of_LGL){
    case(2):
	functionbasis->directional_derivative_matrix = t8dg_directional_derivative_1D_LGL2_matrix;
	break;
  }
  return functionbasis;
}
void t8dg_functionbasis_destroy(t8dg_functionbasis_t **pfunctionbasis){
  T8_FREE(*pfunctionbasis);
  *pfunctionbasis = NULL;
}


