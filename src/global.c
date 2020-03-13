/*
 * global.c
 *
 *  Created on: Mar 13, 2020
 *      Author: lukas
 */
#include "global.h"
#include <t8.h>
#include <sc_containers.h>

void identity_matrix (sc_array_t *dest, const sc_array_t *src, const void *application_data){
  T8_ASSERT (dest->elem_size == src->elem_size);
  T8_ASSERT (dest->elem_count == src->elem_count);
  /* Copy the elem_count*elem_size many bits from src to dest*/
  memcpy (dest->array, src->array, src->elem_count * src->elem_size);
}

/*application_data is faceindex integer*/
void face_vandermonde_1D_linear_LGL (sc_array_t *dest, const sc_array_t *src, const void *application_data){
  T8_ASSERT (dest->elem_size == src->elem_size);
  T8_ASSERT (dest->elem_count == 1);
  T8_ASSERT (src->elem_count == 2);
  int faceindex = *((int*)  application_data);
  T8_ASSERT (faceindex >=0 && faceindex <=1);
  *((double *)t8_sc_array_index_locidx(dest,0))= *((double *)t8_sc_array_index_locidx(src,faceindex));
  dest->array[0] = src->array [faceindex];
}

/*application_data: faceindex integer*/
void face_vandermonde_transpose_1D_linear_LGL (sc_array_t *dest, const sc_array_t *src, const void *application_data){
  T8_ASSERT (dest->elem_size == src->elem_size);
  T8_ASSERT (dest->elem_count == 2);
  T8_ASSERT (src->elem_count == 1);
  int faceindex = *((int*)  application_data);
  T8_ASSERT (faceindex >=0 && faceindex <=1);

  /* TODO: Falsch!!!  */
  dest->array[faceindex] = src->array [0];
  dest->array[1 - faceindex] = 0;
}

void directional_derivative_1D_LGL2_matrix(sc_array_t *dest, const sc_array_t *src, const void *application_data){
  T8_ASSERT (dest->elem_size == src->elem_size);
  T8_ASSERT (dest->elem_count == 2);
  T8_ASSERT (src->elem_count == 2);

  /* TODO: Falsch!!!  */
  dest->array[0]=src->array[0]+src->array[1];
  dest->array[1]=-(src->array[0]+src->array[1]);
}

/*application data is flow_velocity*/
double upwind_flux_1D(const double u_minus,const double u_plus, const void *application_data){
  T8_ASSERT(application_data!=NULL);
  double *c = (double*) application_data;
  if(*c > 0)return u_minus;
  return u_plus;
}

t8dg_quadrature_t * t8dg_1D_LGL_quadrature(int number_of_LGL)
{
  T8_ASSERT(number_of_LGL >=2);
  t8dg_quadrature_t *rquad = T8_ALLOC(t8dg_quadrature_t,1);
  rquad->tensorflag =0;
  rquad->tensor1=NULL;
  rquad->tensor2=NULL;
  rquad->number_of_vertices = number_of_LGL;
  rquad->number_of_faces = 2;
  rquad->number_of_facevertices[0]=1;
  rquad->number_of_facevertices[1]=1;/*TODO: Fill possible other values*/
  rquad->vertices = sc_array_new_count(sizeof(double),rquad->number_of_vertices);
  rquad->weights = sc_array_new_count(sizeof(double),rquad->number_of_vertices);

  T8_ASSERT(number_of_LGL == 2);
  *((double *)t8_sc_array_index_locidx(rquad->vertices,0))=0;
  *((double *)t8_sc_array_index_locidx(rquad->vertices,1))=1;
  *((double *)t8_sc_array_index_locidx(rquad->weights,0))=0.5;
  *((double *)t8_sc_array_index_locidx(rquad->weights,1))=0.5;

  return rquad;
}
void t8_dg_quadrature_destroy(t8dg_quadrature_t **pquadrature){
  sc_array_destroy((*pquadrature)->vertices);
  sc_array_destroy((*pquadrature)->weights);
  T8_FREE(*pquadrature);
  *pquadrature = NULL;
}
t8dg_functionbasis_t * t8dg_1D_LGL_functionbasis(int number_of_LGL){

}
void t8dg_functionbasis_destroy(t8dg_functionbasis_t **pfunctionbasis);


