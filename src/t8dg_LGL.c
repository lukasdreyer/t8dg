/*
 * global.c
 *
 *  Created on: Mar 13, 2020
 *      Author: lukas
 */
#include "t8dg.h"
#include "t8dg_LGL.h"
#include <sc_containers.h>


void t8dg_LGL_vertex_set_get_3D_vertex(double reference_vertex[3], t8dg_LGL_vertex_set_t *vertex_set,int idof){
  double *vertex;
  int idim;
  vertex = (double*) t8_sc_array_index_locidx(vertex_set->vertices,idof);
  for(idim = 0 ; idim < vertex_set->dim; idim++){
      reference_vertex[idim] = vertex[idim];
  }
  for(idim = vertex_set->dim ; idim < DIM3 ; idim++){
      reference_vertex[idim] = 0;
  }
}


#if 0
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
#endif


double t8dg_LGL_quadrature_get_weight(t8dg_LGL_quadrature_t *quadrature, int iquad){
  return *(double *)t8_sc_array_index_locidx(quadrature->weights,iquad);
}


static t8dg_LGL_vertex_set_t *
t8dg_LGL_vertex_set_new_1D(int number_of_LGL_vertices){
  T8DG_ASSERT(number_of_LGL_vertices >=1 && number_of_LGL_vertices<=4);
  int iface;
  t8dg_LGL_vertex_set_t *vertices = T8DG_ALLOC(t8dg_LGL_vertex_set_t, 1);
  vertices->dim = 1;
  vertices->number_of_faces = 2;
  vertices->number_of_vertices = number_of_LGL_vertices;
  for(iface = 0; iface < vertices->number_of_faces; iface++){
      vertices->number_of_facevertices[iface]=1;
      vertices->facevertex_indices[iface] = sc_array_new_count(sizeof(int),vertices->number_of_facevertices[iface]);
  }
  vertices->number_of_facevertices[0]=1;
  vertices->vertices = sc_array_new_count(sizeof(double) * vertices->dim, vertices->number_of_vertices);

  double  *vertex_array;
  vertex_array = (double *) sc_array_index(vertices->vertices,0);

  switch(vertices->number_of_vertices){
    case(1):
	vertex_array[0]= 1;
	break;
    case(2):
	vertex_array[0] = 0;
	vertex_array[1] = 1;
	break;
    case(3):
	vertex_array[0] = 0;
	vertex_array[1] = 1./2;
	vertex_array[2] = 1;
	break;
    case(4):
	vertex_array[0] = 0;
	vertex_array[1] = (1-sqrt(5))/2;
	vertex_array[2] = (1+sqrt(5))/2;
	vertex_array[3] = 1;
	break;
    default:
      printf("Not yet implemented!\n");
      T8DG_ASSERT(0);
  }

  return vertices;
}
static t8dg_LGL_quadrature_t *
t8dg_LGL_quadrature_new(t8dg_LGL_vertex_set_t *vertex_set)
{
  T8DG_ASSERT(vertex_set->dim == 1);
  t8dg_LGL_quadrature_t *rquad = T8DG_ALLOC(t8dg_LGL_quadrature_t,1);
  rquad->vertices = vertex_set;
  rquad->number_of_quadrature_points = vertex_set->number_of_vertices;

  rquad->weights = sc_array_new_count(sizeof(double),vertex_set->number_of_vertices);

  double  *weights;
  weights = (double *)t8_sc_array_index_locidx(rquad->weights,0);

  switch(vertex_set->number_of_vertices){
    case(1):
	weights[0]= 1;
	break;
    case(2):
	weights[0] = 0.5;//Auf Referenzelement [0,1]
	weights[1] = 0.5;
	break;
    case(3):
	weights[0] = 1./6;
	weights[1] = 4./6;
	weights[2] = 1./6;
	break;
    case(4):
	weights[0] = 1./12;
	weights[1] = 5./12;
	weights[2] = 5./12;
	weights[3] = 1./12;
	break;
    default:
      printf("Not yet implemented!\n");
      T8DG_ASSERT(0);
  }
  return rquad;
}
static t8dg_LGL_functionbasis_t *
t8dg_LGL_functionbasis_new(t8dg_LGL_vertex_set_t *vertex_set){
  T8DG_ASSERT(vertex_set->dim == 1);
  t8dg_LGL_functionbasis_t *rfunctionbasis = T8DG_ALLOC(t8dg_LGL_functionbasis_t,1);
  rfunctionbasis->vertices = vertex_set;
  rfunctionbasis->number_of_dof = rfunctionbasis->vertices->number_of_vertices;

  /* TODO: check or implement generally! */
#if 0
  T8DG_ASSERT(vertex_set->number_of_vertices ==2);
  rfunctionbasis->directional_derivative_matrix = t8dg_directional_derivative_1D_LGL2_matrix;
#endif

  return rfunctionbasis;
}


void t8dg_LGL_quadrature_and_functionbasis_new_1D(t8dg_LGL_quadrature_t **pquadrature,t8dg_LGL_functionbasis_t **pfunctionbasis,int number_of_LGL_vertices){
  t8dg_LGL_vertex_set_t *vertices;
  t8dg_LGL_quadrature_t *quadrature;
  t8dg_LGL_functionbasis_t *functionbasis;
  vertices = t8dg_LGL_vertex_set_new_1D(number_of_LGL_vertices);
  quadrature = t8dg_LGL_quadrature_new(vertices);
  functionbasis = t8dg_LGL_functionbasis_new(vertices);
  *pquadrature = quadrature;
  *pfunctionbasis = functionbasis;
}

void t8dg_LGL_vertex_set_destroy(t8dg_LGL_vertex_set_t **pvertex_set){
  int iface = 0;
  t8dg_LGL_vertex_set_t *vertex_set = *pvertex_set;
  for(iface = 0; iface < vertex_set->number_of_faces; iface++){
      sc_array_destroy(vertex_set->facevertex_indices[iface]);
  }
  sc_array_destroy(vertex_set->vertices);
  T8DG_FREE(vertex_set);
}
void t8dg_LGL_quadrature_and_functionbasis_destroy(t8dg_LGL_quadrature_t **pquadrature,t8dg_LGL_functionbasis_t **pfunctionbasis){
  T8DG_ASSERT((*pquadrature)->vertices == (*pfunctionbasis)->vertices);
  t8dg_LGL_vertex_set_destroy(&(*pquadrature)->vertices);
  sc_array_destroy((*pquadrature)->weights);
  T8DG_FREE(*pquadrature);
  T8DG_FREE(*pfunctionbasis);
  *pquadrature = NULL;
  *pfunctionbasis = NULL;
}


