#ifndef SRC_GLOBAL_H_
#define SRC_GLOBAL_H_

#include <sc_containers.h>
#include <sc_dmatrix.h>
#include <t8_forest.h>
#include <example/common/t8_example_common.h>

#include "t8dg.h"


struct t8dg_LGL_quadrature
{
  int 			number_of_quadrature_points;
  t8dg_LGL_vertex_set_t	*vertices;
  sc_array_t		*weights;
};

struct t8dg_LGL_functionbasis
{
  int 					number_of_dof;
  t8dg_matrix_application		directional_derivative_matrix;
  t8dg_LGL_vertex_set_t			*vertices;
};

struct t8dg_LGL_vertex_set{
  int			dim;
/*  int 			tensorflag;
  t8dg_quadrature_t 	*tensor1;
  t8dg_quadrature_t 	*tensor2;*/
  int			number_of_vertices;
  int			number_of_faces;
  int			number_of_facevertices[MAX_FACES];
  sc_array_t		*vertices; /* dim * number_of_vertices, make access available via function and allocate only if !tensor? */
  sc_array_t		*facevertex_indices[MAX_FACES];
};

void t8dg_LGL_transform_element_dof_to_face_quad(sc_array_t * face_quad_array, const sc_array_t * element_dof_array,t8dg_LGL_vertex_set_t *vertices);
void t8dg_LGL_transform_face_quad_to_element_dof(sc_array_t * element_dof_array, const sc_array_t * face_quad_array ,t8dg_LGL_vertex_set_t *vertices);

void t8dg_LGL_quadrature_and_functionbasis_new_1D(t8dg_LGL_quadrature_t **pquadrature,t8dg_LGL_functionbasis_t **pfunctionbasis,int number_of_LGL);
void t8dg_LGL_quadrature_and_functionbasis_destroy(t8dg_LGL_quadrature_t **pquadrature,t8dg_LGL_functionbasis_t **pfunctionbasis);



#endif /* SRC_GLOBAL_H_ */
