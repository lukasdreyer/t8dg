/*
 * global.c
 *
 *  Created on: Mar 13, 2020
 *      Author: lukas
 */
#include "t8dg.h"
#include "t8dg_LGL.h"
#include <sc_containers.h>

/** Vertex set used for the quadrature and functionbasis
 * Currently only for 1D
 */
typedef struct t8dg_LGL_vertexset
{
  int                 dim;      /**< Dimension of the reference element*/
  int                 number_of_vertices;       /**< Number of element vertices*/
  int                 number_of_faces;          /**< Number of faces */
  int                 number_of_facevertices[MAX_FACES];        /**< For each face, the number of face vertices*/
  sc_array_t         *vertices; /**< size: dim * number_of_vertices, make access available via function and allocate only if not tensor? */
  sc_array_t         *facevertex_indices[MAX_FACES];    /**< Lookup table, for each facequadindex save the elementquadindex*/
#if 0
  int                 tensorflag;
  t8dg_quadrature_t  *tensor1;
  t8dg_quadrature_t  *tensor2;
                     */
#endif
} t8dg_LGL_vertexset_t;

/** Additionally to the LGL vertices save the quadrature weights*/
struct t8dg_LGL_quadrature
{
  int                 number_of_quadrature_points;      /**< Number of element quadrature points*/
  t8dg_LGL_vertexset_t *vertices;                       /**< LGL quadrature vertices*/
  sc_array_t         *weights;                          /**< LGL quadrature weights*/
};

/** The functionbasis provides the the directional derivative_matrix */
struct t8dg_LGL_functionbasis
{
  int                 number_of_dof;                    /**< Number of degrees of freedom*/
  t8dg_matrix_application directional_derivative_matrix;        /**< TODO: change to function that gets functionbasis and direction_idx as input*/
  t8dg_LGL_vertexset_t *vertices;                       /**< LGL vertices used as basis nodes for Lagrange nodal basis*/
};

static void
t8dg_vertexset_set_facevertex_index (t8dg_LGL_vertexset_t * vertices, int iface, t8dg_quad_idx_t ifacequad, t8dg_quad_idx_t ielemquad)
{
  *(t8dg_quad_idx_t *) t8dg_sc_array_index_quadidx (vertices->facevertex_indices[iface], ifacequad) = ielemquad;
}

static              t8dg_quad_idx_t
t8dg_LGL_vertexset_facequadidx_to_elementquadidx (t8dg_LGL_quadrature_t * quadrature, int iface, t8dg_quad_idx_t ifacequad)
{
  return *(t8dg_quad_idx_t *) t8dg_sc_array_index_quadidx (quadrature->vertices->facevertex_indices[iface], ifacequad);
}

static t8dg_LGL_vertexset_t *
t8dg_LGL_vertexset_new_1D (int number_of_LGL_vertices)
{
  T8DG_ASSERT (number_of_LGL_vertices >= 1 && number_of_LGL_vertices <= 4);     /*Larger not implemented */

  int                 iface;
  t8dg_LGL_vertexset_t *vertices;
  vertices = T8DG_ALLOC (t8dg_LGL_vertexset_t, 1);

  vertices->dim = 1;
  vertices->number_of_faces = 2;
  vertices->number_of_vertices = number_of_LGL_vertices;
  for (iface = 0; iface < vertices->number_of_faces; iface++) {
    vertices->number_of_facevertices[iface] = 1;
    vertices->facevertex_indices[iface] = sc_array_new_count (sizeof (int), vertices->number_of_facevertices[iface]);
  }
  t8dg_vertexset_set_facevertex_index (vertices, 0, 0, 0);      /* The facequad vertex of the left face is the 0th elementquad vertex */
  t8dg_vertexset_set_facevertex_index (vertices, 1, 0, number_of_LGL_vertices - 1);     /*The facequad vertex of the right face is the last elementquad vertex */

  vertices->vertices = sc_array_new_count (sizeof (double) * vertices->dim, vertices->number_of_vertices);

  double             *vertex_array;
  vertex_array = (double *) sc_array_index (vertices->vertices, 0);

  /*LGL vertices on [0,1] */
  switch (vertices->number_of_vertices) {
  case (1):
    vertex_array[0] = 1;
    break;
  case (2):
    vertex_array[0] = 0;
    vertex_array[1] = 1;
    break;
  case (3):
    vertex_array[0] = 0;
    vertex_array[1] = 1. / 2;
    vertex_array[2] = 1;
    break;
  case (4):
    vertex_array[0] = 0;
    vertex_array[1] = (1 - sqrt (5)) / 2;
    vertex_array[2] = (1 + sqrt (5)) / 2;
    vertex_array[3] = 1;
    break;
  default:
    printf ("Not yet implemented!\n");
    T8DG_ASSERT (0);
  }

  return vertices;
}

static t8dg_LGL_quadrature_t *
t8dg_LGL_quadrature_new (t8dg_LGL_vertexset_t * vertexset)
{
  T8DG_ASSERT (vertexset->dim == 1);    /*Other not yet implemented */

  t8dg_LGL_quadrature_t *rquad = T8DG_ALLOC (t8dg_LGL_quadrature_t, 1);

  rquad->vertices = vertexset;
  rquad->number_of_quadrature_points = vertexset->number_of_vertices;

  rquad->weights = sc_array_new_count (sizeof (double), vertexset->number_of_vertices);

  double             *weights;
  weights = (double *) t8_sc_array_index_locidx (rquad->weights, 0);

  /* sum of element weights = 1 = Vol([0,1]) */
  switch (vertexset->number_of_vertices) {
  case (1):
    weights[0] = 1;
    break;
  case (2):
    weights[0] = 0.5;
    weights[1] = 0.5;
    break;
  case (3):
    weights[0] = 1. / 6;
    weights[1] = 4. / 6;
    weights[2] = 1. / 6;
    break;
  case (4):
    weights[0] = 1. / 12;
    weights[1] = 5. / 12;
    weights[2] = 5. / 12;
    weights[3] = 1. / 12;
    break;
  default:
    printf ("Not yet implemented!\n");
    T8DG_ASSERT (0);
  }
  return rquad;
}

static t8dg_LGL_functionbasis_t *
t8dg_LGL_functionbasis_new (t8dg_LGL_vertexset_t * vertexset)
{
  t8dg_LGL_functionbasis_t *rfunctionbasis = T8DG_ALLOC (t8dg_LGL_functionbasis_t, 1);
  rfunctionbasis->vertices = vertexset;
  rfunctionbasis->number_of_dof = rfunctionbasis->vertices->number_of_vertices;
  /*TODO: directional derivative matrix */
  return rfunctionbasis;
}

void
t8dg_LGL_quadrature_and_functionbasis_new_1D (t8dg_LGL_quadrature_t ** pquadrature,
                                              t8dg_LGL_functionbasis_t ** pfunctionbasis, int number_of_LGL_vertices)
{
  t8dg_LGL_vertexset_t *vertices;
  t8dg_LGL_quadrature_t *quadrature;
  t8dg_LGL_functionbasis_t *functionbasis;
  vertices = t8dg_LGL_vertexset_new_1D (number_of_LGL_vertices);
  quadrature = t8dg_LGL_quadrature_new (vertices);
  functionbasis = t8dg_LGL_functionbasis_new (vertices);
  *pquadrature = quadrature;
  *pfunctionbasis = functionbasis;
}

void
t8dg_LGL_vertexset_destroy (t8dg_LGL_vertexset_t ** pvertexset)
{
  int                 iface = 0;
  t8dg_LGL_vertexset_t *vertexset = *pvertexset;
  for (iface = 0; iface < vertexset->number_of_faces; iface++) {
    sc_array_destroy (vertexset->facevertex_indices[iface]);
  }
  sc_array_destroy (vertexset->vertices);
  T8DG_FREE (vertexset);
}

void
t8dg_LGL_quadrature_and_functionbasis_destroy (t8dg_LGL_quadrature_t ** pquadrature, t8dg_LGL_functionbasis_t ** pfunctionbasis)
{
  T8DG_ASSERT ((*pquadrature)->vertices == (*pfunctionbasis)->vertices);
  t8dg_LGL_vertexset_destroy (&(*pquadrature)->vertices);
  sc_array_destroy ((*pquadrature)->weights);
  T8DG_FREE (*pquadrature);
  T8DG_FREE (*pfunctionbasis);
  *pquadrature = NULL;
  *pfunctionbasis = NULL;
}

void
t8dg_LGL_functionbasis_apply_derivative_matrix_transpose (sc_array_t * dof_values,
                                                          sc_array_t * derivative_dof_values, t8dg_LGL_functionbasis_t * functionbasis)
{
  T8_ASSERT (functionbasis->vertices->dim == 1);
  double             *dof_array;
  double             *derivative_array;
  dof_array = (double *) sc_array_index_int (dof_values, 0);
  derivative_array = (double *) sc_array_index_int (derivative_dof_values, 0);
  switch (functionbasis->vertices->number_of_vertices) {
  case (1):
    dof_array[0] = 0;
    break;
  case (2):
    dof_array[0] = 1 * derivative_array[0] - 1 * derivative_array[1];
    dof_array[1] = -1 * derivative_array[0] + 1 * derivative_array[1];
    break;
  case (3):
    dof_array[0] = -3 * derivative_array[0] + 4 * derivative_array[1] - 1 * derivative_array[2];
    dof_array[1] = -1 * derivative_array[0] + 0 * derivative_array[1] + 1 * derivative_array[2];
    dof_array[2] = 1 * derivative_array[0] - 4 * derivative_array[1] + 3 * derivative_array[2];
    break;
  default:
    T8DG_ASSERT (0);
  }
}

void
t8dg_LGL_vertexset_get_3D_vertex (double reference_vertex[3], t8dg_LGL_vertexset_t * vertexset, int ivertex)
{
  T8DG_ASSERT (vertexset->dim == 1);    /* other not yet implemented */
  double             *vertex;
  int                 idim;
  vertex = (double *) sc_array_index_int (vertexset->vertices, ivertex);
  for (idim = 0; idim < vertexset->dim; idim++) {
    reference_vertex[idim] = vertex[idim];
  }
  for (idim = vertexset->dim; idim < DIM3; idim++) {
    reference_vertex[idim] = 0;
  }
}

/*TODO: ASSERTS for all get functions!!! */

t8dg_quad_idx_t
t8dg_LGL_functionbasis_get_num_dof (t8dg_LGL_functionbasis_t * functionbasis)
{
  return functionbasis->number_of_dof;
}

void
t8dg_LGL_functionbasis_get_vertex (double vertex[3], t8dg_LGL_functionbasis_t * functionbasis, int idof)
{
  t8dg_LGL_vertexset_get_3D_vertex (vertex, functionbasis->vertices, idof);
}

int
t8dg_LGL_quadrature_get_num_faces (t8dg_LGL_quadrature_t * quadrature)
{
  return quadrature->vertices->number_of_faces;
}

t8dg_quad_idx_t
t8dg_LGL_quadrature_get_num_element_vertices (t8dg_LGL_quadrature_t * quadrature)
{
  return quadrature->number_of_quadrature_points;
}

t8dg_quad_idx_t
t8dg_LGL_quadrature_get_num_face_vertices (t8dg_LGL_quadrature_t * quadrature, int iface)
{
  return quadrature->vertices->number_of_facevertices[iface];
}

void
t8dg_LGL_quadrature_get_element_vertex (double vertex[3], t8dg_LGL_quadrature_t * quadrature, t8dg_quad_idx_t iquad)
{
  t8dg_LGL_vertexset_get_3D_vertex (vertex, quadrature->vertices, iquad);
}

double
t8dg_LGL_quadrature_get_element_weight (t8dg_LGL_quadrature_t * quadrature, t8dg_quad_idx_t iquad)
{
  return *(double *) t8dg_sc_array_index_quadidx (quadrature->weights, iquad);
}

void
t8dg_LGL_quadrature_get_face_vertex (double vertex[3], t8dg_LGL_quadrature_t * quadrature, int iface, t8dg_quad_idx_t ifacequad)
{
  t8dg_quad_idx_t     ielemquad;
  ielemquad = t8dg_LGL_vertexset_facequadidx_to_elementquadidx (quadrature, iface, ifacequad);
  t8dg_LGL_vertexset_get_3D_vertex (vertex, quadrature->vertices, ielemquad);
}

double
t8dg_LGL_quadrature_get_face_weight (t8dg_LGL_quadrature_t * quadrature, int iface, t8dg_quad_idx_t ifacequad)
{
  T8DG_ASSERT (quadrature->vertices->dim == 1);
  return 1;
}

void                t8dg_LGL_transform_element_dof_to_face_quad
  (sc_array_t * face_quad_array, const sc_array_t * element_dof_array, int iface, t8dg_LGL_quadrature_t * quadrature,
   t8dg_LGL_functionbasis_t * functionbasis)
{
  T8DG_ASSERT (quadrature->vertices == functionbasis->vertices);
  T8DG_ASSERT (quadrature->vertices->dim == 1);
  T8DG_ASSERT (element_dof_array->elem_size == sizeof (double));
  T8DG_ASSERT (face_quad_array->elem_size == sizeof (double));
  t8dg_quad_idx_t     iquad, ifacequad;
  for (ifacequad = 0; ifacequad < t8dg_LGL_quadrature_get_num_face_vertices (quadrature, iface); ifacequad++) {
    iquad = *(int *) t8dg_sc_array_index_quadidx (quadrature->vertices->facevertex_indices[iface], ifacequad);
    *(double *) t8dg_sc_array_index_quadidx (face_quad_array, ifacequad) =
      *(double *) t8dg_sc_array_index_quadidx (element_dof_array, iquad);
  }
}

void
t8dg_LGL_transform_face_quad_to_element_dof (sc_array_t * element_dof_array,
                                             const sc_array_t * face_quad_array,
                                             int iface, t8dg_LGL_quadrature_t * quadrature, t8dg_LGL_functionbasis_t * functionbasis)
{
  T8DG_ASSERT (quadrature->vertices == functionbasis->vertices);
  T8DG_ASSERT (quadrature->vertices->dim == 1);
  T8DG_ASSERT (element_dof_array->elem_size == sizeof (double));
  T8DG_ASSERT (face_quad_array->elem_size == sizeof (double));
  t8dg_quad_idx_t     iquad, ifacequad;
  for (iquad = 0; iquad < t8dg_LGL_quadrature_get_num_element_vertices (quadrature); iquad++) {
    *(double *) t8dg_sc_array_index_quadidx (element_dof_array, iquad) = 0;
  }
  for (ifacequad = 0; ifacequad < t8dg_LGL_quadrature_get_num_face_vertices (quadrature, iface); ifacequad++) {
    iquad = *(int *) t8dg_sc_array_index_quadidx (quadrature->vertices->facevertex_indices[iface], ifacequad);
    *(double *) t8dg_sc_array_index_quadidx (element_dof_array, iquad) =
      *(double *) t8dg_sc_array_index_quadidx (face_quad_array, ifacequad);
  }
}
