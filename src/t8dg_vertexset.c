/*
 *  This file is part of t8dg.
 *  t8dg is a discontinuous galerkin solver on adaptive meshes
 *  that uses the adaptive meshing library t8code to distribute
 *  the workload equally among parallel processes to achieve good scaling.
 *
 *  Copyright (C) 2020 the developers
 *
 *  t8dg is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  t8dg is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with t8dg.  If not, see <http://www.gnu.org/licenses/>.
 */

#include <t8_eclass.h>
#include "t8dg_vertexset.h"
#include "t8dg.h"
#include "t8dg_refcount.h"

/** Vertex set used for the quadrature and functionbasis
 * Currently only for 1D
 */
struct t8dg_vertexset
{
  t8dg_refcount_t     rc;

  t8dg_vertexset_type_t type;
  t8_eclass_t         element_class;

  int                 embedded_dimensions;
  int                 num_children;

  int                 number_of_vertices;       /**< Number of element vertices*/
  sc_array_t         *vertices; /**< size: dim * number_of_vertices, make access available via function and allocate only if not tensor? */
};

t8dg_vertexset_type_t
t8dg_vertexset_get_type (const t8dg_vertexset_t * vertexset)
{
  T8DG_ASSERT (vertexset != NULL);
  return vertexset->type;
}

int
t8dg_vertexset_get_embedded_dim (const t8dg_vertexset_t * vertexset)
{
  T8DG_ASSERT (vertexset != NULL);
  return vertexset->embedded_dimensions;
}

int
t8dg_vertexset_get_eclass_dim (const t8dg_vertexset_t * vertexset)
{
  T8DG_ASSERT (vertexset != NULL);
  return t8_eclass_to_dimension[vertexset->element_class];
}

int
t8dg_vertexset_get_num_vertices (const t8dg_vertexset_t * vertexset)
{
  T8DG_ASSERT (vertexset != NULL);
  return vertexset->number_of_vertices;
}

t8_eclass_t
t8dg_vertexset_get_eclass (const t8dg_vertexset_t * vertexset)
{
  T8DG_ASSERT (vertexset != NULL);
  return vertexset->element_class;
}

/*fills a 3dim vertex at the correct dimensions, unused dims need to be set to 0*/
void
t8dg_vertexset_fill_vertex3D (const t8dg_vertexset_t * vertexset, const int ivertex, const int startdim, double reference_vertex[3])
{
  double             *vertex;
  int                 idim, enddim;

  T8DG_ASSERT (vertexset != NULL);
  T8DG_ASSERT (ivertex >= 0 && ivertex < vertexset->number_of_vertices);
  T8DG_ASSERT (startdim >= 0);
  T8DG_ASSERT (startdim + vertexset->embedded_dimensions <= 3);

  vertex = (double *) sc_array_index_int (vertexset->vertices, ivertex);
  enddim = startdim + vertexset->embedded_dimensions;
  for (idim = startdim; idim < enddim; idim++) {
    reference_vertex[idim] = vertex[idim - startdim];
  }
}

double
t8dg_vertexset_get_first_coordinate (const t8dg_vertexset_t * vertexset, const int ivertex)
{
  return *(double *) sc_array_index_int (vertexset->vertices, ivertex);
}

t8dg_vertexset_t   *
t8dg_vertexset_new_1D_LGL (const int number_of_LGL_vertices)
{
  t8dg_vertexset_t   *vertexset;
  vertexset = T8DG_ALLOC_ZERO (t8dg_vertexset_t, 1);

  t8dg_refcount_init (&vertexset->rc);
  vertexset->type = T8DG_VERT_LGL;
  vertexset->element_class = T8_ECLASS_LINE;
  vertexset->embedded_dimensions = 1;
  vertexset->num_children = 2;

  vertexset->number_of_vertices = number_of_LGL_vertices;
  vertexset->vertices = sc_array_new_count (sizeof (double), vertexset->number_of_vertices);

  double             *vertex_array;
  vertex_array = (double *) sc_array_index (vertexset->vertices, 0);

  /*LGL vertices on [0,1] */
  switch (vertexset->number_of_vertices) {
  case (1):
    vertex_array[0] = .5;
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
    vertex_array[1] = (1 - sqrt (1. / 5)) / 2;
    vertex_array[2] = (1 + sqrt (1. / 5)) / 2;
    vertex_array[3] = 1;
    break;
  case (5):
    vertex_array[0] = 0;
    vertex_array[1] = (1 - sqrt (3. / 7)) / 2;
    vertex_array[2] = 0.5;
    vertex_array[3] = (1 + sqrt (3. / 7)) / 2;
    vertex_array[4] = 1;
    break;
  case (6):
    vertex_array[0] = 0;
    vertex_array[1] = (1 - sqrt (1. / 3 + 2 * sqrt (7) / 21)) / 2;
    vertex_array[2] = (1 - sqrt (1. / 3 - 2 * sqrt (7) / 21)) / 2;
    vertex_array[3] = (1 + sqrt (1. / 3 - 2 * sqrt (7) / 21)) / 2;
    vertex_array[4] = (1 + sqrt (1. / 3 + 2 * sqrt (7) / 21)) / 2;
    vertex_array[5] = 1;
    break;
  case (7):
    vertex_array[0] = 0;
    vertex_array[1] = (1 - sqrt ((5. + 2 * sqrt (5. / 3)) / 11)) / 2;
    vertex_array[2] = (1 - sqrt ((5. - 2 * sqrt (5. / 3)) / 11)) / 2;
    vertex_array[3] = 0.5;
    vertex_array[4] = (1 + sqrt ((5. - 2 * sqrt (5. / 3)) / 11)) / 2;
    vertex_array[5] = (1 + sqrt ((5. + 2 * sqrt (5. / 3)) / 11)) / 2;
    vertex_array[6] = 1;
    break;

  default:
    T8DG_ABORT ("Not yet implemented!");
  }

  return vertexset;
}

t8dg_vertexset_t   *
t8dg_vertexset_new_lgl_facevertexset (const t8dg_vertexset_t * vertexset, int iface)
{
  T8DG_ASSERT (vertexset->type == T8DG_VERT_LGL);
  t8dg_vertexset_t   *face_vertexset;
  face_vertexset = T8DG_ALLOC_ZERO (t8dg_vertexset_t, 1);
  face_vertexset->embedded_dimensions = vertexset->embedded_dimensions;
  face_vertexset->type = T8DG_VERT_LGL;
  t8dg_refcount_init (&face_vertexset->rc);
  switch (vertexset->element_class) {
  case T8_ECLASS_LINE:
    face_vertexset->element_class = T8_ECLASS_VERTEX;
    face_vertexset->num_children = 0;   /*or 1? */
    face_vertexset->number_of_vertices = 1;
    face_vertexset->vertices = sc_array_new_count (sizeof (double), 1);
    *(double *) sc_array_index (face_vertexset->vertices, 0) =
      *(double *) sc_array_index (vertexset->vertices, iface ? vertexset->number_of_vertices - 1 : 0);
    break;
  case T8_ECLASS_TRIANGLE:
    T8DG_ABORT ("Not yet implemented!");
    break;
  default:
    T8DG_ABORT ("Not yet implemented!");
    break;
  }
  return face_vertexset;
}

t8dg_vertexset_t   *
t8dg_vertexset_new_childvertexset (const t8dg_vertexset_t * vertexset, int ichild)
{
  T8DG_ASSERT (ichild >= 0 || ichild <= vertexset->num_children);

  int                 ivertex;

  t8dg_vertexset_t   *child_vertexset;
  child_vertexset = T8DG_ALLOC_ZERO (t8dg_vertexset_t, 1);

  t8dg_refcount_init (&child_vertexset->rc);
  child_vertexset->type = vertexset->type;
  child_vertexset->element_class = vertexset->element_class;
  child_vertexset->number_of_vertices = t8dg_vertexset_get_num_vertices (vertexset);
  child_vertexset->embedded_dimensions = vertexset->embedded_dimensions;
  child_vertexset->num_children = vertexset->num_children;

  child_vertexset->vertices = sc_array_new_count (sizeof (double), child_vertexset->number_of_vertices);

  switch (child_vertexset->element_class) {
  case T8_ECLASS_LINE:
    for (ivertex = 0; ivertex < child_vertexset->number_of_vertices; ivertex++) {
      *(double *) sc_array_index_int (child_vertexset->vertices, ivertex) =
        ichild / 2.0 + t8dg_vertexset_get_first_coordinate (vertexset, ivertex) / 2.0;
    }
    break;
  case T8_ECLASS_TRIANGLE:
    T8DG_ABORT ("Not yet implemented!");
    break;

  default:
    T8DG_ABORT ("Not yet implemented!");
    break;
  }
  return child_vertexset;
}

void
t8dg_vertexset_reset (t8dg_vertexset_t ** pvertexset)
{
  t8dg_vertexset_t   *vertexset;
  T8DG_ASSERT (pvertexset != NULL);
  vertexset = *pvertexset;
  T8DG_ASSERT (vertexset != NULL);
  T8DG_ASSERT (vertexset->rc.refcount == 0);

  if (vertexset->vertices != NULL) {
    sc_array_destroy (vertexset->vertices);
    vertexset->vertices = NULL;
  }
  T8DG_FREE (vertexset);
  *pvertexset = NULL;
}

void
t8dg_vertexset_ref (t8dg_vertexset_t * vertexset)
{
  T8DG_ASSERT (vertexset != NULL);
  t8dg_refcount_ref (&vertexset->rc);
}

void
t8dg_vertexset_unref (t8dg_vertexset_t ** pvertexset)
{
  t8dg_vertexset_t   *vertexset;

  T8DG_ASSERT (pvertexset != NULL);
  vertexset = *pvertexset;
  T8DG_ASSERT (vertexset != NULL);
  T8DG_ASSERT (t8dg_refcount_is_active (&vertexset->rc));

  if (t8dg_refcount_unref (&vertexset->rc)) {
    t8dg_vertexset_reset (pvertexset);
  }
}

void
t8dg_vertexset_destroy (t8dg_vertexset_t ** pvertexset)
{
  t8dg_vertexset_t   *vertexset;

  T8DG_ASSERT (pvertexset != NULL);
  vertexset = *pvertexset;

  T8DG_ASSERT (vertexset != NULL);
  T8DG_ASSERT (t8dg_refcount_is_last (&vertexset->rc));

  t8dg_refcount_unref (&vertexset->rc);
  t8dg_vertexset_reset (pvertexset);
}
