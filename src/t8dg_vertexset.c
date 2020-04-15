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
  t8_eclass_t         element_class;    /* TODO: replace dim and number_of_faces by element_class */

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
t8dg_vertexset_get_dim (const t8dg_vertexset_t * vertexset)
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

  vertex = (double *) sc_array_index_int (vertexset->vertices, ivertex);
  enddim = startdim + t8dg_vertexset_get_dim (vertexset);
  for (idim = startdim; idim < enddim; idim++) {
    reference_vertex[idim] = vertex[idim - startdim];
  }
}

t8dg_vertexset_t   *
t8dg_vertexset_new_1D_LGL (const int number_of_LGL_vertices)
{
  SC_CHECK_ABORT (number_of_LGL_vertices >= 1 && number_of_LGL_vertices <= 4, "Only up to 4 LGL vertices implemented"); /*Larger not implemented */

  t8dg_vertexset_t   *vertices;
  vertices = T8DG_ALLOC_ZERO (t8dg_vertexset_t, 1);

  t8dg_refcount_init (&vertices->rc);
  vertices->type = T8DG_VERT_LGL;
  vertices->element_class = T8_ECLASS_LINE;

  vertices->number_of_vertices = number_of_LGL_vertices;
  vertices->vertices = sc_array_new_count (sizeof (double), vertices->number_of_vertices);

  double             *vertex_array;
  vertex_array = (double *) sc_array_index (vertices->vertices, 0);

  /*LGL vertices on [0,1] */
  switch (vertices->number_of_vertices) {
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
  default:
    printf ("Not yet implemented!\n");
    T8DG_ASSERT (0);
  }

  return vertices;
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
