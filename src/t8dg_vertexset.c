
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
  int                 dim;      /**< Dimension of the reference element*/
  int                 number_of_vertices;       /**< Number of element vertices*/
  int                 number_of_faces;          /**< Number of faces */
  int                 number_of_facevertices[MAX_FACES];        /**< For each face, the number of face vertices*/
  sc_array_t         *vertices; /**< size: dim * number_of_vertices, make access available via function and allocate only if not tensor? */
  sc_array_t         *facevertex_indices[MAX_FACES];    /**< Lookup table, for each facequadindex save the elementquadindex*/
#if 0
  int                 tensorflag;
  t8dg_vertexset_t   *tensor1;
  t8dg_vertexset_t   *tensor2;
#endif
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
  return vertexset->dim;
}

int
t8dg_vertexset_get_num_element_vertices (const t8dg_vertexset_t * vertexset)
{
  T8DG_ASSERT (vertexset != NULL);
  return vertexset->number_of_vertices;
}

int
t8dg_vertexset_get_num_faces (const t8dg_vertexset_t * vertexset)
{
  T8DG_ASSERT (vertexset != NULL);
  return vertexset->number_of_faces;
}

int
t8dg_vertexset_get_num_face_vertices (const t8dg_vertexset_t * vertexset, const int iface)
{
  T8DG_ASSERT (vertexset != NULL);
  T8DG_ASSERT (iface >= 0 && iface < vertexset->number_of_faces);
  return vertexset->number_of_facevertices[iface];
}

int
t8dg_vertexset_get_LGL_facevertex_element_index (const t8dg_vertexset_t * vertexset, const int iface, const int ifacevertex)
{
  T8DG_ASSERT (vertexset != NULL);
  T8DG_ASSERT (t8dg_vertexset_get_type (vertexset) == T8DG_LGL);
  T8DG_ASSERT (iface >= 0 && iface < t8dg_vertexset_get_num_faces (vertexset));
  T8DG_ASSERT (ifacevertex >= 0 && ifacevertex < t8dg_vertexset_get_num_face_vertices (vertexset, iface));

  return *(int *) sc_array_index_int (vertexset->facevertex_indices[iface], ifacevertex);
}

void
t8dg_vertexset_get_vertex (double reference_vertex[3], const t8dg_vertexset_t * vertexset, const int ivertex)
{
  double             *vertex;
  int                 idim;

  T8DG_ASSERT (vertexset != NULL);
  T8DG_ASSERT (ivertex >= 0 && ivertex < vertexset->number_of_vertices);

  vertex = (double *) sc_array_index_int (vertexset->vertices, ivertex);
  for (idim = 0; idim < vertexset->dim; idim++) {
    reference_vertex[idim] = vertex[idim];
  }
  for (idim = vertexset->dim; idim < DIM3; idim++) {
    reference_vertex[idim] = 0;
  }
}

static void
t8dg_vertexset_set_facevertex_index (t8dg_vertexset_t * vertices, const int iface, const int ifacevertex, const int ielemvertex)
{
  T8DG_ASSERT (vertices != NULL);
  T8DG_ASSERT (iface >= 0 && iface < t8dg_vertexset_get_num_faces (vertices));
  T8DG_ASSERT (ifacevertex >= 0 && ifacevertex < t8dg_vertexset_get_num_face_vertices (vertices, iface));
  T8DG_ASSERT (ielemvertex >= 0 && ielemvertex < t8dg_vertexset_get_num_element_vertices (vertices));

  *(int *) sc_array_index_int (vertices->facevertex_indices[iface], ifacevertex) = ielemvertex;
}

t8dg_vertexset_t   *
t8dg_vertexset_new_1D_LGL (const int number_of_LGL_vertices)
{
  SC_CHECK_ABORT (number_of_LGL_vertices >= 1 && number_of_LGL_vertices <= 4, "Only up to 4 LGL vertices implemented"); /*Larger not implemented */

  int                 iface;
  t8dg_vertexset_t   *vertices;
  vertices = T8DG_ALLOC_ZERO (t8dg_vertexset_t, 1);

  t8dg_refcount_init (&vertices->rc);
  vertices->type = T8DG_LGL;

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

void
t8dg_vertexset_reset (t8dg_vertexset_t ** pvertexset)
{
  int                 iface = 0;
  t8dg_vertexset_t   *vertexset;
  T8DG_ASSERT (pvertexset != NULL);
  vertexset = *pvertexset;
  T8DG_ASSERT (vertexset != NULL);
  T8DG_ASSERT (vertexset->rc.refcount == 0);

  if (t8dg_vertexset_get_type (vertexset) == T8DG_LGL) {
    for (iface = 0; iface < vertexset->number_of_faces; iface++) {
      sc_array_destroy (vertexset->facevertex_indices[iface]);
      vertexset->facevertex_indices[iface] = NULL;
    }
  }
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
