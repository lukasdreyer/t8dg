#include "t8dg.h"
#include "t8dg_quadrature.h"
#include "t8dg_vertexset.h"
#include "t8dg_refcount.h"
#include "t8dg_tensor.h"
#include <t8_eclass.h>

typedef struct t8dg_face_quadrature t8dg_face_quadrature_t;
/** Additionally to the vertices save the quadrature weights*/
struct t8dg_quadrature
{
  t8dg_refcount_t     rc;
  t8_eclass_t         element_class;
  t8dg_quadrature_type_t type;
  int                 embedded_dimension;

  t8dg_quadrature_t **face_quadrature;
  int                 num_face_quads; /**< the number of face quadratures */

  /*Those are only allocated for line, tri and tet and pyramid */
  t8dg_vertexset_t   *vertexset;                       /**< quadrature vertices, for LGL also gives lookup table for face vertices*/
  sc_array_t         *weights;                          /**< quadrature weights*/

  /*Otherwise, the element is in tensorform */
  int                 tensor;
  t8dg_quadrature_t  *tensor_first_quad;
  t8dg_quadrature_t  *tensor_second_quad;
};

static void
t8dg_quadrature_fill_0D_weights (t8dg_quadrature_t * quadrature)
{
  *(double *) sc_array_index (quadrature->weights, 0) = 1;
}

static void
t8dg_quadrature_fill_LGL_1D_weights (t8dg_quadrature_t * quadrature)
{
  T8DG_ASSERT (t8dg_vertexset_get_eclass (quadrature->vertexset) == T8_ECLASS_LINE);
  T8DG_ASSERT (t8dg_vertexset_get_type (quadrature->vertexset) == T8DG_VERT_LGL);
  double             *weights;
  weights = (double *) sc_array_index (quadrature->weights, 0);

  /* sum of element weights = 1 = Vol([0,1]) */
  switch (t8dg_vertexset_get_num_vertices (quadrature->vertexset)) {
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
  case (5):
    weights[0] = 1. / 20;
    weights[1] = 49. / 180;
    weights[2] = 16. / 45;
    weights[3] = 49. / 180;
    weights[4] = 1. / 20;
    break;
  case (6):
    weights[0] = 1. / 30;
    weights[1] = (14 - sqrt (7)) / 60;
    weights[2] = (14 + sqrt (7)) / 60;
    weights[3] = (14 + sqrt (7)) / 60;
    weights[4] = (14 - sqrt (7)) / 60;
    weights[5] = 1. / 30;
    break;
  case (7):
    weights[0] = 1. / 42;
    weights[1] = (124 - 7 * sqrt (15)) / 700;
    weights[2] = (124 + 7 * sqrt (15)) / 700;
    weights[3] = 128. / 525;
    weights[4] = (124 + 7 * sqrt (15)) / 700;
    weights[5] = (124 - 7 * sqrt (15)) / 700;
    weights[6] = 1. / 42;
    break;
  default:
    T8DG_ABORT ("Not yet implemented!\n");
  }
}

t8dg_quadrature_t  *
t8dg_quadrature_new_vertexset (t8dg_vertexset_t * vertexset, int create_face_quadrature)
{
  T8DG_ASSERT (vertexset != NULL);
  T8DG_ASSERT (t8dg_vertexset_get_type (vertexset) == T8DG_VERT_LGL);

  t8dg_quadrature_t  *quadrature = T8DG_ALLOC_ZERO (t8dg_quadrature_t, 1);

  t8dg_refcount_init (&quadrature->rc);

  quadrature->vertexset = vertexset;
  t8dg_vertexset_ref (vertexset);
  quadrature->embedded_dimension = t8dg_vertexset_get_embedded_dim (vertexset);

  quadrature->element_class = t8dg_vertexset_get_eclass (vertexset);
  quadrature->tensor = 0;

  quadrature->weights = sc_array_new_count (sizeof (double), t8dg_vertexset_get_num_vertices (vertexset));

  switch (quadrature->element_class) {
  case T8_ECLASS_VERTEX:
    t8dg_quadrature_fill_0D_weights (quadrature);
    break;
  case T8_ECLASS_LINE:
    t8dg_quadrature_fill_LGL_1D_weights (quadrature);
    break;
  case T8_ECLASS_TRIANGLE:
    T8DG_ABORT ("Not yet implemented");
//    t8dg_quadrature_fill_LGL_triangle_weights(quadrature);
    break;
  default:
    T8DG_ABORT ("Not yet implemented");
    break;
  }

  /*TODO: triangle */
  if (create_face_quadrature) {
    quadrature->num_face_quads = 2;
    quadrature->face_quadrature = T8DG_ALLOC (t8dg_quadrature_t *, quadrature->num_face_quads);
    if (t8dg_vertexset_get_type (vertexset) == T8DG_VERT_LGL) {
      t8dg_vertexset_t   *left_face_vertexset;
      t8dg_vertexset_t   *right_face_vertexset;
      left_face_vertexset = t8dg_vertexset_new_lgl_facevertexset (vertexset, 0);
      right_face_vertexset = t8dg_vertexset_new_lgl_facevertexset (vertexset, 1);
      quadrature->face_quadrature[0] = t8dg_quadrature_new_vertexset (left_face_vertexset, 0);
      quadrature->face_quadrature[1] = t8dg_quadrature_new_vertexset (right_face_vertexset, 0);
      t8dg_vertexset_unref (&left_face_vertexset);
      t8dg_vertexset_unref (&right_face_vertexset);
    }
    else {
      T8DG_ABORT ("Not yet implemented");
    }
  }
  else {
    quadrature->num_face_quads = 0;
    quadrature->face_quadrature = NULL;
  }
  return quadrature;
}

t8dg_quadrature_t  *
t8dg_quadrature_new_tensor (t8dg_quadrature_t * tensor_first_quad, t8dg_quadrature_t * tensor_second_quad, int create_face_quadrature)
{
  t8dg_quadrature_t  *quadrature = T8DG_ALLOC_ZERO (t8dg_quadrature_t, 1);
  t8dg_refcount_init (&quadrature->rc);

  quadrature->element_class = t8dg_tensor_eclass (tensor_first_quad->element_class, tensor_second_quad->element_class);
  quadrature->tensor = 1;

  if (tensor_first_quad->type == T8DG_QUAD_LGL && tensor_second_quad->type == T8DG_QUAD_LGL) {
    quadrature->type = T8DG_QUAD_LGL;
  }
  else {
    T8DG_ABORT ("Not yet implemented");
  }

  quadrature->tensor_first_quad = tensor_first_quad;
  quadrature->tensor_second_quad = tensor_second_quad;
  t8dg_quadrature_ref (tensor_first_quad);
  t8dg_quadrature_ref (tensor_second_quad);
  quadrature->embedded_dimension = tensor_first_quad->embedded_dimension + tensor_second_quad->embedded_dimension;

  if (create_face_quadrature) {
    int                 iface;
    T8DG_ASSERT (tensor_first_quad->num_face_quads > 0 && tensor_second_quad->num_face_quads > 0);
    quadrature->num_face_quads = tensor_first_quad->num_face_quads + tensor_second_quad->num_face_quads;
    quadrature->face_quadrature = T8DG_ALLOC (t8dg_quadrature_t *, quadrature->num_face_quads);
    for (iface = 0; iface < tensor_first_quad->num_face_quads; iface++) {
      quadrature->face_quadrature[iface] = t8dg_quadrature_new_tensor (tensor_first_quad->face_quadrature[iface], tensor_second_quad, 0);
    }
    for (iface = 0; iface < tensor_second_quad->num_face_quads; iface++) {
      quadrature->face_quadrature[iface + tensor_first_quad->num_face_quads] =
        t8dg_quadrature_new_tensor (tensor_first_quad, tensor_second_quad->face_quadrature[iface], 0);
    }
  }
  else {
    quadrature->num_face_quads = 0;
    quadrature->face_quadrature = NULL;
  }

  return quadrature;
}

t8dg_quadrature_t  *
t8dg_quadrature_new_hypercube (int dim, t8dg_vertexset_t * vertexset1D, int create_face_quadrature)
{
  T8DG_ASSERT (vertexset1D != NULL);
  T8DG_ASSERT (t8dg_vertexset_get_eclass (vertexset1D) == T8_ECLASS_LINE);

  t8dg_quadrature_t  *quadrature1D;
  quadrature1D = t8dg_quadrature_new_vertexset (vertexset1D, create_face_quadrature);
  if (dim == 1) {
    return quadrature1D;
  }
  t8dg_quadrature_t  *quadrature2D;
  quadrature2D = t8dg_quadrature_new_tensor (quadrature1D, quadrature1D, create_face_quadrature);
  if (dim == 2) {
    t8dg_quadrature_unref (&quadrature1D);
    return quadrature2D;
  }

  t8dg_quadrature_t  *quadrature3D;
  quadrature3D = t8dg_quadrature_new_tensor (quadrature2D, quadrature1D, create_face_quadrature);

  t8dg_quadrature_unref (&quadrature2D);
  t8dg_quadrature_unref (&quadrature1D);

  return quadrature3D;
}

int
t8dg_quadrature_get_num_face_quadrature (const t8dg_quadrature_t * quadrature)
{
  T8DG_ASSERT (quadrature != NULL);
  return quadrature->num_face_quads;
}

int
t8dg_quadrature_get_num_element_vertices (const t8dg_quadrature_t * quadrature)
{
  T8DG_ASSERT (quadrature != NULL);

  if (quadrature->tensor) {
    return t8dg_quadrature_get_num_element_vertices (quadrature->tensor_first_quad) *
      t8dg_quadrature_get_num_element_vertices (quadrature->tensor_second_quad);
  }
  return t8dg_vertexset_get_num_vertices (quadrature->vertexset);
}

int
t8dg_quadrature_get_dim (const t8dg_quadrature_t * quadrature)
{
  T8DG_ASSERT (quadrature != NULL);
  return t8_eclass_to_dimension[quadrature->element_class];
}

static void
t8dg_quadrature_fill_element_vertex (const t8dg_quadrature_t * quadrature, const int iquad, const int startdim, double vertex[3])
{
  T8DG_ASSERT (iquad >= 0 && iquad < t8dg_quadrature_get_num_element_vertices (quadrature));
  if (quadrature->tensor) {
    t8dg_quadrature_fill_element_vertex (quadrature->tensor_first_quad,
                                         iquad % t8dg_quadrature_get_num_element_vertices (quadrature->tensor_first_quad), startdim,
                                         vertex);
    t8dg_quadrature_fill_element_vertex (quadrature->tensor_second_quad,
                                         iquad / t8dg_quadrature_get_num_element_vertices (quadrature->tensor_first_quad),
                                         startdim + quadrature->tensor_first_quad->embedded_dimension, vertex);
  }
  else {
    t8dg_vertexset_fill_vertex3D (quadrature->vertexset, iquad, startdim, vertex);
  }
}

void
t8dg_quadrature_get_element_vertex (const t8dg_quadrature_t * quadrature, const int iquad, double vertex[3])
{
  T8DG_ASSERT (iquad >= 0 && iquad < t8dg_quadrature_get_num_element_vertices (quadrature));
  int                 startdim = 0;

  vertex[0] = vertex[1] = vertex[2] = 0;

  t8dg_quadrature_fill_element_vertex (quadrature, iquad, startdim, vertex);
}

double
t8dg_quadrature_get_element_weight (const t8dg_quadrature_t * quadrature, const int iquad)
{
  T8DG_ASSERT (quadrature != NULL);

  if (quadrature->tensor) {
    return t8dg_quadrature_get_element_weight (quadrature->tensor_first_quad,
                                               iquad % t8dg_quadrature_get_num_element_vertices (quadrature->tensor_first_quad)) *
      t8dg_quadrature_get_element_weight (quadrature->tensor_second_quad,
                                          iquad / t8dg_quadrature_get_num_element_vertices (quadrature->tensor_first_quad));
  }
  return *(double *) sc_array_index_int (quadrature->weights, iquad);
}

int
t8dg_quadrature_get_num_face_vertices (const t8dg_quadrature_t * quadrature, const int iface)
{
  T8DG_ASSERT (quadrature != NULL);
  T8DG_ASSERT (iface >= 0 && iface < quadrature->num_face_quads);
  return t8dg_quadrature_get_num_element_vertices (quadrature->face_quadrature[iface]);
}

void
t8dg_quadrature_get_face_vertex (const t8dg_quadrature_t * quadrature, const int iface, const int ifacequad, double vertex[3])
{
  T8DG_ASSERT (quadrature != NULL);
  T8DG_ASSERT (iface >= 0 && iface < quadrature->num_face_quads);
  t8dg_quadrature_get_element_vertex (quadrature->face_quadrature[iface], ifacequad, vertex);
}

double
t8dg_quadrature_get_face_weight (const t8dg_quadrature_t * quadrature, const int iface, const int ifacequad)
{
  T8DG_ASSERT (iface >= 0 && iface < quadrature->num_face_quads);
  return t8dg_quadrature_get_element_weight (quadrature->face_quadrature[iface], ifacequad);
}

t8dg_quadrature_type_t
t8dg_quadrature_get_type (const t8dg_quadrature_t * quadrature)
{
  return quadrature->type;
}

t8_eclass_t
t8dg_quadrature_get_eclass (const t8dg_quadrature_t * quadrature)
{
  return quadrature->element_class;
}

double
t8dg_quadrature_integrate_reference_element (t8dg_quadrature_t * quadrature, t8dg_scalar_function_3d_fn integrand_fn, void *integrand_data)
{
  int                 num_elem_quad, iquad;
  double              integral = 0;
  double              vertex[3];
  num_elem_quad = t8dg_quadrature_get_num_element_vertices (quadrature);
  for (iquad = 0; iquad < num_elem_quad; iquad++) {
    t8dg_quadrature_get_element_vertex (quadrature, iquad, vertex);
    integral += t8dg_quadrature_get_element_weight (quadrature, iquad) * integrand_fn (vertex, integrand_data);
  }
  return integral;
}

void
t8dg_quadrature_ref (t8dg_quadrature_t * quadrature)
{
  T8DG_ASSERT (quadrature != NULL);
  t8dg_refcount_ref (&quadrature->rc);
}

void
t8dg_quadrature_unref (t8dg_quadrature_t ** pquadrature)
{
  t8dg_quadrature_t  *quadrature;

  T8DG_ASSERT (pquadrature != NULL);
  quadrature = *pquadrature;
  T8DG_ASSERT (quadrature != NULL);
  T8DG_ASSERT (t8dg_refcount_is_active (&quadrature->rc));

  if (t8dg_refcount_unref (&quadrature->rc)) {
    t8dg_quadrature_reset (pquadrature);
  }
}

void
t8dg_quadrature_destroy (t8dg_quadrature_t ** pquadrature)
{
  t8dg_quadrature_t  *quadrature;

  T8DG_ASSERT (pquadrature != NULL);
  quadrature = *pquadrature;

  T8DG_ASSERT (quadrature != NULL);
  T8DG_ASSERT (t8dg_refcount_is_last (&quadrature->rc));

  t8dg_refcount_unref (&quadrature->rc);
  t8dg_quadrature_reset (pquadrature);
}

void
t8dg_quadrature_reset (t8dg_quadrature_t ** pquadrature)
{
  int                 iface;
  T8DG_ASSERT (pquadrature != NULL);
  t8dg_quadrature_t  *quadrature = *pquadrature;
  T8DG_ASSERT (quadrature != NULL);

  for (iface = 0; iface < quadrature->num_face_quads; iface++) {
    t8dg_quadrature_destroy (&quadrature->face_quadrature[iface]);
  }
  if (quadrature->num_face_quads > 0) {
    T8DG_FREE (quadrature->face_quadrature);
  }

  if (quadrature->tensor) {
    t8dg_quadrature_unref (&quadrature->tensor_first_quad);
    t8dg_quadrature_unref (&quadrature->tensor_second_quad);
  }
  else {
    sc_array_destroy (quadrature->weights);
    t8dg_vertexset_unref (&quadrature->vertexset);
  }

  T8_FREE (quadrature);
  *pquadrature = NULL;
}
