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

  /*Those are only allocated for line, tri and tet and pyramid */
  t8dg_vertexset_t   *vertexset;                       /**< quadrature vertices, for LGL also gives lookup table for face vertices*/
  sc_array_t         *weights;                          /**< quadrature weights*/
//  t8dg_face_quadrature_t *face_quadrature[MAX_FACES];

  /*Otherwise, the element is in tensorform */
  int                 num_tensor;
  t8dg_quadrature_t  *tensor_quad[DIM3];
  int                 tensor_num_vertices[DIM3];

  /* Information for Lookuptable LGL element->face, only filled for LGL */
  int                 beginning_index[MAX_FACES];
  int                 stride[3];
};

#if 0
struct t8dg_face_quadrature
{
  sc_array_t         *ifacequad_to_ielemquad;
  t8dg_vertexset_t   *vertexset;
  sc_array_t         *weights;
};
#endif

void               *
t8dg_sc_array_index_quadidx (const sc_array_t * array, t8dg_quad_idx_t iz)
{
  T8DG_ASSERT (array != NULL);
  T8DG_ASSERT (iz >= 0 && iz < (t8dg_quad_idx_t) array->elem_count);

  return (void *) (array->array + (array->elem_size * iz));
}

static void
t8dg_quadrature_fill_LGL_1D_weights (t8dg_quadrature_t * quadrature)
{
  T8DG_ASSERT (t8dg_vertexset_get_dim (quadrature->vertexset) == 1);
  T8DG_ASSERT (t8dg_vertexset_get_type (quadrature->vertexset) == T8DG_VERT_LGL);
  double             *weights;
  weights = (double *) t8_sc_array_index_locidx (quadrature->weights, 0);

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
t8dg_quadrature_new_vertexset (t8dg_vertexset_t * vertexset)
{
  T8DG_ASSERT (vertexset != NULL);
  T8DG_ASSERT (t8dg_vertexset_get_type (vertexset) == T8DG_VERT_LGL);

  t8dg_quadrature_t  *quadrature = T8DG_ALLOC_ZERO (t8dg_quadrature_t, 1);

  t8dg_refcount_init (&quadrature->rc);

  quadrature->vertexset = vertexset;
  t8dg_vertexset_ref (vertexset);

  quadrature->element_class = t8dg_vertexset_get_eclass (vertexset);
  quadrature->num_tensor = 1;

  quadrature->weights = sc_array_new_count (sizeof (double), t8dg_vertexset_get_num_vertices (vertexset));
  switch (t8dg_vertexset_get_type (vertexset)) {
  case T8DG_VERT_LGL:
    t8dg_quadrature_fill_LGL_1D_weights (quadrature);
    quadrature->beginning_index[0] = 0;
    quadrature->beginning_index[1] = t8dg_vertexset_get_num_vertices (vertexset) - 1;
    break;
  default:
    T8DG_ABORT ("Not yet implemented");
  }
  /* TODO: Facequadrature if not LGL */

  return quadrature;
}

t8dg_quadrature_t  *
t8dg_quadrature_new_tensor (int num_tensor, t8dg_quadrature_t * quad_tensors[3])
{
  int                 itensor;
  t8dg_quadrature_t  *quadrature = T8DG_ALLOC_ZERO (t8dg_quadrature_t, 1);
  t8dg_refcount_init (&quadrature->rc);
  if (num_tensor == 2) {
    quadrature->element_class = T8_ECLASS_QUAD;
    //could also be prism
  }
  else {
    quadrature->element_class = T8_ECLASS_HEX;
  }
  quadrature->num_tensor = num_tensor;
  quadrature->type = quad_tensors[0]->type;
  for (itensor = 0; itensor < num_tensor; itensor++) {
    quadrature->tensor_num_vertices[itensor] = t8dg_quadrature_get_num_element_vertices (quad_tensors[itensor]);
    quadrature->tensor_quad[itensor] = quad_tensors[itensor];
    t8dg_quadrature_ref (quad_tensors[itensor]);
    if (quadrature->type != quad_tensors[itensor]->type) {
      quadrature->type = T8_QUAD_UNKNOWN;
    }
  }
  if (quadrature->type == T8DG_QUAD_LGL) {
    int                 stride = 1;
    for (itensor = 0; itensor < num_tensor; itensor++) {
      quadrature->beginning_index[2 * itensor] = 0;
      quadrature->beginning_index[2 * itensor + 1] = stride * (quadrature->tensor_num_vertices[itensor] - 1);
      quadrature->stride[itensor] = stride;
      stride *= quadrature->tensor_num_vertices[itensor];
    }
  }
  return quadrature;
}

t8dg_quadrature_t  *
t8dg_quadrature_new_hypercube (int dim, t8dg_vertexset_t * vertexset1D)
{
  T8DG_ASSERT (vertexset1D != NULL);
  T8DG_ASSERT (t8dg_vertexset_get_eclass (vertexset1D) == T8_ECLASS_LINE);

  t8dg_quadrature_t  *quadrature1D;
  quadrature1D = t8dg_quadrature_new_vertexset (vertexset1D);
  if (dim == 1)
    return quadrature1D;

  t8dg_quadrature_t  *quadrature_hypercube;
  t8dg_quadrature_t  *quad_tensors[3] = { NULL, NULL, NULL };
  int                 itensor;
  for (itensor = 0; itensor < dim; itensor++) {
    quad_tensors[itensor] = quadrature1D;
  }

  quadrature_hypercube = t8dg_quadrature_new_tensor (dim, quad_tensors);
  t8dg_quadrature_unref (&quadrature1D);

  return quadrature_hypercube;
}

t8dg_quad_idx_t
t8dg_quadrature_lgl_facequadix_lookup (const t8dg_quadrature_t * quadrature, const int iface, const t8dg_quad_idx_t ifacequad)
{
  int                 ifirst, isecond;
  switch (iface / 2) {
  case 0:
    ifirst = 1;
    isecond = 2;
    break;
  case 1:
    ifirst = 0;
    isecond = 2;
    break;
  case 2:
    ifirst = 0;
    isecond = 1;
    break;
  default:
    T8DG_ABORT ("Facenumber too big");
  }

  switch (quadrature->num_tensor) {
  case 1:
    return quadrature->beginning_index[iface];
  case 2:
    return quadrature->beginning_index[iface] + quadrature->stride[ifirst] * (ifacequad % quadrature->tensor_num_vertices[ifirst]);
  case 3:
    return quadrature->beginning_index[iface] + quadrature->stride[ifirst] * (ifacequad % quadrature->tensor_num_vertices[ifirst])
      + quadrature->stride[isecond] * (ifacequad / quadrature->tensor_num_vertices[ifirst]);
  default:
    T8DG_ABORT ("Not supported");
  }
}

int
t8dg_quadrature_get_num_faces (const t8dg_quadrature_t * quadrature)
{
  T8DG_ASSERT (quadrature != NULL);
  return t8_eclass_num_faces[quadrature->element_class];
}

t8dg_quad_idx_t
t8dg_quadrature_get_num_element_vertices (const t8dg_quadrature_t * quadrature)
{
  int                 itensor;
  int                 num_vertices;
  T8DG_ASSERT (quadrature != NULL);
  if (quadrature->num_tensor == 1) {
    return t8dg_vertexset_get_num_vertices (quadrature->vertexset);
  }
  else {
    num_vertices = 1;
    for (itensor = 0; itensor < quadrature->num_tensor; itensor++) {
      num_vertices *= t8dg_quadrature_get_num_element_vertices (quadrature->tensor_quad[itensor]);
    }
    return num_vertices;
  }
}

int
t8dg_quadrature_get_dim (const t8dg_quadrature_t * quadrature)
{
  T8DG_ASSERT (quadrature != NULL);
  return t8_eclass_to_dimension[quadrature->element_class];
}

void
t8dg_quadrature_get_element_vertex (double vertex[3], const t8dg_quadrature_t * quadrature, const t8dg_quad_idx_t iquad)
{
  vertex[0] = 0;
  vertex[1] = 0;
  vertex[2] = 0;
  T8DG_ASSERT (quadrature != NULL);
  int                 itensor;
  int                 iquadtensor[DIM3];
  int                 startdim = 0;
  if (quadrature->num_tensor == 1) {
    t8dg_vertexset_fill_vertex3D (quadrature->vertexset, iquad, 0, vertex);
  }
  else {
    t8dg_tensor_transform_tensoridx (iquad, quadrature->tensor_num_vertices, iquadtensor);
    for (itensor = 0; itensor < quadrature->num_tensor; itensor++) {
      t8dg_vertexset_fill_vertex3D (quadrature->tensor_quad[itensor]->vertexset, iquadtensor[itensor], startdim, vertex);
      startdim += t8dg_quadrature_get_dim (quadrature->tensor_quad[itensor]);
    }
  }
}

double
t8dg_quadrature_get_element_weight (const t8dg_quadrature_t * quadrature, const t8dg_quad_idx_t iquad)
{
  double              weight = 1;
  int                 itensor;
  int                 iquadtensor[3];

  T8DG_ASSERT (quadrature != NULL);
  if (quadrature->num_tensor == 1) {
    weight = *(double *) t8dg_sc_array_index_quadidx (quadrature->weights, iquad);
  }
  else {
    t8dg_tensor_transform_tensoridx (iquad, quadrature->tensor_num_vertices, iquadtensor);
    for (itensor = 0; itensor < quadrature->num_tensor; itensor++) {
      weight *= t8dg_quadrature_get_element_weight (quadrature->tensor_quad[itensor], iquadtensor[itensor]);
    }
  }
  return weight;
}

t8dg_quad_idx_t
t8dg_quadrature_get_num_face_vertices (const t8dg_quadrature_t * quadrature, const int iface)
{
   /*TODO*/ T8DG_ASSERT (quadrature != NULL);
  if (t8dg_quadrature_get_dim (quadrature) == 1) {
    return 1;
  }
  return 0;
}

void
t8dg_quadrature_get_face_vertex (double vertex[3], const t8dg_quadrature_t * quadrature, const int iface, const t8dg_quad_idx_t ifacequad)
{
  T8DG_ASSERT (quadrature != NULL);
  T8DG_CHECK_ABORT (t8dg_quadrature_get_type (quadrature) == T8DG_QUAD_LGL, "Not yet implemented");
 /*TODO*/}

double
t8dg_quadrature_get_face_weight (const t8dg_quadrature_t * quadrature, const int iface, const t8dg_quad_idx_t ifacequad)
{
  /*TODO: tensor */
  return 1;
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
  t8dg_quad_idx_t     num_elem_quad, iquad;
  double              integral = 0;
  double              vertex[3];
  num_elem_quad = t8dg_quadrature_get_num_element_vertices (quadrature);
  for (iquad = 0; iquad < num_elem_quad; iquad++) {
    t8dg_quadrature_get_element_vertex (vertex, quadrature, iquad);
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
  int                 itensor;
  T8DG_ASSERT (pquadrature != NULL);
  t8dg_quadrature_t  *quadrature = *pquadrature;
  T8DG_ASSERT (quadrature != NULL);

  if (quadrature->num_tensor == 1) {
    sc_array_destroy (quadrature->weights);
    quadrature->weights = NULL;
    t8dg_vertexset_unref (&quadrature->vertexset);
    quadrature->vertexset = NULL;
  }
  else {
    for (itensor = 0; itensor < quadrature->num_tensor; itensor++) {
      t8dg_quadrature_unref (quadrature->tensor_quad + itensor);
    }
  }

  /*TODO: Face quadrature */

  T8_FREE (quadrature);
  *pquadrature = NULL;
}
