#include "t8dg.h"
#include "t8dg_functionbasis.h"
#include "t8dg_dof.h"
#include "t8dg_vertexset.h"
#include "t8dg_dmatrix.h"
#include "t8dg_refcount.h"
#include "t8dg_tensor.h"
#include "t8dg_sc_array.h"

#include <sc_containers.h>

/** The functionbasis provides the the directional derivative_matrix and interpolation/projection */
struct t8dg_functionbasis
{
  t8dg_refcount_t     rc;                               /**< Reference count */
  int                 number_of_dof;                    /**< Number of degrees of freedom*/
  t8dg_functionbasis_type_t type;                       /**< enum of possible types */

  t8_eclass_t         element_class;                    /**< Element class of the reference element */
  int                 embedded_dimension;
  int                 num_children;

  void               *data;
  t8dg_functionbasis_t **face_functionbasis;
  int                 num_face_functionbasis;
  sc_array_t        **lgl_faceidx_to_elemidx;

  int                 tensor;
};

typedef struct t8dg_functionbasis_0D_lagrange_data
{
  t8dg_vertexset_t   *vertexset;                        /**< Vertexset used for 1D Lagrange Basis*/
  t8dg_dmatrix_t     *identity;     /**< precalculated interpolation Matrix */
} t8dg_functionbasis_0D_lagrange_data_t;

typedef struct t8dg_functionbasis_1D_lagrange_data
{
  t8dg_vertexset_t   *vertexset;                        /**< Vertexset used for 1D Lagrange Basis*/
  sc_array_t         *barycentric_weights_1D;           /**<  barycentric weights used for 1D Lagrange Basis */
  t8dg_dmatrix_t     *derivative_matrix;                                /**< precalculated derivative Matrix */
  t8dg_dmatrix_t     *interpolate_to_child_matrix[2];     /**< precalculated interpolation Matrix */
} t8dg_functionbasis_1D_lagrange_data_t;

typedef struct t8dg_functionbasis_triangle_lagrange_data
{
  t8dg_vertexset_t   *vertexset;                        /**< Vertexset used for 1D Lagrange Basis*/
  t8dg_dmatrix_t     *derivative_matrix[2];                                /**< precalculated derivative Matrix */
  t8dg_dmatrix_t     *interpolate_to_child_matrix[4];     /**< precalculated interpolation Matrix */
} t8dg_functionbasis_triangle_lagrange_data_t;

typedef struct t8dg_functionbasis_tensor_data
{
  t8dg_functionbasis_t *tensor_first_functionbasis;                        /**<  Pointer to a tensorfunctionbasis, can be same*/
  t8dg_functionbasis_t *tensor_second_functionbasis;
} t8dg_functionbasis_tensor_data_t;

static void
t8dg_functionbasis_fill_barycentric_weights_1D (t8dg_functionbasis_t * functionbasis)
{
  T8DG_ASSERT (t8dg_functionbasis_is_valid (functionbasis));
  T8DG_ASSERT (t8dg_functionbasis_get_dim (functionbasis) == 1);
  T8DG_ASSERT (t8dg_functionbasis_is_lagrange (functionbasis));
  T8DG_ASSERT (functionbasis->data != NULL);
  t8dg_functionbasis_1D_lagrange_data_t *lagrange_data = (t8dg_functionbasis_1D_lagrange_data_t *) functionbasis->data;
  T8DG_ASSERT (lagrange_data->barycentric_weights_1D->elem_count == (size_t) t8dg_vertexset_get_num_vertices (lagrange_data->vertexset));
  int                 j, k;     /* Indices to enumerate the barycentric weights */
  double              dist;
  sc_array_t         *barycentric_weights = lagrange_data->barycentric_weights_1D;

  for (j = 0; j < functionbasis->number_of_dof; j++) {
    *(double *) sc_array_index_int (barycentric_weights, j) = 1;
  }

  for (j = 1; j < functionbasis->number_of_dof; j++) {
    for (k = 0; k < j; k++) {
      dist = t8dg_vertexset_get_first_coordinate (lagrange_data->vertexset, k)
        - t8dg_vertexset_get_first_coordinate (lagrange_data->vertexset, j);
      *(double *) sc_array_index_int (barycentric_weights, k) *= dist;
      *(double *) sc_array_index_int (barycentric_weights, j) *= (-dist);
    }
  }

  for (j = 0; j < functionbasis->number_of_dof; j++) {
    *(double *) sc_array_index_int (barycentric_weights, j) = 1. / *(double *) sc_array_index_int (barycentric_weights, j);
  }
}

static t8dg_dmatrix_t *
t8dg_functionbasis_1D_lagrange_interpolation_matrix (t8dg_functionbasis_t * functionbasis, t8dg_vertexset_t * interpolation_vertices)
{
  T8DG_ASSERT (t8dg_functionbasis_is_valid (functionbasis));
  T8DG_ASSERT (t8dg_functionbasis_is_lagrange (functionbasis));
  T8DG_CHECK_ABORT (t8dg_functionbasis_get_dim (functionbasis) == 1, "Not yet implemented");

  int                 number_interpolation_vertices, number_lagrange_vertices;
  int                 iinterpolation, ilagrange;
  int                 row_has_match;
  double              sum;
  double              entry;
  double              dist;
  double              x_lagrange, x_interpolation;
  t8dg_dmatrix_t     *interpolation_matrix;
  t8dg_functionbasis_1D_lagrange_data_t *lagrange_data;

  lagrange_data = (t8dg_functionbasis_1D_lagrange_data_t *) functionbasis->data;

  number_lagrange_vertices = t8dg_functionbasis_get_num_dof (functionbasis);
  number_interpolation_vertices = t8dg_vertexset_get_num_vertices (interpolation_vertices);

  interpolation_matrix = t8dg_dmatrix_new_zero (number_interpolation_vertices, number_lagrange_vertices);

  for (iinterpolation = 0; iinterpolation < number_interpolation_vertices; iinterpolation++) {
    row_has_match = 0;
    for (ilagrange = 0; ilagrange < number_lagrange_vertices; ilagrange++) {
      x_lagrange = t8dg_vertexset_get_first_coordinate (lagrange_data->vertexset, ilagrange);
      x_interpolation = t8dg_vertexset_get_first_coordinate (interpolation_vertices, iinterpolation);
      if (t8dg_almost_equal (x_lagrange, x_interpolation)) {
        row_has_match = 1;
        t8dg_dmatrix_set_at (interpolation_matrix, iinterpolation, ilagrange, 1);
      }
    }
    if (!row_has_match) {
      sum = 0;
      for (ilagrange = 0; ilagrange < number_lagrange_vertices; ilagrange++) {
        x_lagrange = t8dg_vertexset_get_first_coordinate (lagrange_data->vertexset, ilagrange);
        x_interpolation = t8dg_vertexset_get_first_coordinate (interpolation_vertices, iinterpolation);
        dist = x_interpolation - x_lagrange;
        entry = *(double *) sc_array_index_int (lagrange_data->barycentric_weights_1D, ilagrange);
        entry /= dist;
        t8dg_dmatrix_set_at (interpolation_matrix, iinterpolation, ilagrange, entry);
        sum += entry;
      }
      t8dg_dmatrix_scale_row (interpolation_matrix, iinterpolation, 1. / sum);
    }
  }

  return interpolation_matrix;
}

/* possible optimization: only calculate half the entries*/
static t8dg_dmatrix_t *
t8dg_functionbasis_1D_lagrange_derivative_matrix (t8dg_functionbasis_t * functionbasis)
{
  T8DG_ASSERT (t8dg_functionbasis_is_valid (functionbasis));
  T8DG_ASSERT (t8dg_functionbasis_is_lagrange (functionbasis));
  T8DG_ASSERT (t8dg_functionbasis_get_dim (functionbasis) == 1);

  int                 irow, icolumn;
  double              x_row, x_column, entry, bary_weight_column, bary_weight_row;
  int                 num_dof;
  t8dg_dmatrix_t     *derivative_matrix;
  t8dg_functionbasis_1D_lagrange_data_t *lagrange_data;

  lagrange_data = (t8dg_functionbasis_1D_lagrange_data_t *) functionbasis->data;
  num_dof = t8dg_functionbasis_get_num_dof (functionbasis);
  derivative_matrix = t8dg_dmatrix_new_zero (num_dof, num_dof);

  for (irow = 0; irow < num_dof; irow++) {
    for (icolumn = 0; icolumn < num_dof; icolumn++) {
      if (irow != icolumn) {
        bary_weight_row = *(double *) sc_array_index_int (lagrange_data->barycentric_weights_1D, irow);
        bary_weight_column = *(double *) sc_array_index_int (lagrange_data->barycentric_weights_1D, icolumn);
        x_row = t8dg_vertexset_get_first_coordinate (lagrange_data->vertexset, irow);
        x_column = t8dg_vertexset_get_first_coordinate (lagrange_data->vertexset, icolumn);
        entry = bary_weight_column / bary_weight_row / (x_row - x_column);
        t8dg_dmatrix_set_at (derivative_matrix, irow, icolumn, entry);
        /*Calculate diagonal entry using negative sum trick */
        t8dg_dmatrix_set_at (derivative_matrix, irow, irow, t8dg_dmatrix_at (derivative_matrix, irow, irow) - entry);
      }
    }
  }
  return derivative_matrix;
}

static t8dg_functionbasis_t *
t8dg_functionbasis_new_0D (t8dg_vertexset_t * vertex)
{
  t8dg_functionbasis_t *functionbasis;
  t8dg_functionbasis_0D_lagrange_data_t *lagrange_data;
  functionbasis = T8DG_ALLOC_ZERO (t8dg_functionbasis_t, 1);
  lagrange_data = T8DG_ALLOC_ZERO (t8dg_functionbasis_0D_lagrange_data_t, 1);

  lagrange_data->vertexset = vertex;
  lagrange_data->identity = t8dg_dmatrix_new (1, 1);
  t8dg_dmatrix_set_at (lagrange_data->identity, 0, 0, 1);

  functionbasis->data = lagrange_data;

  t8dg_refcount_init (&functionbasis->rc);

  functionbasis->element_class = T8_ECLASS_VERTEX;

  functionbasis->num_children = 1;
  functionbasis->type = T8DG_FB_LAGRANGE_LGL;

  functionbasis->number_of_dof = 1;
  functionbasis->embedded_dimension = t8dg_vertexset_get_embedded_dim (vertex);
  return functionbasis;
}

t8dg_functionbasis_t *
t8dg_functionbasis_new_1D_lagrange (t8dg_vertexset_t * vertexset, int create_face_functionbasis)
{
  T8DG_ASSERT (vertexset != NULL);
  T8DG_ASSERT (t8dg_vertexset_get_eclass (vertexset) == T8_ECLASS_LINE);

  t8dg_functionbasis_t *functionbasis;
  t8dg_functionbasis_1D_lagrange_data_t *lagrange_data;

  functionbasis = T8DG_ALLOC_ZERO (t8dg_functionbasis_t, 1);

  t8dg_refcount_init (&functionbasis->rc);
  functionbasis->element_class = T8_ECLASS_LINE;
  functionbasis->embedded_dimension = t8dg_vertexset_get_embedded_dim (vertexset);
  functionbasis->number_of_dof = t8dg_vertexset_get_num_vertices (vertexset);

  if (t8dg_vertexset_get_type (vertexset) == T8DG_VERT_LGL) {
    functionbasis->type = T8DG_FB_LAGRANGE_LGL;
  }
  else {
    T8DG_ABORT ("Not yet implemented");
  }

  lagrange_data = T8DG_ALLOC_ZERO (t8dg_functionbasis_1D_lagrange_data_t, 1);
  functionbasis->data = lagrange_data;
  t8dg_vertexset_ref (vertexset);
  lagrange_data->vertexset = vertexset;

  lagrange_data->barycentric_weights_1D = sc_array_new_count (sizeof (double), functionbasis->number_of_dof);
  t8dg_functionbasis_fill_barycentric_weights_1D (functionbasis);

  functionbasis->num_children = 2;
  t8dg_vertexset_t   *left_vertexset;
  t8dg_vertexset_t   *right_vertexset;

  left_vertexset = t8dg_vertexset_new_childvertexset (vertexset, 0);
  right_vertexset = t8dg_vertexset_new_childvertexset (vertexset, 1);

  lagrange_data->interpolate_to_child_matrix[0] = t8dg_functionbasis_1D_lagrange_interpolation_matrix (functionbasis, left_vertexset);
  lagrange_data->interpolate_to_child_matrix[1] = t8dg_functionbasis_1D_lagrange_interpolation_matrix (functionbasis, right_vertexset);

  t8dg_vertexset_destroy (&left_vertexset);
  t8dg_vertexset_destroy (&right_vertexset);

  lagrange_data->derivative_matrix = t8dg_functionbasis_1D_lagrange_derivative_matrix (functionbasis);

  if (create_face_functionbasis) {
    functionbasis->num_face_functionbasis = 2;
    functionbasis->face_functionbasis = T8DG_ALLOC (t8dg_functionbasis_t *, functionbasis->num_face_functionbasis);
    functionbasis->lgl_faceidx_to_elemidx = T8DG_ALLOC (sc_array_t *, functionbasis->num_face_functionbasis);
    if (t8dg_vertexset_get_type (vertexset) == T8DG_VERT_LGL) {
      t8dg_vertexset_t   *left_face_vertexset;
      t8dg_vertexset_t   *right_face_vertexset;
      left_face_vertexset = t8dg_vertexset_new_lgl_facevertexset (vertexset, 0);
      right_face_vertexset = t8dg_vertexset_new_lgl_facevertexset (vertexset, 1);
      functionbasis->face_functionbasis[0] = t8dg_functionbasis_new_0D (left_face_vertexset);
      functionbasis->face_functionbasis[1] = t8dg_functionbasis_new_0D (right_face_vertexset);

      functionbasis->lgl_faceidx_to_elemidx[0] = sc_array_new_count (sizeof (int), 1);
      functionbasis->lgl_faceidx_to_elemidx[1] = sc_array_new_count (sizeof (int), 1);
      *(int *) sc_array_index_int (functionbasis->lgl_faceidx_to_elemidx[0], 0) = 0;
      *(int *) sc_array_index_int (functionbasis->lgl_faceidx_to_elemidx[1], 0) = functionbasis->number_of_dof - 1;
    }
    else {
      T8DG_ABORT ("Not yet implemented");
    }
  }
  T8DG_ASSERT (t8dg_functionbasis_is_valid (functionbasis));
  return functionbasis;
}

t8dg_functionbasis_t *
t8dg_functionbasis_new_tensor (t8dg_functionbasis_t * tensor_first_functionbasis, t8dg_functionbasis_t * tensor_second_functionbasis,
                               int create_face_functionbasis)
{
  t8dg_functionbasis_t *functionbasis_tensor;
  t8dg_functionbasis_tensor_data_t *tensor_data;

  functionbasis_tensor = T8DG_ALLOC_ZERO (t8dg_functionbasis_t, 1);

  tensor_data = T8DG_ALLOC_ZERO (t8dg_functionbasis_tensor_data_t, 1);
  tensor_data->tensor_first_functionbasis = tensor_first_functionbasis;
  tensor_data->tensor_second_functionbasis = tensor_second_functionbasis;
  functionbasis_tensor->data = tensor_data;

  t8dg_refcount_init (&functionbasis_tensor->rc);

  functionbasis_tensor->element_class =
    t8dg_tensor_eclass (tensor_first_functionbasis->element_class, tensor_second_functionbasis->element_class);

  functionbasis_tensor->tensor = 1;
  functionbasis_tensor->number_of_dof = tensor_first_functionbasis->number_of_dof * tensor_second_functionbasis->number_of_dof;
  functionbasis_tensor->embedded_dimension =
    tensor_first_functionbasis->embedded_dimension + tensor_second_functionbasis->embedded_dimension;
  functionbasis_tensor->num_children = tensor_first_functionbasis->num_children * tensor_second_functionbasis->num_children;

  if (tensor_first_functionbasis->type == T8DG_FB_LAGRANGE_LGL && tensor_second_functionbasis->type == T8DG_FB_LAGRANGE_LGL) {
    functionbasis_tensor->type = T8DG_FB_LAGRANGE_LGL;
  }
  t8dg_functionbasis_ref (tensor_first_functionbasis);
  t8dg_functionbasis_ref (tensor_second_functionbasis);

  if (create_face_functionbasis) {
    int                 iface, ifacedof, num_face_dof, tensor1_num_face_dof, tensor1_num_dof, tensor1_faceidx_lookup;
    T8DG_ASSERT (tensor_first_functionbasis->num_face_functionbasis > 0 && tensor_second_functionbasis->num_face_functionbasis > 0);
    functionbasis_tensor->num_face_functionbasis =
      tensor_first_functionbasis->num_face_functionbasis + tensor_second_functionbasis->num_face_functionbasis;

    functionbasis_tensor->face_functionbasis = T8DG_ALLOC (t8dg_functionbasis_t *, functionbasis_tensor->num_face_functionbasis);
    functionbasis_tensor->lgl_faceidx_to_elemidx = T8DG_ALLOC (sc_array_t *, functionbasis_tensor->num_face_functionbasis);

    tensor1_num_dof = tensor_first_functionbasis->number_of_dof;

    /*Create faces corresponding to the first tensor fb */
    for (iface = 0; iface < tensor_first_functionbasis->num_face_functionbasis; iface++) {
      functionbasis_tensor->face_functionbasis[iface] =
        t8dg_functionbasis_new_tensor (tensor_first_functionbasis->face_functionbasis[iface], tensor_second_functionbasis, 0);

      num_face_dof = functionbasis_tensor->face_functionbasis[iface]->number_of_dof;
      functionbasis_tensor->lgl_faceidx_to_elemidx[iface] = sc_array_new_count (sizeof (int), num_face_dof);

      tensor1_num_face_dof = tensor_first_functionbasis->face_functionbasis[iface]->number_of_dof;

      for (ifacedof = 0; ifacedof < num_face_dof; ifacedof++) {
        tensor1_faceidx_lookup =
          *(int *) sc_array_index_int (tensor_first_functionbasis->lgl_faceidx_to_elemidx[iface], ifacedof % tensor1_num_face_dof);

        *(int *) sc_array_index_int (functionbasis_tensor->lgl_faceidx_to_elemidx[iface], ifacedof) = (ifacedof / tensor1_num_face_dof) * tensor1_num_dof + tensor1_faceidx_lookup;     /*TODO: Check */
      }
    }
    /*Create faces corresponding to the second tensor fb */
    for (iface = 0; iface < tensor_second_functionbasis->num_face_functionbasis; iface++) {
      functionbasis_tensor->face_functionbasis[iface + tensor_first_functionbasis->num_face_functionbasis] =
        t8dg_functionbasis_new_tensor (tensor_first_functionbasis, tensor_second_functionbasis->face_functionbasis[iface], 0);

      num_face_dof = functionbasis_tensor->face_functionbasis[iface + tensor_first_functionbasis->num_face_functionbasis]->number_of_dof;
      functionbasis_tensor->lgl_faceidx_to_elemidx[iface + tensor_first_functionbasis->num_face_functionbasis] =
        sc_array_new_count (sizeof (int), num_face_dof);

      tensor1_num_face_dof = tensor_first_functionbasis->face_functionbasis[iface]->number_of_dof;

      for (ifacedof = 0; ifacedof < num_face_dof; ifacedof++) {
        tensor1_faceidx_lookup =
          *(int *) sc_array_index_int (tensor_first_functionbasis->lgl_faceidx_to_elemidx[iface], ifacedof / tensor1_num_dof);

        *(int *) sc_array_index_int (functionbasis_tensor->lgl_faceidx_to_elemidx[iface + tensor_first_functionbasis->num_face_functionbasis], ifacedof) = (ifacedof % tensor1_num_dof) + tensor1_num_dof * tensor1_faceidx_lookup;     /*TODO: Check */
      }

    }
  }
  else {
    functionbasis_tensor->num_face_functionbasis = 0;
    functionbasis_tensor->face_functionbasis = NULL;
  }

  return functionbasis_tensor;
}

t8dg_functionbasis_t *
t8dg_functionbasis_new_hypercube_lagrange (int dim, t8dg_vertexset_t * vertexset, int create_face_functionbasis)
{
  T8DG_ASSERT (vertexset != NULL);
  T8DG_ASSERT (t8dg_vertexset_get_eclass (vertexset) == T8_ECLASS_LINE);
  T8DG_ASSERT (dim >= 1 && dim <= 3);

  t8dg_functionbasis_t *functionbasis1D;
  t8dg_functionbasis_t *functionbasis2D;
  t8dg_functionbasis_t *functionbasis3D;

  functionbasis1D = t8dg_functionbasis_new_1D_lagrange (vertexset, create_face_functionbasis);
  if (dim == 1) {
    return functionbasis1D;
  }
  functionbasis2D = t8dg_functionbasis_new_tensor (functionbasis1D, functionbasis1D, create_face_functionbasis);
  if (dim == 2) {
    t8dg_functionbasis_unref (&functionbasis1D);
    return functionbasis2D;
  }
  functionbasis3D = t8dg_functionbasis_new_tensor (functionbasis2D, functionbasis1D, create_face_functionbasis);
  t8dg_functionbasis_unref (&functionbasis2D);
  t8dg_functionbasis_unref (&functionbasis1D);
  return functionbasis3D;
}

static t8dg_vertexset_t *
t8dg_functionbasis_get_lagrange_vertexset (const t8dg_functionbasis_t * functionbasis)
{
  T8DG_ASSERT (t8dg_functionbasis_is_valid (functionbasis));
  T8DG_ASSERT (t8dg_functionbasis_is_lagrange (functionbasis));
  switch (functionbasis->element_class) {
  case T8_ECLASS_VERTEX:
    {
      t8dg_functionbasis_0D_lagrange_data_t *lagrange0D_data;
      lagrange0D_data = (t8dg_functionbasis_0D_lagrange_data_t *) functionbasis->data;
      return lagrange0D_data->vertexset;
      break;
    }
  case T8_ECLASS_LINE:
    {
      t8dg_functionbasis_1D_lagrange_data_t *lagrange1D_data;
      lagrange1D_data = (t8dg_functionbasis_1D_lagrange_data_t *) functionbasis->data;
      return lagrange1D_data->vertexset;
      break;
    }
  case T8_ECLASS_TRIANGLE:
    {
      t8dg_functionbasis_triangle_lagrange_data_t *lagrange_triangle_data;
      lagrange_triangle_data = (t8dg_functionbasis_triangle_lagrange_data_t *) functionbasis->data;
      return lagrange_triangle_data->vertexset;
      break;
    }

  default:
    T8DG_ABORT ("Not yet implemented");
    break;
  }
}

static void
t8dg_functionbasis_fill_lagrange_vertex (const t8dg_functionbasis_t * functionbasis, const int idof, const int startdim, double vertex[3])
{
  T8DG_ASSERT (idof >= 0 && idof < t8dg_functionbasis_get_num_dof (functionbasis));
  if (functionbasis->tensor) {
    t8dg_functionbasis_tensor_data_t *tensor_data = (t8dg_functionbasis_tensor_data_t *) functionbasis->data;
    t8dg_functionbasis_fill_lagrange_vertex (tensor_data->tensor_first_functionbasis,
                                             idof % t8dg_functionbasis_get_num_dof (tensor_data->tensor_first_functionbasis), startdim,
                                             vertex);
    t8dg_functionbasis_fill_lagrange_vertex (tensor_data->tensor_second_functionbasis,
                                             idof / t8dg_functionbasis_get_num_dof (tensor_data->tensor_first_functionbasis),
                                             startdim + tensor_data->tensor_first_functionbasis->embedded_dimension, vertex);
  }
  else {
    t8dg_vertexset_fill_vertex3D (t8dg_functionbasis_get_lagrange_vertexset (functionbasis), idof, startdim, vertex);
  }
}

void
t8dg_functionbasis_get_lagrange_vertex (const t8dg_functionbasis_t * functionbasis, const int idof, double vertex[3])
{
  T8DG_ASSERT (t8dg_functionbasis_is_valid (functionbasis));
  T8DG_ASSERT (t8dg_functionbasis_is_lagrange (functionbasis));
  T8DG_ASSERT (idof >= 0 && idof < t8dg_functionbasis_get_num_dof (functionbasis));
  int                 startdim = 0;

  vertex[0] = vertex[1] = vertex[2] = 0;

  t8dg_functionbasis_fill_lagrange_vertex (functionbasis, idof, startdim, vertex);
}

void
t8dg_functionbasis_interpolate_scalar_fn (const t8dg_functionbasis_t * functionbasis,
                                          t8dg_scalar_function_3d_fn function, void *scalar_fn_data, t8dg_element_dof_values_t * dof_values)
{
  T8DG_ASSERT (t8dg_functionbasis_is_valid (functionbasis));
  T8DG_ASSERT ((size_t) functionbasis->number_of_dof == dof_values->elem_count);
  int                 idof;
  double              vertex[3];
  if (t8dg_functionbasis_is_lagrange (functionbasis)) {
    for (idof = 0; idof < functionbasis->number_of_dof; idof++) {
      t8dg_functionbasis_get_lagrange_vertex (functionbasis, idof, vertex);
      t8dg_element_dof_values_set_value (dof_values, idof, function (vertex, scalar_fn_data));
    }
  }
  else {
    T8DG_ABORT ("Not yet implemented");
  }
}

t8dg_dmatrix_t     *
t8dg_functionbasis_get_lagrange_derivative_matrix (const t8dg_functionbasis_t * functionbasis, const int direction_idx)
{
  switch (functionbasis->element_class) {
  case T8_ECLASS_VERTEX:
    {
      T8DG_ABORT ("Is this really needed?");
      T8DG_ASSERT (direction_idx == 0);
      t8dg_functionbasis_0D_lagrange_data_t *lagrange0d_data;
      lagrange0d_data = (t8dg_functionbasis_0D_lagrange_data_t *) functionbasis->data;
      return lagrange0d_data->identity;
      break;
    }
  case T8_ECLASS_LINE:
    {
      T8DG_ASSERT (direction_idx == 0);
      t8dg_functionbasis_1D_lagrange_data_t *lagrange1d_data;
      lagrange1d_data = (t8dg_functionbasis_1D_lagrange_data_t *) functionbasis->data;
      return lagrange1d_data->derivative_matrix;
      break;
    }
  case T8_ECLASS_TRIANGLE:
    {
      t8dg_functionbasis_triangle_lagrange_data_t *lagrange_triangle_data;
      lagrange_triangle_data = (t8dg_functionbasis_triangle_lagrange_data_t *) functionbasis->data;
      return lagrange_triangle_data->derivative_matrix[direction_idx];
      break;
    }
  default:
    T8DG_ABORT ("Not yet implemented");
    break;
  }

}

/*TODO: make transpose char*/
void
t8dg_functionbasis_apply_derivative_matrix_transpose (const t8dg_functionbasis_t * functionbasis, int direction_idx,
                                                      t8dg_element_dof_values_t * derivative_dof_values,
                                                      t8dg_element_dof_values_t * dof_values)
{
  SC_CHECK_ABORT (t8dg_functionbasis_is_lagrange (functionbasis), "Not yet implemented");

  T8DG_ASSERT (direction_idx >= 0 && direction_idx < t8dg_functionbasis_get_dim (functionbasis));

  int                 num_vectors, ivector;
  int                 vector_length;
  int                 stride;
  sc_array_t         *vector_derivative_dof_values;
  sc_array_t         *vector_dof_values;

  if (functionbasis->tensor) {
    t8dg_functionbasis_tensor_data_t *tensor_data;
    t8dg_functionbasis_t *functionbasis_apply;
    tensor_data = (t8dg_functionbasis_tensor_data_t *) functionbasis->data;
    if (direction_idx < t8dg_functionbasis_get_dim (tensor_data->tensor_first_functionbasis)) {
      num_vectors = t8dg_functionbasis_get_num_dof (tensor_data->tensor_second_functionbasis);
      vector_length = t8dg_functionbasis_get_num_dof (tensor_data->tensor_first_functionbasis);
      stride = 1;
      functionbasis_apply = tensor_data->tensor_first_functionbasis;
    }
    else {
      num_vectors = t8dg_functionbasis_get_num_dof (tensor_data->tensor_first_functionbasis);
      vector_length = t8dg_functionbasis_get_num_dof (tensor_data->tensor_second_functionbasis);
      stride = num_vectors;
      direction_idx -= t8dg_functionbasis_get_dim (tensor_data->tensor_first_functionbasis);    /*TODO: or embedded dim? */
      functionbasis_apply = tensor_data->tensor_second_functionbasis;
    }
    vector_derivative_dof_values = sc_array_new_count (sizeof (double), vector_length);
    vector_dof_values = sc_array_new_count (sizeof (double), vector_length);
    for (ivector = 0; ivector < num_vectors; ivector++) {
      t8dg_tensor_array_extract_vector (derivative_dof_values, ivector, stride, vector_derivative_dof_values);
      t8dg_tensor_array_extract_vector (dof_values, ivector, stride, vector_dof_values);

      t8dg_functionbasis_apply_derivative_matrix_transpose (functionbasis_apply, direction_idx, vector_derivative_dof_values,
                                                            vector_dof_values);

      t8dg_tensor_array_inject_vector (vector_dof_values, ivector, stride, dof_values);
    }
    sc_array_destroy (vector_derivative_dof_values);
    sc_array_destroy (vector_dof_values);
  }
  else {
    t8dg_dmatrix_transpose_mult_sc_array (t8dg_functionbasis_get_lagrange_derivative_matrix (functionbasis, direction_idx),
                                          derivative_dof_values, dof_values);
  }
}

/*TODO: make transpose char*/
void
t8dg_functionbasis_apply_derivative_matrix (const t8dg_functionbasis_t * functionbasis, int direction_idx,
                                            t8dg_element_dof_values_t * dof_values, t8dg_element_dof_values_t * derivative_dof_values)
{
  SC_CHECK_ABORT (t8dg_functionbasis_is_lagrange (functionbasis), "Not yet implemented");

  int                 num_vectors, ivector;
  int                 vector_length;
  int                 stride;
  sc_array_t         *vector_derivative_dof_values;
  sc_array_t         *vector_dof_values;

  if (functionbasis->tensor) {
    t8dg_functionbasis_tensor_data_t *tensor_data;
    t8dg_functionbasis_t *functionbasis_apply;
    tensor_data = (t8dg_functionbasis_tensor_data_t *) functionbasis->data;
    if (direction_idx < t8dg_functionbasis_get_dim (tensor_data->tensor_first_functionbasis)) {
      num_vectors = t8dg_functionbasis_get_num_dof (tensor_data->tensor_second_functionbasis);
      vector_length = t8dg_functionbasis_get_num_dof (tensor_data->tensor_first_functionbasis);
      stride = 1;
      functionbasis_apply = tensor_data->tensor_first_functionbasis;
    }
    else {
      num_vectors = t8dg_functionbasis_get_num_dof (tensor_data->tensor_first_functionbasis);
      vector_length = t8dg_functionbasis_get_num_dof (tensor_data->tensor_second_functionbasis);
      stride = num_vectors;
      direction_idx -= t8dg_functionbasis_get_dim (tensor_data->tensor_first_functionbasis);    /*TODO: or embedded dim? */
      functionbasis_apply = tensor_data->tensor_second_functionbasis;
    }
    vector_derivative_dof_values = sc_array_new_count (sizeof (double), vector_length);
    vector_dof_values = sc_array_new_count (sizeof (double), vector_length);
    for (ivector = 0; ivector < num_vectors; ivector++) {
      t8dg_tensor_array_extract_vector (derivative_dof_values, ivector, stride, vector_derivative_dof_values);
      t8dg_tensor_array_extract_vector (dof_values, ivector, stride, vector_dof_values);

      t8dg_functionbasis_apply_derivative_matrix (functionbasis_apply, direction_idx, vector_dof_values, vector_derivative_dof_values);

      t8dg_tensor_array_inject_vector (vector_derivative_dof_values, ivector, stride, derivative_dof_values);
    }
    sc_array_destroy (vector_derivative_dof_values);
    sc_array_destroy (vector_dof_values);
  }
  else {
    t8dg_dmatrix_mult_sc_array (t8dg_functionbasis_get_lagrange_derivative_matrix (functionbasis, direction_idx), dof_values,
                                derivative_dof_values);
  }
}

void
t8dg_functionbasis_apply_child_interpolation_matrix (const t8dg_functionbasis_t * functionbasis, const int ichild,
                                                     t8dg_element_dof_values_t * element_dof, t8dg_element_dof_values_t * child_dof)
{
  int                 ivector, num_vectors, vector_length, stride;
  sc_array_t         *vector_element_dof;
  sc_array_t         *vector_child_dof;
  if (functionbasis->tensor) {
    t8dg_functionbasis_tensor_data_t *tensor_data;
    tensor_data = (t8dg_functionbasis_tensor_data_t *) functionbasis->data;

    /*Interpolate first tensor */
    num_vectors = t8dg_functionbasis_get_num_dof (tensor_data->tensor_second_functionbasis);
    vector_length = t8dg_functionbasis_get_num_dof (tensor_data->tensor_first_functionbasis);
    stride = 1;

    vector_element_dof = sc_array_new_count (sizeof (double), vector_length);
    vector_child_dof = sc_array_new_count (sizeof (double), vector_length);

    for (ivector = 0; ivector < num_vectors; ivector++) {
      t8dg_tensor_array_extract_vector (element_dof, ivector, stride, vector_element_dof);

      T8DG_ASSERT (tensor_data->tensor_first_functionbasis->num_children);
      t8dg_functionbasis_apply_child_interpolation_matrix (tensor_data->tensor_first_functionbasis,
                                                           ichild % tensor_data->tensor_first_functionbasis->num_children,
                                                           vector_element_dof, vector_child_dof);

      t8dg_tensor_array_inject_vector (vector_child_dof, ivector, stride, child_dof);
    }
    sc_array_destroy (vector_element_dof);
    sc_array_destroy (vector_child_dof);

    num_vectors = t8dg_functionbasis_get_num_dof (tensor_data->tensor_first_functionbasis);
    vector_length = t8dg_functionbasis_get_num_dof (tensor_data->tensor_second_functionbasis);
    stride = num_vectors;

    vector_element_dof = sc_array_new_count (sizeof (double), vector_length);
    vector_child_dof = sc_array_new_count (sizeof (double), vector_length);

    for (ivector = 0; ivector < num_vectors; ivector++) {
      /*Extract from the interpolation in first direction */
      t8dg_tensor_array_extract_vector (child_dof, ivector, stride, vector_element_dof);

      t8dg_functionbasis_apply_child_interpolation_matrix (tensor_data->tensor_second_functionbasis,
                                                           ichild / tensor_data->tensor_first_functionbasis->num_children,
                                                           vector_element_dof, vector_child_dof);

      t8dg_tensor_array_inject_vector (vector_child_dof, ivector, stride, child_dof);
    }
    sc_array_destroy (vector_element_dof);
    sc_array_destroy (vector_child_dof);
  }
  else {
    t8dg_dmatrix_t     *interpolation_matrix;
    interpolation_matrix = t8dg_functionbasis_get_lagrange_child_interpolation_matrix (functionbasis, ichild);
    t8dg_dmatrix_mult_sc_array (interpolation_matrix, element_dof, child_dof);
  }
}

void
t8dg_functionbasis_apply_child_interpolation_matrix_transpose (const t8dg_functionbasis_t * functionbasis, const int ichild,
                                                               t8dg_element_dof_values_t * child_dof,
                                                               t8dg_element_dof_values_t * element_dof)
{
  int                 ivector, num_vectors, vector_length, stride;
  sc_array_t         *vector_element_dof;
  sc_array_t         *vector_child_dof;
  if (functionbasis->tensor) {
    t8dg_functionbasis_tensor_data_t *tensor_data;
    tensor_data = (t8dg_functionbasis_tensor_data_t *) functionbasis->data;

    /*Interpolate first tensor */
    num_vectors = t8dg_functionbasis_get_num_dof (tensor_data->tensor_second_functionbasis);
    vector_length = t8dg_functionbasis_get_num_dof (tensor_data->tensor_first_functionbasis);
    stride = 1;

    vector_element_dof = sc_array_new_count (sizeof (double), vector_length);
    vector_child_dof = sc_array_new_count (sizeof (double), vector_length);

    for (ivector = 0; ivector < num_vectors; ivector++) {
      t8dg_tensor_array_extract_vector (child_dof, ivector, stride, vector_child_dof);

      t8dg_functionbasis_apply_child_interpolation_matrix_transpose (tensor_data->tensor_first_functionbasis,
                                                                     ichild % tensor_data->tensor_first_functionbasis->num_children,
                                                                     vector_child_dof, vector_element_dof);

      t8dg_tensor_array_inject_vector (vector_element_dof, ivector, stride, element_dof);
    }
    sc_array_destroy (vector_element_dof);
    sc_array_destroy (vector_child_dof);

    num_vectors = t8dg_functionbasis_get_num_dof (tensor_data->tensor_first_functionbasis);
    vector_length = t8dg_functionbasis_get_num_dof (tensor_data->tensor_second_functionbasis);
    stride = num_vectors;

    vector_element_dof = sc_array_new_count (sizeof (double), vector_length);
    vector_child_dof = sc_array_new_count (sizeof (double), vector_length);

    for (ivector = 0; ivector < num_vectors; ivector++) {
      /*Extract from the interpolation in first direction */
      t8dg_tensor_array_extract_vector (element_dof, ivector, stride, vector_child_dof);

      t8dg_functionbasis_apply_child_interpolation_matrix_transpose (tensor_data->tensor_second_functionbasis,
                                                                     ichild / tensor_data->tensor_first_functionbasis->num_children,
                                                                     vector_child_dof, vector_element_dof);

      t8dg_tensor_array_inject_vector (vector_element_dof, ivector, stride, element_dof);
    }
    sc_array_destroy (vector_element_dof);
    sc_array_destroy (vector_child_dof);
  }
  else {
    t8dg_dmatrix_t     *interpolation_matrix;
    interpolation_matrix = t8dg_functionbasis_get_lagrange_child_interpolation_matrix (functionbasis, ichild);
    t8dg_dmatrix_transpose_mult_sc_array (interpolation_matrix, child_dof, element_dof);
  }
}

static int
t8dg_functionbasis_lgl_facedof_idx_lookup (const t8dg_functionbasis_t * functionbasis, const int iface, const int ifacedof)
{
  T8DG_ASSERT (functionbasis->type == T8DG_FB_LAGRANGE_LGL);
  T8DG_ASSERT (iface >= 0 && iface < functionbasis->num_face_functionbasis);
  T8DG_ASSERT (ifacedof >= 0 && ifacedof < functionbasis->face_functionbasis[iface]->number_of_dof);

  return *(int *) sc_array_index_int (functionbasis->lgl_faceidx_to_elemidx[iface], ifacedof);
}

void
t8dg_functionbasis_transform_element_dof_to_face_dof (const t8dg_functionbasis_t * functionbasis, const int iface,
                                                      t8dg_element_dof_values_t * element_dof_array,
                                                      t8dg_face_dof_values_t * face_dof_array)
{
  T8DG_CHECK_ABORT (functionbasis->type == T8DG_FB_LAGRANGE_LGL, "Not implemented");

  int                 idof, ifacedof, num_face_dof;

  num_face_dof = t8dg_functionbasis_get_num_face_dof (functionbasis, iface);

  for (ifacedof = 0; ifacedof < num_face_dof; ifacedof++) {
    idof = t8dg_functionbasis_lgl_facedof_idx_lookup (functionbasis, iface, ifacedof);
    *(double *) sc_array_index_int (face_dof_array, ifacedof) = *(double *) sc_array_index_int (element_dof_array, idof);
  }
}

void
t8dg_functionbasis_transform_face_dof_to_element_dof (const t8dg_functionbasis_t * functionbasis, const int iface,
                                                      t8dg_face_dof_values_t * face_dof_array,
                                                      t8dg_element_dof_values_t * element_dof_array)
{
  T8DG_CHECK_ABORT (functionbasis->type == T8DG_FB_LAGRANGE_LGL, "Not implemented");

  int                 idof, ifacedof, num_face_dof;

  t8dg_element_dof_values_set_zero (element_dof_array);

  num_face_dof = t8dg_functionbasis_get_num_face_dof (functionbasis, iface);

  for (ifacedof = 0; ifacedof < num_face_dof; ifacedof++) {
    idof = t8dg_functionbasis_lgl_facedof_idx_lookup (functionbasis, iface, ifacedof);
    *(double *) sc_array_index_int (element_dof_array, idof) = *(double *) sc_array_index_int (face_dof_array, ifacedof);
  }
}

int
t8dg_functionbasis_is_valid (const t8dg_functionbasis_t * functionbasis)
{
  return (functionbasis != NULL && t8dg_refcount_is_active (&functionbasis->rc));
  /* ADD more checks */
}

int
t8dg_functionbasis_is_lagrange (const t8dg_functionbasis_t * functionbasis)
{
  T8DG_ASSERT (t8dg_functionbasis_is_valid (functionbasis));
  return functionbasis->type == T8DG_FB_LAGRANGE_LGL;
}

int
t8dg_functionbasis_is_tensor (const t8dg_functionbasis_t * functionbasis)
{
  T8DG_ASSERT (t8dg_functionbasis_is_valid (functionbasis));
  return functionbasis->tensor;
}

int
t8dg_functionbasis_get_num_dof (const t8dg_functionbasis_t * functionbasis)
{
  T8DG_ASSERT (t8dg_functionbasis_is_valid (functionbasis));
  return functionbasis->number_of_dof;
}

t8dg_functionbasis_type_t
t8dg_functionbasis_get_type (const t8dg_functionbasis_t * functionbasis)
{
  T8DG_ASSERT (t8dg_functionbasis_is_valid (functionbasis));
  return functionbasis->type;
}

int
t8dg_functionbasis_get_dim (const t8dg_functionbasis_t * functionbasis)
{
  T8DG_ASSERT (t8dg_functionbasis_is_valid (functionbasis));
  return t8_eclass_to_dimension[functionbasis->element_class];
}

t8_eclass_t
t8dg_functionbasis_get_eclass (const t8dg_functionbasis_t * functionbasis)
{
  T8DG_ASSERT (t8dg_functionbasis_is_valid (functionbasis));
  return functionbasis->element_class;
}

int
t8dg_functionbasis_get_num_children (t8dg_functionbasis_t * functionbasis)
{
  T8DG_ASSERT (t8dg_functionbasis_is_valid (functionbasis));
  return 1 << t8dg_functionbasis_get_dim (functionbasis);
}

t8dg_functionbasis_t *
t8dg_functionbasis_get_face_functionbasis (t8dg_functionbasis_t * functionbasis, int iface)
{
  T8DG_ASSERT (t8dg_functionbasis_is_valid (functionbasis));
  T8DG_ASSERT (iface >= 0 && iface < functionbasis->num_face_functionbasis);
  return functionbasis->face_functionbasis[iface];
}

int
t8dg_functionbasis_get_num_face_functionbasis (const t8dg_functionbasis_t * functionbasis)
{
  T8DG_ASSERT (t8dg_functionbasis_is_valid (functionbasis));
  return functionbasis->num_face_functionbasis;
}

int
t8dg_functionbasis_get_num_face_children (const t8dg_functionbasis_t * functionbasis, int iface)
{
  T8DG_ASSERT (t8dg_functionbasis_is_valid (functionbasis));
  T8DG_ASSERT (iface >= 0 && iface < functionbasis->num_face_functionbasis);
  return t8dg_functionbasis_get_num_children (functionbasis->face_functionbasis[iface]);
}

int
t8dg_functionbasis_get_num_face_dof (const t8dg_functionbasis_t * functionbasis, const int iface)
{
  T8DG_ASSERT (t8dg_functionbasis_is_valid (functionbasis));
  T8DG_ASSERT (iface >= 0 && iface < functionbasis->num_face_functionbasis);
  return t8dg_functionbasis_get_num_dof (functionbasis->face_functionbasis[iface]);
}

void
t8dg_functionbasis_get_lagrange_face_vertex (const t8dg_functionbasis_t * functionbasis, const int iface, const int idof, double vertex[3])
{
  T8DG_ASSERT (t8dg_functionbasis_is_valid (functionbasis));
  T8DG_ASSERT (iface >= 0 && iface < functionbasis->num_face_functionbasis);
  t8dg_functionbasis_get_lagrange_vertex (functionbasis->face_functionbasis[iface], idof, vertex);
}

t8dg_dmatrix_t     *
t8dg_functionbasis_get_lagrange_child_interpolation_matrix (const t8dg_functionbasis_t * functionbasis, const int ichild)
{
  T8DG_ASSERT (ichild >= 0 && ichild < functionbasis->num_children);
  switch (functionbasis->element_class) {
  case T8_ECLASS_VERTEX:
    {
      t8dg_functionbasis_0D_lagrange_data_t *lagrange0d_data;
      lagrange0d_data = (t8dg_functionbasis_0D_lagrange_data_t *) functionbasis->data;
      return lagrange0d_data->identity;
      break;
    }
  case T8_ECLASS_LINE:
    {
      t8dg_functionbasis_1D_lagrange_data_t *lagrange1d_data;
      lagrange1d_data = (t8dg_functionbasis_1D_lagrange_data_t *) functionbasis->data;
      return lagrange1d_data->interpolate_to_child_matrix[ichild];
      break;
    }
  case T8_ECLASS_TRIANGLE:
    {
      t8dg_functionbasis_triangle_lagrange_data_t *lagrange_triangle_data;
      lagrange_triangle_data = (t8dg_functionbasis_triangle_lagrange_data_t *) functionbasis->data;
      return lagrange_triangle_data->interpolate_to_child_matrix[ichild];
      break;
    }
  default:
    T8DG_ABORT ("Not yet implemented");
    break;
  }
}

void
t8dg_functionbasis_ref (t8dg_functionbasis_t * functionbasis)
{
  T8DG_ASSERT (t8dg_functionbasis_is_valid (functionbasis));
  t8dg_refcount_ref (&functionbasis->rc);
}

static void
t8dg_functionbasis_reset (t8dg_functionbasis_t * functionbasis)
{
  T8DG_ASSERT (functionbasis != NULL);
  int                 iface;
  for (iface = 0; iface < functionbasis->num_face_functionbasis; iface++) {
    t8dg_functionbasis_destroy (&functionbasis->face_functionbasis[iface]);
    sc_array_destroy (functionbasis->lgl_faceidx_to_elemidx[iface]);
  }
  if (functionbasis->num_face_functionbasis > 0) {
    T8DG_FREE (functionbasis->face_functionbasis);
    T8DG_FREE (functionbasis->lgl_faceidx_to_elemidx);
  }

  if (t8dg_functionbasis_is_tensor (functionbasis)) {
    t8dg_functionbasis_tensor_data_t *tensor_data;
    tensor_data = (t8dg_functionbasis_tensor_data_t *) functionbasis->data;

    t8dg_functionbasis_unref (&tensor_data->tensor_first_functionbasis);
    t8dg_functionbasis_unref (&tensor_data->tensor_second_functionbasis);
    T8DG_FREE (tensor_data);
  }
  else if (t8dg_functionbasis_is_lagrange (functionbasis)) {
    switch (functionbasis->element_class) {
    case T8_ECLASS_VERTEX:
      {
        t8dg_functionbasis_0D_lagrange_data_t *lagrange0d_data;
        lagrange0d_data = (t8dg_functionbasis_0D_lagrange_data_t *) functionbasis->data;
        t8dg_dmatrix_destroy (&lagrange0d_data->identity);
        t8dg_vertexset_unref (&lagrange0d_data->vertexset);
        T8DG_FREE (lagrange0d_data);
        break;
      }
    case T8_ECLASS_LINE:
      {
        t8dg_functionbasis_1D_lagrange_data_t *lagrange1d_data;
        lagrange1d_data = (t8dg_functionbasis_1D_lagrange_data_t *) functionbasis->data;
        sc_array_destroy (lagrange1d_data->barycentric_weights_1D);
        t8dg_dmatrix_destroy (&lagrange1d_data->derivative_matrix);
        t8dg_dmatrix_destroy (&lagrange1d_data->interpolate_to_child_matrix[0]);
        t8dg_dmatrix_destroy (&lagrange1d_data->interpolate_to_child_matrix[1]);
        t8dg_vertexset_unref (&lagrange1d_data->vertexset);
        T8DG_FREE (lagrange1d_data);
        break;
      }
    case T8_ECLASS_TRIANGLE:
      {
        t8dg_functionbasis_triangle_lagrange_data_t *lagrange_triangle_data;
        lagrange_triangle_data = (t8dg_functionbasis_triangle_lagrange_data_t *) functionbasis->data;
        t8dg_dmatrix_destroy (&lagrange_triangle_data->derivative_matrix[0]);
        t8dg_dmatrix_destroy (&lagrange_triangle_data->derivative_matrix[1]);
        t8dg_dmatrix_destroy (&lagrange_triangle_data->interpolate_to_child_matrix[0]);
        t8dg_dmatrix_destroy (&lagrange_triangle_data->interpolate_to_child_matrix[1]);
        t8dg_dmatrix_destroy (&lagrange_triangle_data->interpolate_to_child_matrix[2]);
        t8dg_dmatrix_destroy (&lagrange_triangle_data->interpolate_to_child_matrix[3]);
        t8dg_vertexset_unref (&lagrange_triangle_data->vertexset);
        T8DG_FREE (lagrange_triangle_data);
        break;
      }
    default:
      T8DG_ABORT ("Not yet implemented");
      break;
    }
  }
  else {
    T8DG_ABORT ("Not yet implemented");
  }
  t8dg_refcount_unref (&functionbasis->rc);
  T8DG_FREE (functionbasis);
}

void
t8dg_functionbasis_unref (t8dg_functionbasis_t ** pfunctionbasis)
{
  t8dg_functionbasis_t *functionbasis;

  T8DG_ASSERT (pfunctionbasis != NULL);
  functionbasis = *pfunctionbasis;
  T8DG_ASSERT (functionbasis != NULL);
  T8DG_ASSERT (t8dg_refcount_is_active (&functionbasis->rc));

  if (t8dg_refcount_is_last (&functionbasis->rc)) {
    t8dg_functionbasis_reset (functionbasis);
  }
  else {
    t8dg_refcount_unref (&functionbasis->rc);
  }
  *pfunctionbasis = NULL;
}

void
t8dg_functionbasis_destroy (t8dg_functionbasis_t ** pfunctionbasis)
{
  t8dg_functionbasis_t *functionbasis;

  T8DG_ASSERT (pfunctionbasis != NULL);
  functionbasis = *pfunctionbasis;

  T8DG_ASSERT (functionbasis != NULL);
  T8DG_ASSERT (t8dg_refcount_is_last (&functionbasis->rc));

  t8dg_functionbasis_reset (functionbasis);
  *pfunctionbasis = NULL;
}
