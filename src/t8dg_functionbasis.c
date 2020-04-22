#include "t8dg.h"
#include "t8dg_functionbasis.h"
#include "t8dg_vertexset.h"
#include <sc_containers.h>
#include "t8dg_dmatrix.h"
#include "t8dg_refcount.h"

/** The functionbasis provides the the directional derivative_matrix and interpolation/projection */
struct t8dg_functionbasis
{
  t8dg_refcount_t     rc;
  int                 number_of_dof;                    /**< Number of degrees of freedom*/
  t8dg_functionbasis_type_t type;

  t8_eclass_t         element_class;

  /** Vertexset used for 1D Lagrange Basis*/
  t8dg_vertexset_t   *vertexset;
  /** barycentric weights used for 1D Lagrange Basis */
  sc_array_t         *barycentric_weights_1D;

  /*precalculated values for Triangle/tet ? */

  /** Only explicitely saved for non tensor*/
  t8dg_dmatrix_t     *derivative_matrix;
  t8dg_dmatrix_t     *interpolate_to_child_matrix[MAX_SUBELEMENTS];

  /** Tensor Information */
  int                 num_tensor;
                        /**< if == 1, line, tri or tet (or pyramid) */
  t8dg_functionbasis_t *tensor_fb[DIM3];
  int                 tensor_num_dof[DIM3];
};

int
t8dg_functionbasis_is_lagrange (const t8dg_functionbasis_t * functionbasis)
{
  return functionbasis->type == T8DG_LAGRANGE_LGL_1D || functionbasis->type == T8DG_LAGRANGE_GL_1D;
}

static void
t8dg_functionbasis_fill_barycentric_weights_1D (sc_array_t * barycentric_weights, t8dg_vertexset_t * vertexset)
{
  T8DG_ASSERT (t8dg_vertexset_get_dim (vertexset) == 1);
  T8DG_ASSERT (barycentric_weights->elem_count == (size_t) t8dg_vertexset_get_num_vertices (vertexset));
  int                 poly_degree = barycentric_weights->elem_count - 1;
  int                 j, k;     /* Indices to enumerate the barycentric weights */
  double              dist;

  for (j = 0; j <= poly_degree; j++) {
    *(double *) sc_array_index_int (barycentric_weights, j) = 1;
  }

  for (j = 1; j <= poly_degree; j++) {
    for (k = 0; k < j; k++) {
      dist = t8dg_vertexset_get_first_coordinate (vertexset, k) - t8dg_vertexset_get_first_coordinate (vertexset, j);
      *(double *) sc_array_index_int (barycentric_weights, k) *= dist;
      *(double *) sc_array_index_int (barycentric_weights, j) *= (-dist);
    }
  }

  for (j = 0; j <= poly_degree; j++) {
    *(double *) sc_array_index_int (barycentric_weights, j) = 1. / *(double *) sc_array_index_int (barycentric_weights, j);
  }
}

t8dg_functionbasis_t *
t8dg_functionbasis_new_1D_Lagrange (t8dg_vertexset_t * vertexset)
{
  T8DG_ASSERT (vertexset != NULL);
  T8DG_ASSERT (t8dg_vertexset_get_eclass (vertexset) == T8_ECLASS_LINE);

  t8dg_functionbasis_t *functionbasis;
  functionbasis = T8DG_ALLOC_ZERO (t8dg_functionbasis_t, 1);

  t8dg_refcount_init (&functionbasis->rc);

  t8dg_vertexset_ref (vertexset);
  functionbasis->vertexset = vertexset;
  functionbasis->element_class = t8dg_vertexset_get_eclass (vertexset);
  functionbasis->num_tensor = 1;
  functionbasis->number_of_dof = t8dg_vertexset_get_num_vertices (vertexset);

  functionbasis->barycentric_weights_1D = sc_array_new_count (sizeof (double), functionbasis->number_of_dof);
  t8dg_functionbasis_fill_barycentric_weights_1D (functionbasis->barycentric_weights_1D, functionbasis->vertexset);

  t8dg_vertexset_t   *left_vertexset;
  t8dg_vertexset_t   *right_vertexset;

  left_vertexset = t8dg_vertexset_new_childvertexset_1D (vertexset, 0);
  right_vertexset = t8dg_vertexset_new_childvertexset_1D (vertexset, 1);

  functionbasis->interpolate_to_child_matrix[0] = t8dg_functionbasis_Lagrange_interpolation_matrix (functionbasis, left_vertexset);
  functionbasis->interpolate_to_child_matrix[1] = t8dg_functionbasis_Lagrange_interpolation_matrix (functionbasis, right_vertexset);

  t8dg_vertexset_destroy (&left_vertexset);
  t8dg_vertexset_destroy (&right_vertexset);

  functionbasis->derivative_matrix = t8dg_functionbasis_Lagrange_derivative_matrix (functionbasis);

  switch (t8dg_vertexset_get_type (vertexset)) {
  case (T8DG_VERT_LGL):
    functionbasis->type = T8DG_LAGRANGE_LGL_1D;
    break;
  case (T8DG_VERT_GL):
    functionbasis->type = T8DG_LAGRANGE_GL_1D;
    break;
  default:
    T8DG_ABORT ("Not yet implemented");
  }
  return functionbasis;
}

int
t8dg_functionbasis_get_num_dof (const t8dg_functionbasis_t * functionbasis)
{
  T8DG_ASSERT (functionbasis != NULL);
  return functionbasis->number_of_dof;
}

static void
t8dg_functionbasis_get_Lagrange_vertex (const t8dg_functionbasis_t * functionbasis, const int idof, double vertex[3])
{
  T8DG_ASSERT (functionbasis != NULL);
  T8DG_ASSERT (t8dg_functionbasis_is_lagrange (functionbasis));
  T8DG_ASSERT (idof >= 0 && idof < t8dg_functionbasis_get_num_dof (functionbasis));
  int                 itensor;
  int                 idoftensor[DIM3];
  int                 startdim = 0;

  vertex[0] = vertex[1] = vertex[2] = 0;

  if (functionbasis->num_tensor == 1) {
    t8dg_vertexset_fill_vertex3D (functionbasis->vertexset, idof, 0, vertex);
  }
  else {
    t8dg_transform_3tensoridx (idof, functionbasis->tensor_num_dof, idoftensor);
    for (itensor = 0; itensor < functionbasis->num_tensor; itensor++) {
      t8dg_vertexset_fill_vertex3D (functionbasis->tensor_fb[itensor]->vertexset, idoftensor[itensor], startdim, vertex);
      startdim += t8dg_functionbasis_get_dim (functionbasis->tensor_fb[itensor]);
    }
  }
}

t8dg_functionbasis_type_t
t8dg_functionbasis_get_type (const t8dg_functionbasis_t * functionbasis)
{
  return functionbasis->type;
}

int
t8dg_functionbasis_get_dim (const t8dg_functionbasis_t * functionbasis)
{
  return t8_eclass_to_dimension[functionbasis->element_class];
}

t8_eclass_t
t8dg_functionbasis_get_eclass (const t8dg_functionbasis_t * functionbasis)
{
  return functionbasis->element_class;
}

static double
t8dg_functionbasis_get_barycentric_weight (t8dg_functionbasis_t * functionbasis, int idof)
{
  T8DG_ASSERT (t8dg_functionbasis_is_lagrange (functionbasis));
  T8DG_ASSERT (t8dg_functionbasis_get_dim (functionbasis) == 1);
  return *(double *) sc_array_index_int (functionbasis->barycentric_weights_1D, idof);
}

t8dg_dmatrix_t     *
t8dg_functionbasis_Lagrange_interpolation_matrix (t8dg_functionbasis_t * functionbasis, t8dg_vertexset_t * interpolation_vertices)
{
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

  number_lagrange_vertices = t8dg_functionbasis_get_num_dof (functionbasis);
  number_interpolation_vertices = t8dg_vertexset_get_num_vertices (interpolation_vertices);

  interpolation_matrix = t8dg_dmatrix_new_zero (number_interpolation_vertices, number_lagrange_vertices);

  for (iinterpolation = 0; iinterpolation < number_interpolation_vertices; iinterpolation++) {
    row_has_match = 0;
    for (ilagrange = 0; ilagrange < number_lagrange_vertices; ilagrange++) {
      x_lagrange = t8dg_vertexset_get_first_coordinate (functionbasis->vertexset, ilagrange);
      x_interpolation = t8dg_vertexset_get_first_coordinate (interpolation_vertices, iinterpolation);
      if (t8dg_almost_equal (x_lagrange, x_interpolation)) {
        row_has_match = 1;
        t8dg_dmatrix_set_at (interpolation_matrix, iinterpolation, ilagrange, 1);
      }
    }
    if (!row_has_match) {
      sum = 0;
      for (ilagrange = 0; ilagrange < number_lagrange_vertices; ilagrange++) {
        x_lagrange = t8dg_vertexset_get_first_coordinate (functionbasis->vertexset, ilagrange);
        x_interpolation = t8dg_vertexset_get_first_coordinate (interpolation_vertices, iinterpolation);
        dist = x_interpolation - x_lagrange;
        entry = t8dg_functionbasis_get_barycentric_weight (functionbasis, ilagrange) / dist;
        t8dg_dmatrix_set_at (interpolation_matrix, iinterpolation, ilagrange, entry);
        sum += entry;
      }
      t8dg_dmatrix_scale_row (interpolation_matrix, iinterpolation, 1. / sum);
    }
  }

  return interpolation_matrix;
}

/* possible optimization: only calculate half the entries*/
t8dg_dmatrix_t     *
t8dg_functionbasis_Lagrange_derivative_matrix (t8dg_functionbasis_t * functionbasis)
{
  T8DG_ASSERT (t8dg_functionbasis_is_lagrange (functionbasis));
  T8DG_CHECK_ABORT (t8dg_functionbasis_get_dim (functionbasis) == 1, "Not yet implemented");

  int                 irow, icolumn;
  double              x_row, x_column, entry, bary_weight_column, bary_weight_row;
  int                 num_dof;
  t8dg_dmatrix_t     *derivative_matrix;

  num_dof = t8dg_functionbasis_get_num_dof (functionbasis);
  derivative_matrix = t8dg_dmatrix_new_zero (num_dof, num_dof);

  for (irow = 0; irow < num_dof; irow++) {
    for (icolumn = 0; icolumn < num_dof; icolumn++) {
      if (irow != icolumn) {
        bary_weight_row = *(double *) sc_array_index_int (functionbasis->barycentric_weights_1D, irow);
        bary_weight_column = *(double *) sc_array_index_int (functionbasis->barycentric_weights_1D, icolumn);
        x_row = t8dg_vertexset_get_first_coordinate (functionbasis->vertexset, irow);
        x_column = t8dg_vertexset_get_first_coordinate (functionbasis->vertexset, icolumn);
        entry = bary_weight_column / bary_weight_row / (x_row - x_column);
        t8dg_dmatrix_set_at (derivative_matrix, irow, icolumn, entry);
        /*Calculate diagonal entry using negative sum trick */
        t8dg_dmatrix_set_at (derivative_matrix, irow, irow, t8dg_dmatrix_at (derivative_matrix, irow, irow) - entry);
      }
    }
  }
  return derivative_matrix;
}

void
t8dg_functionbasis_interpolate_scalar_fn (const t8dg_functionbasis_t * functionbasis,
                                          t8dg_scalar_function_3d_fn function, void *scalar_fn_data, sc_array_t * dof_values)
{
  T8DG_ASSERT ((size_t) functionbasis->number_of_dof == dof_values->elem_count);
  int                 idof;
  if (t8dg_functionbasis_is_lagrange (functionbasis)) {
    for (idof = 0; idof < functionbasis->number_of_dof; idof++) {
      double              vertex[3] = { 0, 0, 0 };
      t8dg_functionbasis_get_Lagrange_vertex (functionbasis, idof, vertex);
      *(double *) sc_array_index_int (dof_values, idof) = function (vertex, scalar_fn_data);
    }
  }
  else {
    T8DG_ABORT ("Not yet implemented");
  }
}

t8dg_dmatrix_t     *
t8dg_functionbasis_get_child_interpolation_matrix (t8dg_functionbasis_t * functionbasis, int ichild)
{
  return functionbasis->interpolate_to_child_matrix[ichild];
}

t8dg_dmatrix_t     *
t8dg_functionbasis_get_derivative_matrix (t8dg_functionbasis_t * functionbasis)
{
  return functionbasis->derivative_matrix;
}

void
t8dg_functionbasis_ref (t8dg_functionbasis_t * functionbasis)
{
  T8DG_ASSERT (functionbasis != NULL);
  t8dg_refcount_ref (&functionbasis->rc);
}

void
t8dg_functionbasis_unref (t8dg_functionbasis_t ** pfunctionbasis)
{
  t8dg_functionbasis_t *functionbasis;

  T8DG_ASSERT (pfunctionbasis != NULL);
  functionbasis = *pfunctionbasis;
  T8DG_ASSERT (functionbasis != NULL);
  T8DG_ASSERT (t8dg_refcount_is_active (&functionbasis->rc));

  if (t8dg_refcount_unref (&functionbasis->rc)) {
    t8dg_functionbasis_reset (pfunctionbasis);
  }
}

void
t8dg_functionbasis_destroy (t8dg_functionbasis_t ** pfunctionbasis)
{
  t8dg_functionbasis_t *functionbasis;

  T8DG_ASSERT (pfunctionbasis != NULL);
  functionbasis = *pfunctionbasis;

  T8DG_ASSERT (functionbasis != NULL);
  T8DG_ASSERT (t8dg_refcount_is_last (&functionbasis->rc));

  t8dg_refcount_unref (&functionbasis->rc);
  t8dg_functionbasis_reset (pfunctionbasis);
}

void
t8dg_functionbasis_reset (t8dg_functionbasis_t ** pfunctionbasis)
{
  int                 itensor;

  T8DG_ASSERT (pfunctionbasis != NULL);
  t8dg_functionbasis_t *functionbasis = *pfunctionbasis;
  T8DG_ASSERT (functionbasis != NULL);

  if (functionbasis->num_tensor == 1) {
    if (t8dg_functionbasis_get_dim (functionbasis) == 1) {
      sc_array_destroy (functionbasis->barycentric_weights_1D);
    }
    if (t8dg_functionbasis_is_lagrange (functionbasis)) {
      t8dg_vertexset_unref (&functionbasis->vertexset);
      functionbasis->vertexset = NULL;
    }
    t8dg_dmatrix_destroy (&functionbasis->derivative_matrix);
    t8dg_dmatrix_destroy (&functionbasis->interpolate_to_child_matrix[0]);
    t8dg_dmatrix_destroy (&functionbasis->interpolate_to_child_matrix[1]);

  }
  else {
    for (itensor = 0; itensor < functionbasis->num_tensor; itensor++) {
      t8dg_functionbasis_unref (functionbasis->tensor_fb + itensor);
    }
  }
  functionbasis->number_of_dof = -1;
  T8DG_FREE (functionbasis);
  *pfunctionbasis = NULL;

}
