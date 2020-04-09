#include "t8dg.h"
#include "t8dg_global_precomputed_values.h"
#include "t8dg_vertexset.h"
#include "t8dg_quadrature.h"
#include "t8dg_functionbasis.h"
#include "t8dg_sc_array.h"
//#include <sc_dmatrix.h>

/**precomputed values that need to be precalculated for each element type*/
struct t8dg_global_precomputed_values
{
  int                 dim;
  int                 number_of_faces;
  int                 number_of_element_quad_points;
//  int                 number_of_face_quad_points;
  int                 number_of_dof;
  t8dg_quadrature_t  *quadrature;
  t8dg_functionbasis_t *functionbasis;
#if 0
  sc_dmatrix_t        vandermonde;      /*if LGL for fb and quad, does not need to be allocated */
  sc_dmatrix_t        face_vandermonde[MAX_FACES];      /* if LGL for fb and quad use lookuptable instead */
  sc_dmatrix_t        interpolate_to_subelement[MAX_SUB_ELEMENTS];
#endif
};

t8dg_global_precomputed_values_t *
t8dg_global_precomputed_values_new_1D_LGL (const int number_of_LGL_vertices)
{
  t8dg_global_precomputed_values_t *values;
  t8dg_vertexset_t   *vertices;
  t8dg_quadrature_t  *quadrature;
  t8dg_functionbasis_t *functionbasis;
  vertices = t8dg_vertexset_new_1D_LGL (number_of_LGL_vertices);
  quadrature = t8dg_quadrature_new (vertices);
  functionbasis = t8dg_functionbasis_new_Lagrange (vertices);
  t8dg_vertexset_unref (&vertices);
  values = T8DG_ALLOC (t8dg_global_precomputed_values_t, 1);
  values->quadrature = quadrature;
  values->functionbasis = functionbasis;
  values->number_of_faces = t8dg_quadrature_get_num_faces (values->quadrature);
  values->number_of_element_quad_points = t8dg_quadrature_get_num_element_vertices (values->quadrature);
//  values->number_of_face_quad_points = t8dg_quadrature_get_num_face_vertices(values->quadrature);
  values->number_of_dof = t8dg_functionbasis_get_num_dof (values->functionbasis);
  values->dim = 1;
  return values;
}

void
t8dg_global_precomputed_values_destroy (t8dg_global_precomputed_values_t ** pvalues)
{
  T8DG_ASSERT (pvalues != NULL);
  t8dg_global_precomputed_values_t *values;
  values = *pvalues;
  T8DG_ASSERT (values != NULL);

  t8dg_functionbasis_destroy (&values->functionbasis);
  t8dg_quadrature_destroy (&values->quadrature);
  T8DG_FREE (*pvalues);
  *pvalues = NULL;
}

void
t8dg_global_precomputed_values_transform_element_dof_to_element_quad (const t8dg_global_precomputed_values_t * values,
                                                                      sc_array_t * element_dof_array, sc_array_t * element_quad_array)
{
  T8DG_CHECK_ABORT (t8dg_quadrature_get_type (values->quadrature) == T8DG_LGL &&
                    t8dg_functionbasis_get_type (values->functionbasis) == T8DG_LAGRANGE_LGL, "Not yet implemented");
  if (t8dg_quadrature_get_type (values->quadrature) == T8DG_LGL && t8dg_functionbasis_get_type (values->functionbasis) == T8DG_LAGRANGE_LGL) {
    t8dg_sc_array_copy (element_dof_array, element_quad_array);
    return;
  }
}

void
t8dg_global_precomputed_values_transform_element_quad_to_element_dof (const t8dg_global_precomputed_values_t * values,
                                                                      sc_array_t * element_quad_array, sc_array_t * element_dof_array)
{
  T8DG_CHECK_ABORT (t8dg_quadrature_get_type (values->quadrature) == T8DG_LGL &&
                    t8dg_functionbasis_get_type (values->functionbasis) == T8DG_LAGRANGE_LGL, "Not yet implemented");
  if (t8dg_quadrature_get_type (values->quadrature) == T8DG_LGL && t8dg_functionbasis_get_type (values->functionbasis) == T8DG_LAGRANGE_LGL) {
    t8dg_sc_array_copy (element_quad_array, element_dof_array);
    return;
  }
}

void
t8dg_global_precomputed_values_transform_element_dof_to_face_quad (const t8dg_global_precomputed_values_t * values,
                                                                   const int iface,
                                                                   const sc_array_t * element_dof_array, sc_array_t * face_quad_array)
{
  T8DG_CHECK_ABORT (t8dg_quadrature_get_type (values->quadrature) == T8DG_LGL &&
                    t8dg_functionbasis_get_type (values->functionbasis) == T8DG_LAGRANGE_LGL &&
                    t8dg_quadrature_get_num_element_vertices (values->quadrature) == t8dg_functionbasis_get_num_dof (values->functionbasis),
                    "Not yet implemented");

  T8DG_ASSERT (element_dof_array->elem_size == sizeof (double));
  T8DG_ASSERT (face_quad_array->elem_size == sizeof (double));
  T8DG_ASSERT (element_dof_array->elem_count == (size_t) t8dg_functionbasis_get_num_dof (values->functionbasis));
  T8DG_ASSERT (face_quad_array->elem_count == (size_t) t8dg_quadrature_get_num_face_vertices (values->quadrature, iface));

  t8dg_quad_idx_t     iquad, ifacequad, num_face_vertices;
  t8dg_vertexset_t   *vertexset = t8dg_quadrature_get_vertexset (values->quadrature);

  num_face_vertices = t8dg_quadrature_get_num_face_vertices (values->quadrature, iface);

  for (ifacequad = 0; ifacequad < num_face_vertices; ifacequad++) {
    iquad = t8dg_vertexset_get_LGL_facevertex_element_index (vertexset, iface, ifacequad);
    *(double *) t8dg_sc_array_index_quadidx (face_quad_array, ifacequad) =
      *(double *) t8dg_sc_array_index_quadidx (element_dof_array, iquad);
  }
}

void
t8dg_global_precomputed_values_transform_face_quad_to_element_dof (t8dg_global_precomputed_values_t * values,
                                                                   const int iface,
                                                                   const sc_array_t * face_quad_array, sc_array_t * element_dof_array)
{
  T8DG_CHECK_ABORT (t8dg_quadrature_get_type (values->quadrature) == T8DG_LGL &&
                    t8dg_functionbasis_get_type (values->functionbasis) == T8DG_LAGRANGE_LGL &&
                    t8dg_quadrature_get_num_element_vertices (values->quadrature) == t8dg_functionbasis_get_num_dof (values->functionbasis),
                    "Not yet implemented");

  T8DG_ASSERT (element_dof_array->elem_size == sizeof (double));
  T8DG_ASSERT (face_quad_array->elem_size == sizeof (double));
  T8DG_ASSERT (element_dof_array->elem_count == (size_t) t8dg_functionbasis_get_num_dof (values->functionbasis));
  T8DG_ASSERT (face_quad_array->elem_count == (size_t) t8dg_quadrature_get_num_face_vertices (values->quadrature, iface));

  t8dg_quad_idx_t     iquad, ifacequad, num_element_vertices, num_face_vertices;
  t8dg_vertexset_t   *vertexset = t8dg_quadrature_get_vertexset (values->quadrature);

  num_element_vertices = t8dg_quadrature_get_num_element_vertices (values->quadrature);
  num_face_vertices = t8dg_quadrature_get_num_face_vertices (values->quadrature, iface);

  for (iquad = 0; iquad < num_element_vertices; iquad++) {
    *(double *) t8dg_sc_array_index_quadidx (element_dof_array, iquad) = 0;
  }
  for (ifacequad = 0; ifacequad < num_face_vertices; ifacequad++) {
    iquad = t8dg_vertexset_get_LGL_facevertex_element_index (vertexset, iface, ifacequad);
    *(double *) t8dg_sc_array_index_quadidx (element_dof_array, iquad) =
      *(double *) t8dg_sc_array_index_quadidx (face_quad_array, ifacequad);
  }
}

/*derivative_dof_values could be const aswell if not for sc_array*/
void
t8dg_global_precomputed_values_element_apply_derivative_matrix_transpose (const t8dg_global_precomputed_values_t * global_values,
                                                                          sc_array_t * derivative_dof_values, sc_array_t * dof_values)
{
  SC_CHECK_ABORT (t8dg_global_precomputed_values_get_dim (global_values) == 1 &&
                  t8dg_functionbasis_get_type (global_values->functionbasis) == T8DG_LAGRANGE_LGL, "Not yet implemented");

  double             *dof_array;
  double             *derivative_array;
  dof_array = (double *) sc_array_index (dof_values, 0);
  derivative_array = (double *) sc_array_index (derivative_dof_values, 0);

  switch (global_values->number_of_dof) {
  case (1):
    dof_array[0] = 0;
    break;
  case (2):
    dof_array[0] = -1 * derivative_array[0] - 1 * derivative_array[1];
    dof_array[1] = 1 * derivative_array[0] + 1 * derivative_array[1];
    break;
  case (3):
    dof_array[0] = -3 * derivative_array[0] - 1 * derivative_array[1] + 1 * derivative_array[2];
    dof_array[1] = 4 * derivative_array[0] + 0 * derivative_array[1] - 4 * derivative_array[2];
    dof_array[2] = -1 * derivative_array[0] + 1 * derivative_array[1] + 3 * derivative_array[2];
    break;
  default:
    SC_ABORT ("derivative_matrix not yet implemented");
  }
}

int
t8dg_global_precomputed_values_get_num_dof (const t8dg_global_precomputed_values_t * values)
{
  return values->number_of_dof;
}

t8dg_quad_idx_t
t8dg_global_precomputed_values_get_num_elem_quad (const t8dg_global_precomputed_values_t * values)
{
  return values->number_of_element_quad_points;
}

#if 0
t8dg_quad_idx_t
t8dg_global_precomputed_values_get_num_face_quad (t8dg_global_precomputed_values_t * values)
{
  return values->number_of_face_quad_points;
}
#endif
int
t8dg_global_precomputed_values_get_num_faces (const t8dg_global_precomputed_values_t * values)
{
  return values->number_of_faces;
}

t8dg_functionbasis_t *
t8dg_global_precomputed_values_get_functionbasis (const t8dg_global_precomputed_values_t * values)
{
  return values->functionbasis;
}

t8dg_quadrature_t  *
t8dg_global_precomputed_values_get_quadrature (const t8dg_global_precomputed_values_t * values)
{
  return values->quadrature;
}

int
t8dg_global_precomputed_values_get_dim (const t8dg_global_precomputed_values_t * values)
{
  return values->dim;
}
