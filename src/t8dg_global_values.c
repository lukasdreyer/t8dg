#include "t8dg.h"
#include "t8dg_global_values.h"
#include "t8dg_quadrature.h"
#include "t8dg_functionbasis.h"
#include "t8dg_sc_array.h"

/**precomputed values that need to be precalculated for each element type*/
struct t8dg_global_values
{
  t8_eclass_t         element_class;
  t8dg_quadrature_t  *quadrature;
  t8dg_functionbasis_t *functionbasis;
  int                 vandermonde_is_identity;
};

t8dg_global_values_t *
t8dg_global_values_new_1D_LGL (const int number_of_LGL_vertices)
{
  t8dg_global_values_t *values;
  t8dg_vertexset_t   *vertices;
  t8dg_quadrature_t  *quadrature;
  t8dg_functionbasis_t *functionbasis;

  vertices = t8dg_vertexset_new_1D_LGL (number_of_LGL_vertices);
  quadrature = t8dg_quadrature_new_vertexset (vertices, 1);
  functionbasis = t8dg_functionbasis_new_1D_lagrange (vertices, 1);
  t8dg_vertexset_unref (&vertices);

  values = T8DG_ALLOC (t8dg_global_values_t, 1);
  values->quadrature = quadrature;
  values->functionbasis = functionbasis;
  values->vandermonde_is_identity = 1;
  values->element_class = T8_ECLASS_LINE;
  return values;
}

t8dg_global_values_t *
t8dg_global_values_new_hypercube_LGL (const int dim, const int number_of_1D_LGL_vertices)
{
  t8dg_global_values_t *values;
  t8dg_vertexset_t   *vertices;
  t8dg_quadrature_t  *quadrature;
  t8dg_functionbasis_t *functionbasis;

  vertices = t8dg_vertexset_new_1D_LGL (number_of_1D_LGL_vertices);
  quadrature = t8dg_quadrature_new_hypercube (dim, vertices, 1);
  functionbasis = t8dg_functionbasis_new_hypercube_lagrange (dim, vertices, 1);
  t8dg_vertexset_unref (&vertices);

  values = T8DG_ALLOC (t8dg_global_values_t, 1);
  values->quadrature = quadrature;
  values->functionbasis = functionbasis;
  values->vandermonde_is_identity = 1;
  values->element_class = t8dg_quadrature_get_eclass (values->quadrature);
  return values;
}

void
t8dg_global_values_destroy (t8dg_global_values_t ** pvalues)
{
  T8DG_ASSERT (pvalues != NULL);
  t8dg_global_values_t *values;
  values = *pvalues;
  T8DG_ASSERT (values != NULL);

  t8dg_functionbasis_destroy (&values->functionbasis);
  t8dg_quadrature_destroy (&values->quadrature);
  T8DG_FREE (*pvalues);
  *pvalues = NULL;
}

void
t8dg_global_values_element_apply_derivative_matrix_transpose (const t8dg_global_values_t * global_values, int idim,
                                                              sc_array_t * derivative_dof_values, sc_array_t * dof_values)
{
  t8dg_functionbasis_apply_derivative_matrix_transpose (global_values->functionbasis, idim, derivative_dof_values, dof_values);
}

void
t8dg_global_values_transform_element_dof_to_child_dof (const t8dg_global_values_t * global_values,
                                                       sc_array_t * element_dof, sc_array_t * child_dof, const int ichild)
{
  t8dg_functionbasis_apply_child_interpolation_matrix (t8dg_global_values_get_functionbasis (global_values), ichild,
                                                       element_dof, child_dof);
}

/*Getter*/

int
t8dg_global_values_get_max_num_facevalues (const t8dg_global_values_t * global_values)
{
  int                 iface, max_value = 0;
  t8dg_functionbasis_t *functionbasis;
  t8dg_quadrature_t  *quadrature;
  functionbasis = t8dg_global_values_get_functionbasis (global_values);
  quadrature = t8dg_global_values_get_quadrature (global_values);
  for (iface = 0; iface < t8_eclass_num_faces[global_values->element_class]; iface++) {
    max_value = SC_MAX (max_value, t8dg_functionbasis_get_num_face_dof (functionbasis, iface));
    max_value = SC_MAX (max_value, t8dg_quadrature_get_num_face_vertices (quadrature, iface));
  }
  return max_value;
}

int
t8dg_global_values_simplifies (const t8dg_global_values_t * global_values)
{
  return global_values->vandermonde_is_identity;
}

int
t8dg_global_values_get_num_dof (const t8dg_global_values_t * values)
{
  return t8dg_functionbasis_get_num_dof (values->functionbasis);
}

int
t8dg_global_values_get_num_elem_quad (const t8dg_global_values_t * values)
{
  return t8dg_quadrature_get_num_element_vertices (values->quadrature);
}

int
t8dg_global_values_get_num_faces (const t8dg_global_values_t * values)
{
  return t8_eclass_num_faces[values->element_class];
}

t8dg_functionbasis_t *
t8dg_global_values_get_functionbasis (const t8dg_global_values_t * values)
{
  return values->functionbasis;
}

t8dg_quadrature_t  *
t8dg_global_values_get_quadrature (const t8dg_global_values_t * values)
{
  return values->quadrature;
}

int
t8dg_global_values_get_num_face_quad (const t8dg_global_values_t * global_values, int iface)
{
  return t8dg_quadrature_get_num_face_vertices (global_values->quadrature, iface);
}

int
t8dg_global_values_get_num_face_dof (const t8dg_global_values_t * global_values, int iface)
{
  return t8dg_functionbasis_get_num_face_dof (global_values->functionbasis, iface);
}

int
t8dg_global_values_get_dim (const t8dg_global_values_t * values)
{
  return t8_eclass_to_dimension[values->element_class];
}

static int
t8dg_global_values_get_max_num_element_quad (const t8dg_global_values_t * global_values)
{
  return t8dg_quadrature_get_num_element_vertices (global_values->quadrature);
}

static int
t8dg_global_values_get_max_num_element_dof (const t8dg_global_values_t * global_values)
{
  return t8dg_functionbasis_get_num_dof (global_values->functionbasis);
}

int
t8dg_global_values_array_get_max_num_element_quad (t8dg_global_values_t ** global_values_array)
{
  int                 eclass;
  int                 max_num_element_quad = 0;
  for (eclass = 0; eclass < T8_ECLASS_COUNT; eclass++) {
    if (global_values_array[eclass] != NULL) {
      max_num_element_quad = SC_MAX (max_num_element_quad, t8dg_global_values_get_max_num_element_quad (global_values_array[eclass]));
    }
  }
  T8DG_ASSERT (max_num_element_quad > 0);
  return max_num_element_quad;
}

int
t8dg_global_values_array_get_max_num_element_dof (t8dg_global_values_t ** global_values_array)
{
  int                 eclass;
  int                 max_num_element_dof = 0;
  for (eclass = 0; eclass < T8_ECLASS_COUNT; eclass++) {
    if (global_values_array[eclass] != NULL) {
      max_num_element_dof = SC_MAX (max_num_element_dof, t8dg_global_values_get_max_num_element_dof (global_values_array[eclass]));
    }
  }
  T8DG_ASSERT (max_num_element_dof > 0);
  return max_num_element_dof;
}

static              t8dg_quadidx_t
t8dg_global_values_get_max_num_face_quad (t8dg_global_values_t * global_values)
{
  t8dg_quadidx_t      max_num_face_quad = 0;
  int                 iface, num_faces;
  num_faces = t8dg_quadrature_get_num_face_quadrature (global_values->quadrature);
  for (iface = 0; iface < num_faces; iface++) {
    max_num_face_quad = SC_MAX (max_num_face_quad, t8dg_quadrature_get_num_face_vertices (global_values->quadrature, iface));
  }
  return max_num_face_quad;
}

static              t8dg_dofidx_t
t8dg_global_values_get_max_num_face_dof (t8dg_global_values_t * global_values)
{
  t8dg_dofidx_t       max_num_face_dof = 0;
  int                 iface, num_faces;
  num_faces = t8dg_functionbasis_get_num_face_functionbasis (global_values->functionbasis);
  for (iface = 0; iface < num_faces; iface++) {
    max_num_face_dof = SC_MAX (max_num_face_dof, t8dg_functionbasis_get_num_face_dof (global_values->functionbasis, iface));
  }
  return max_num_face_dof;
}

int
t8dg_global_values_array_get_max_num_face_quad (t8dg_global_values_t ** global_values_array)
{
  int                 eclass;
  int                 max_num_face_quad = 0;
  for (eclass = 0; eclass < T8_ECLASS_COUNT; eclass++) {
    if (global_values_array[eclass] != NULL) {
      max_num_face_quad = SC_MAX (max_num_face_quad, t8dg_global_values_get_max_num_face_quad (global_values_array[eclass]));
    }
  }
  return max_num_face_quad;
}

int
t8dg_global_values_array_get_max_num_face_dof (t8dg_global_values_t ** global_values_array)
{
  int                 eclass;
  int                 max_num_face_dof = 0;
  for (eclass = 0; eclass < T8_ECLASS_COUNT; eclass++) {
    if (global_values_array[eclass] != NULL) {
      max_num_face_dof = SC_MAX (max_num_face_dof, t8dg_global_values_get_max_num_face_dof (global_values_array[eclass]));
    }
  }
  return max_num_face_dof;
}

int
t8dg_global_values_get_num_children (const t8dg_global_values_t * global_values)
{
  return t8dg_functionbasis_get_num_children (global_values->functionbasis);
}

/*currently only dependent on the elementtype of the coarse element*/
t8dg_global_values_t *
t8dg_global_values_array_get_global_values (t8dg_global_values_t ** global_values_array, t8_forest_t forest, t8_locidx_t itree,
                                            t8_locidx_t ielement)
{
  t8_eclass_t         eclass;
  eclass = t8dg_forest_get_eclass (forest, itree, ielement);
  if (global_values_array[eclass] == NULL) {
    T8DG_ABORT ("Fail\n");
  }
  return global_values_array[eclass];
}

/*Vandermonde*/
void
t8dg_global_values_transform_element_dof_to_element_quad (const t8dg_global_values_t * values,
                                                          t8dg_element_dof_values_t * element_dof,
                                                          t8dg_element_quad_values_t * element_quad)
{
  if (t8dg_global_values_simplifies (values)) {
    t8dg_element_dof_values_copy (element_dof, element_quad);   /* TODO: write own function */
  }
  else {
    T8DG_ABORT ("Not implemented \n ");
  }
}

/*Vandermonde transpose*/
void
t8dg_global_values_transform_element_quad_to_element_dof (const t8dg_global_values_t * values,
                                                          t8dg_element_quad_values_t * element_quad,
                                                          t8dg_element_dof_values_t * element_dof)
{
  if (t8dg_global_values_simplifies (values)) {
    t8dg_element_dof_values_copy (element_quad, element_dof);   /* TODO: write own function */
  }
  else {
    T8DG_ABORT ("Not implemented \n ");
  }
}
