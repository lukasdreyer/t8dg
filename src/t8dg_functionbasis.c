#include "t8dg.h"
#include "t8dg_functionbasis.h"
#include "t8dg_vertexset.h"
#include <sc_containers.h>
#include "t8dg_dmatrix.h"

/** The functionbasis provides the the directional derivative_matrix and interpolation/projection */
struct t8dg_functionbasis
{
  int                 number_of_dof;                    /**< Number of degrees of freedom*/
  t8dg_functionbasis_type_t type;

  t8_eclass_t         element_class;

  /** Vertexset used for 1D Lagrange Basis*/
  t8dg_vertexset_t   *vertexset;
  /** barycentric weights used for 1D Lagrange Basis */
  sc_array_t         *barycentric_weights_1D;

  /*precalculated values for Triangle ? */

  /** Tensor Information */
  int                 num_tensor;
                        /**< if == 1, line, tri or tet (or pyramid) */
  t8dg_functionbasis_t *tensor_fb[DIM3];
  int                 tensor_num_dof[DIM3];

  t8dg_dmatrix_t     *interpolate_to_children[MAX_SUBELEMENTS];
};

static int
t8dg_functionbasis_is_lagrange (const t8dg_functionbasis_t * functionbasis)
{
  return functionbasis->type == T8DG_LAGRANGE_LGL_1D || functionbasis->type == T8DG_LAGRANGE_GL_1D;
}

#if 0
static void
t8dg_functionbasis_fill_barycentric_weights_1D (sc_array_t * barycentric_weights, t8dg_vertexset_t * vertexset)
{
  T8DG_ABORT ("Not yet implemented");
  return;
}
#endif

t8dg_functionbasis_t *
t8dg_functionbasis_new_1D_Lagrange (t8dg_vertexset_t * vertexset)
{
  T8DG_ASSERT (vertexset != NULL);
  T8DG_ASSERT (t8dg_vertexset_get_eclass (vertexset) == T8_ECLASS_LINE);

  t8dg_functionbasis_t *functionbasis;
  functionbasis = T8DG_ALLOC_ZERO (t8dg_functionbasis_t, 1);
  t8dg_vertexset_ref (vertexset);
  functionbasis->vertexset = vertexset;
  functionbasis->element_class = t8dg_vertexset_get_eclass (vertexset);
  functionbasis->num_tensor = 1;
  functionbasis->number_of_dof = t8dg_vertexset_get_num_vertices (vertexset);
  functionbasis->barycentric_weights_1D = sc_array_new_count (sizeof (double), functionbasis->number_of_dof);
//  t8dg_functionbasis_fill_barycentric_weights_1D(functionbasis->barycentric_weights_1D,functionbasis->vertexset);
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

void
t8dg_functionbasis_destroy (t8dg_functionbasis_t ** pfunctionbasis)
{
  int                 itensor;
  t8dg_functionbasis_t *functionbasis = *pfunctionbasis;
  t8dg_vertexset_unref (&functionbasis->vertexset);
  functionbasis->vertexset = NULL;
  if (functionbasis->num_tensor == 1) {
    if (t8dg_functionbasis_get_dim (functionbasis) == 1) {
      sc_array_destroy (functionbasis->barycentric_weights_1D);
    }
    else {
      for (itensor = 0; itensor < functionbasis->num_tensor; itensor++) {
//            t8dg_functionbasis_unref(functionbasis->tensor_fb + itensor);
      }
    }
  }
  functionbasis->number_of_dof = -1;
  T8DG_FREE (functionbasis);
  *pfunctionbasis = NULL;
}

int
t8dg_functionbasis_get_num_dof (const t8dg_functionbasis_t * functionbasis)
{
  T8DG_ASSERT (functionbasis != NULL);
  return functionbasis->number_of_dof;
}

/*TODO ASSERT*/
static void
t8dg_functionbasis_get_Lagrange_vertex (const t8dg_functionbasis_t * functionbasis, const int idof, double vertex[3])
{
  int                 itensor;
  int                 idoftensor[DIM3];
  int                 startdim = 0;
  T8DG_ASSERT (functionbasis != NULL);
  T8DG_ASSERT (t8dg_functionbasis_is_lagrange (functionbasis));

  vertex[0] = 0;
  vertex[1] = 0;
  vertex[2] = 0;

  if (functionbasis->num_tensor == 1) {
    t8dg_vertexset_fill_vertex3D (functionbasis->vertexset, idof, 0, vertex);
  }
  else {
    t8dg_transform_3tensoridx (idof, functionbasis->tensor_num_dof, idoftensor);
     /*TODO*/ for (itensor = 0; itensor < functionbasis->num_tensor; itensor++) {
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

double
t8dg_functionbasis_evaluate (t8dg_functionbasis_t * functionbasis, double vertex[3])
{
  T8DG_ABORT ("Not yet implemented!");
  return 0;
}

void
t8dg_functionbasis_interpolate_scalar_fn (const t8dg_functionbasis_t * functionbasis,
                                          t8dg_scalar_function_3d_fn function, void *scalar_fn_data, sc_array_t * dof_values)
{
  int                 idof;
  if (t8dg_functionbasis_is_lagrange (functionbasis)) {
    for (idof = 0; idof < functionbasis->number_of_dof; idof++) {
      double              vertex[3] = { 0, 0, 0 };
      t8dg_functionbasis_get_Lagrange_vertex (functionbasis, idof, vertex);
      *(double *) sc_array_index_int (dof_values, idof) = function (vertex, scalar_fn_data);
    }
  }
}

void
t8dg_functionbasis_interpolate_to_child (t8dg_functionbasis_t * functionbasis, int ichild, sc_array_t * dof_values,
                                         sc_array_t * child_dof_values)
{
  if (t8dg_functionbasis_get_eclass (functionbasis) == T8_ECLASS_LINE && t8dg_functionbasis_is_lagrange (functionbasis)) {

  }
  /*here only apply precalculated matrix */
  t8dg_dmatrix_mult_sc_array (functionbasis->interpolate_to_children[ichild], dof_values, child_dof_values);

}
