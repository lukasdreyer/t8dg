#include "t8dg.h"
#include "t8dg_functionbasis.h"
#include "t8dg_vertexset.h"

/** The functionbasis provides the the directional derivative_matrix */
struct t8dg_functionbasis
{
  int                 dim;
  int                 number_of_dof;                    /**< Number of degrees of freedom*/
  t8dg_vertexset_t   *vertexset;                        /**< LGL vertices used as basis nodes for Lagrange nodal basis*/
  t8dg_functionbasis_type_t type;
};

t8dg_functionbasis_t *
t8dg_functionbasis_new_Lagrange (t8dg_vertexset_t * vertexset)
{
  T8DG_ASSERT (vertexset != NULL);
  T8DG_ASSERT (t8dg_vertexset_get_type (vertexset) == T8DG_LGL);
  t8dg_functionbasis_t *functionbasis = T8DG_ALLOC (t8dg_functionbasis_t, 1);
  t8dg_vertexset_ref (vertexset);
  functionbasis->vertexset = vertexset;
  functionbasis->number_of_dof = t8dg_vertexset_get_num_element_vertices (vertexset);
  functionbasis->type = T8DG_LAGRANGE_LGL;
  functionbasis->dim = t8dg_vertexset_get_dim (vertexset);
  return functionbasis;
}

void
t8dg_functionbasis_destroy (t8dg_functionbasis_t ** pfunctionbasis)
{
  t8dg_functionbasis_t *functionbasis = *pfunctionbasis;
  t8dg_vertexset_unref (functionbasis->vertexset);
  functionbasis->dim = -1;
  functionbasis->number_of_dof = -1;
  T8DG_FREE (functionbasis);
  *pfunctionbasis = NULL;
}

/*derivative_dof_values could be const aswell if not for sc_array*/
void
t8dg_functionbasis_apply_derivative_matrix_transpose (sc_array_t * dof_values,
                                                      sc_array_t * derivative_dof_values, const t8dg_functionbasis_t * functionbasis)
{
  SC_CHECK_ABORT (t8dg_functionbasis_get_dim (functionbasis) == 1 &&
                  t8dg_functionbasis_get_type (functionbasis) == T8DG_LAGRANGE_LGL, "Not yet implemented");

  double             *dof_array;
  double             *derivative_array;
  dof_array = (double *) sc_array_index (dof_values, 0);
  derivative_array = (double *) sc_array_index (derivative_dof_values, 0);

  switch (t8dg_functionbasis_get_num_dof (functionbasis)) {
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
t8dg_functionbasis_get_num_dof (const t8dg_functionbasis_t * functionbasis)
{
  T8DG_ASSERT (functionbasis != NULL);
  return functionbasis->number_of_dof;
}

void
t8dg_functionbasis_get_vertex (double vertex[3], const t8dg_functionbasis_t * functionbasis, const int idof)
{
  T8DG_ASSERT (functionbasis != NULL);
  t8dg_vertexset_get_vertex (vertex, functionbasis->vertexset, idof);
}

t8dg_functionbasis_type_t
t8dg_functionbasis_get_type (const t8dg_functionbasis_t * functionbasis)
{
  return functionbasis->type;
}

int
t8dg_functionbasis_get_dim (const t8dg_functionbasis_t * functionbasis)
{
  return functionbasis->dim;
}