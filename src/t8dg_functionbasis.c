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
