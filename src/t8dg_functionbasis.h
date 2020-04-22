/*
 * t8dg_functionbasis.h
 *
 *  Created on: Apr 2, 2020
 *      Author: lukas
 */

#ifndef SRC_T8DG_FUNCTIONBASIS_H_
#define SRC_T8DG_FUNCTIONBASIS_H_

#include "t8dg_vertexset.h"
#include <sc_containers.h>
#include "t8dg_dmatrix.h"

typedef enum t8dg_functionbasis_type
{
  T8DG_LAGRANGE_LGL_1D,
  T8DG_LAGRANGE_GL_1D
} t8dg_functionbasis_type_t;

typedef struct t8dg_functionbasis t8dg_functionbasis_t;

t8dg_dmatrix_t     *t8dg_functionbasis_Lagrange_derivative_matrix (t8dg_functionbasis_t * functionbasis);

t8dg_dmatrix_t     *t8dg_functionbasis_Lagrange_interpolation_matrix (t8dg_functionbasis_t * functionbasis,
                                                                      t8dg_vertexset_t * interpolation_vertices);

void                t8dg_functionbasis_interpolate_scalar_fn (const t8dg_functionbasis_t * functionbasis,
                                                              t8dg_scalar_function_3d_fn function, void *scalar_fn_data,
                                                              sc_array_t * dof_values);

int                 t8dg_functionbasis_is_valid (const t8dg_functionbasis_t * functionbasis);

int                 t8dg_functionbasis_is_lagrange (const t8dg_functionbasis_t * functionbasis);

t8dg_functionbasis_type_t t8dg_functionbasis_get_type (const t8dg_functionbasis_t * functionbasis);

int                 t8dg_functionbasis_get_num_dof (const t8dg_functionbasis_t * functionbasis);

int                 t8dg_functionbasis_get_dim (const t8dg_functionbasis_t * functionbasis);

t8dg_dmatrix_t     *t8dg_functionbasis_get_derivative_matrix (t8dg_functionbasis_t * functionbasis);

t8dg_dmatrix_t     *t8dg_functionbasis_get_child_interpolation_matrix (t8dg_functionbasis_t * functionbasis, int ichild);

/** Creates and allocates a 1D Lagrange functionbasis from a vertexset.
 *
 *
 */
t8dg_functionbasis_t *t8dg_functionbasis_new_1D_Lagrange (t8dg_vertexset_t * vertexset);

/** Destroy a functionbasis.
 * \param [in] pfunctionbasis		On input, the functionbasis pointed to must exist
 *                              with reference count 1. The functionbasis is destroyed and this pointer
 *                              set to NULL.
 */
void                t8dg_functionbasis_destroy (t8dg_functionbasis_t ** pfunctionbasis);

/** Decrease the reference counter of a functionbasis.
 *  When the reference count reaches zero, the functionbasis is destroyed
 *  and in the tensor case the functionbases are dereferenced
 * \param [in] pfunctionbasis		On input, the functionbasis pointed to must exist
 *                              with positive reference count.  If the reference count reaches
 *                              zero, the functionbasis is destroyed and this pointer
 *                              set to NULL.
 *                              Otherwise, the pointer is not changed and
 *                              the functionbasis is not modified in other ways.
 */
void                t8dg_functionbasis_unref (t8dg_functionbasis_t ** pfunctionbasis);

/** Increase the reference counter of a functionbasis
 * \param [in] functionbasis		The functionbasis which gets referenced
 */
void                t8dg_functionbasis_ref (t8dg_functionbasis_t * functionbasis);

/** An estimate on the max value of the derivative on the reference element given dof_values
double                        t8dg_functionbasis_max_derivative(t8dg_functionbasis_t *functionbasis, sc_array_t *dof_values);
*/

#endif /* SRC_T8DG_FUNCTIONBASIS_H_ */
