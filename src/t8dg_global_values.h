/*
 * t8dg_global_values.h
 *
 *  Created on: Apr 2, 2020
 *      Author: lukas
 */
#ifndef SRC_T8DG_GLOBAL_VALUES_H_
#define SRC_T8DG_GLOBAL_VALUES_H_

#include "t8dg_quadrature.h"
#include "t8dg_functionbasis.h"
#include "t8dg_dof.h"
#include "t8dg_quad.h"

T8DG_EXTERN_C_BEGIN ();

typedef struct t8dg_global_values t8dg_global_values_t;

/*Constructors*/

t8dg_global_values_t *t8dg_global_values_new_1D_LGL (const int number_of_LGL_vertices);

t8dg_global_values_t *t8dg_global_values_new_hypercube_LGL (const int dim, const int number_of_1D_LGL_vertices);

/* Add constructors for triangle/prism */

/*Destructor*/

void                t8dg_global_values_destroy (t8dg_global_values_t ** pvalues);

/*Vandermonde*/
void                t8dg_global_values_transform_element_dof_to_element_quad (const t8dg_global_values_t * values,
                                                                              const t8dg_element_dof_values_t * element_dof,
                                                                              t8dg_element_quad_values_t * element_quad);

/*Vandermonde transpose*/
void                t8dg_global_values_transform_element_quad_to_element_dof (const t8dg_global_values_t * values,
                                                                              const sc_array_t * element_quad_array,
                                                                              sc_array_t * element_dof_array);

/*Use from functionbasis*/
void                t8dg_global_values_element_apply_derivative_matrix_transpose
  (const t8dg_global_values_t * global_values, int idim, t8dg_element_dof_values_t * derivative_dof_values,
   t8dg_element_dof_values_t * dof_values);

void                t8dg_global_values_transform_element_dof_to_child_dof (const t8dg_global_values_t *
                                                                           global_values, t8dg_element_dof_values_t * element_dof,
                                                                           t8dg_element_dof_values_t * child_dof, const int ichild);

t8dg_functionbasis_t *t8dg_global_values_get_functionbasis (const t8dg_global_values_t * values);

t8dg_quadrature_t  *t8dg_global_values_get_quadrature (const t8dg_global_values_t * values);

/*Information getter:*/

int                 t8dg_global_values_get_num_dof (const t8dg_global_values_t * values);

int                 t8dg_global_values_get_num_faces (const t8dg_global_values_t * values);

int                 t8dg_global_values_get_num_elem_quad (const t8dg_global_values_t * values);

int                 t8dg_global_values_get_dim (const t8dg_global_values_t * values);

int                 t8dg_global_values_simplifies (const t8dg_global_values_t * global_values);

int                 t8dg_global_values_get_max_num_facevalues (const t8dg_global_values_t * global_values);

int                 t8dg_global_values_array_get_max_num_element_quad (t8dg_global_values_t * global_values_array[T8_ECLASS_COUNT]);
int                 t8dg_global_values_array_get_max_num_element_dof (t8dg_global_values_t * global_values_array[T8_ECLASS_COUNT]);

T8DG_EXTERN_C_END ();

#endif /* SRC_T8DG_GLOBAL_VALUES_H_ */