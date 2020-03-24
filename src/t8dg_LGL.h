#ifndef SRC_GLOBAL_H_
#define SRC_GLOBAL_H_

#include <sc_containers.h>
#include <sc_dmatrix.h>
#include <t8_forest.h>

#include "t8dg.h"

typedef t8dg_locidx_t t8dg_quad_idx_t;
typedef t8dg_locidx_t t8dg_dof_idx_t;

//typedef struct t8dg_LGL_vertex_set t8dg_LGL_vertex_set_t;
typedef struct t8dg_LGL_quadrature t8dg_LGL_quadrature_t;
typedef struct t8dg_LGL_functionbasis t8dg_LGL_functionbasis_t;

void                t8dg_LGL_quadrature_and_functionbasis_new_1D (t8dg_LGL_quadrature_t ** pquadrature,
                                                                  t8dg_LGL_functionbasis_t ** pfunctionbasis, int number_of_LGL);
void                t8dg_LGL_quadrature_and_functionbasis_destroy (t8dg_LGL_quadrature_t ** pquadrature,
                                                                   t8dg_LGL_functionbasis_t ** pfunctionbasis);

/* applies Lookup-table for element-dof to face-quad*/
void                t8dg_LGL_transform_element_dof_to_face_quad
  (sc_array_t * face_quad_array, const sc_array_t * element_dof_array, t8dg_LGL_quadrature_t * quadrature,
   t8dg_LGL_functionbasis_t * functionbasis);
void                t8dg_LGL_transform_face_quad_to_element_dof (sc_array_t * element_dof_array,
                                                                 const sc_array_t * face_quad_array,
                                                                 t8dg_LGL_quadrature_t * quadrature,
                                                                 t8dg_LGL_functionbasis_t * functionbasis);

void                t8dg_LGL_functionbasis_apply_derivative_matrix (sc_array_t * derivative_dof_values,
                                                                    const sc_array_t dof_values, t8dg_LGL_functionbasis_t * functionbasis);
void                t8dg_LGL_functionbasis_apply_derivative_matrix_transpose (sc_array_t * dof_values,
                                                                              const sc_array_t * derivative_dof_values,
                                                                              t8dg_LGL_functionbasis_t * functionbasis);

t8dg_quad_idx_t     t8dg_LGL_functionbasis_get_num_dof (t8dg_LGL_functionbasis_t * functionbasis);
void                t8dg_LGL_functionbasis_get_vertex (double vertex[3], t8dg_LGL_functionbasis_t * functionbasis, t8dg_dof_idx_t idof);

int                 t8dg_LGL_quadrature_get_num_faces (t8dg_LGL_quadrature_t * quadrature);
t8dg_quad_idx_t     t8dg_LGL_quadrature_get_num_element_vertices (t8dg_LGL_quadrature_t * quadrature);
t8dg_quad_idx_t     t8dg_LGL_quadrature_get_num_face_vertices (t8dg_LGL_quadrature_t * quadrature, int iface);
void                t8dg_LGL_quadrature_get_element_vertex (double vertex[3], t8dg_LGL_quadrature_t * quadrature, t8dg_quad_idx_t iquad);
double              t8dg_LGL_quadrature_get_element_weight (t8dg_LGL_quadrature_t * quadrature, t8dg_quad_idx_t iquad);

void                t8dg_LGL_quadrature_get_face_vertex (double vertex[3], t8dg_LGL_quadrature_t * quadrature,
                                                         int iface, t8dg_quad_idx_t iquad);
double              t8dg_LGL_quadrature_get_face_weight (t8dg_LGL_quadrature_t * quadrature, int iface, t8dg_quad_idx_t iquad);

#endif /* SRC_GLOBAL_H_ */
