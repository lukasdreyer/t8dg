/** @file t8dg_quad.h */
#ifndef SRC_T8DG_QUAD_H_
#define SRC_T8DG_QUAD_H_

#include "t8dg.h"
//#include "t8dg_global_values.h"

T8DG_EXTERN_C_BEGIN ();

typedef int         t8dg_quadidx_t;

typedef struct t8dg_quad_values t8dg_quad_values_t;
typedef sc_array_t  t8dg_element_quad_values_t;
typedef sc_array_t  t8dg_face_quad_values_t;

typedef struct t8dg_global_values t8dg_global_values_t; /*forward declared typedef because of circular dependency global_values <-> quad */

t8dg_element_quad_values_t *t8dg_quad_values_new_element_quad_values_view (t8dg_quad_values_t * quad_values, t8_locidx_t itree,
                                                                           t8_locidx_t ielement);

t8dg_face_quad_values_t *t8dg_quad_values_new_face_quad_values_view (t8dg_quad_values_t * quad_values, t8_locidx_t itree,
                                                                     t8_locidx_t ielement);

void                t8dg_element_quad_values_destroy (t8dg_element_quad_values_t ** quad_values);
void                t8dg_face_quad_values_destroy (t8dg_face_quad_values_t ** quad_values);

void                t8dg_quad_values_copy_from_index_to_index (t8dg_quad_values_t * src_quad, t8_locidx_t src_idata,
                                                               t8dg_quad_values_t * dest_quad, t8_locidx_t dest_idata);

double              t8dg_element_quad_values_get_value (t8dg_element_quad_values_t * element_quad_values, t8dg_quadidx_t iquad);

void
               t8dg_element_quad_values_set_value (t8dg_element_quad_values_t * element_quad_values, t8dg_quadidx_t iquad, double value);

double              t8dg_face_quad_values_get_value (t8dg_face_quad_values_t * face_quad_values, t8dg_quadidx_t iquad);

void                t8dg_face_quad_values_set_value (t8dg_face_quad_values_t * face_quad_values, t8dg_quadidx_t iquad, double value);

t8dg_quadidx_t      t8dg_element_quad_values_get_num_element_quad_points (t8dg_element_quad_values_t * element_quad_values);

t8dg_quadidx_t      t8dg_face_quad_values_get_num_face_quad_points (t8dg_face_quad_values_t * face_quad_values);

t8dg_quad_values_t *t8dg_quad_values_new (t8_forest_t forest, t8dg_global_values_t * global_values_array[T8_ECLASS_COUNT]);

t8dg_quad_values_t *t8dg_quad_values_new_local (t8_forest_t forest, t8dg_global_values_t * global_values_array[T8_ECLASS_COUNT]);

void                t8dg_quad_values_destroy (t8dg_quad_values_t ** quad_values);

void                t8dg_quad_values_partition (t8dg_quad_values_t * quad_values_old, t8dg_quad_values_t * quad_values_partition);

T8DG_EXTERN_C_END ();

#endif // SRC_T8DG_QUAD_H_
