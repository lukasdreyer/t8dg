/** @file t8dg_square_matrix.h */
#ifndef SRC_T8DG_DOF_H_
#define SRC_T8DG_DOF_H_

#include "t8dg.h"
//#include "t8dg_global_values.h"

T8DG_EXTERN_C_BEGIN ();

typedef int         t8dg_dofidx_t;

typedef struct t8dg_dof_values t8dg_dof_values_t;
typedef sc_array_t  t8dg_element_dof_values_t;
typedef sc_array_t  t8dg_face_dof_values_t;

typedef struct t8dg_global_values t8dg_global_values_t; /*forward declared typedef because of circular dependency global_values <-> dof */

t8dg_element_dof_values_t *t8dg_dof_values_new_element_dof_values_view (t8dg_dof_values_t * dof_values, t8_locidx_t idata);

t8dg_face_dof_values_t *t8dg_dof_values_new_face_dof_values_view (t8dg_dof_values_t * dof_values, t8_locidx_t idata);

void                t8dg_element_dof_values_destroy (t8dg_element_dof_values_t ** dof_values);
void                t8dg_face_dof_values_destroy (t8dg_face_dof_values_t ** dof_values);

void                t8dg_dof_values_copy_from_index_to_index (t8dg_dof_values_t * src_dof, t8_locidx_t src_idata,
                                                              t8dg_dof_values_t * dest_dof, t8_locidx_t dest_idata);

t8dg_dof_values_t  *t8dg_dof_values_new (t8_forest_t forest, t8dg_global_values_t ** global_values_array);

t8dg_dof_values_t  *t8dg_dof_values_duplicate (t8dg_dof_values_t * src_dof_values);

t8dg_dof_values_t  *t8dg_dof_values_clone (t8dg_dof_values_t * src_dof_values);

int                 t8dg_dof_values_is_valid (t8dg_dof_values_t * dof_values);

void                t8dg_dof_values_axpy (double a, const t8dg_dof_values_t * x, t8dg_dof_values_t * y);

void                t8dg_dof_values_axpyz (double a, const t8dg_dof_values_t * x, const t8dg_dof_values_t * y, t8dg_dof_values_t * z);

void                t8dg_dof_values_debug_print (t8dg_dof_values_t * array);

void                t8dg_dof_values_print (t8dg_dof_values_t * array);

void                t8dg_dof_values_square_values (t8dg_dof_values_t * src, t8dg_dof_values_t * dest);

void                t8dg_dof_values_swap (t8dg_dof_values_t ** parray1, t8dg_dof_values_t ** parray2);

void                t8dg_dof_values_copy (const t8dg_dof_values_t * src, t8dg_dof_values_t * dest);

void                t8dg_dof_values_set_zero (t8dg_dof_values_t * array);

void                t8dg_dof_values_destroy (t8dg_dof_values_t ** dof_values);

void                t8dg_dof_values_partition (t8dg_dof_values_t * dof_values_old, t8dg_dof_values_t * dof_values_partition);

void                t8dg_dof_values_ghost_exchange (t8dg_dof_values_t * dof_values);

double             *t8dg_dof_values_get_double_pointer (const t8dg_dof_values_t * dof_values, t8_locidx_t idata);

t8dg_element_dof_values_t *t8dg_element_dof_values_duplicate (const t8dg_element_dof_values_t * element_dof_values);

t8dg_element_dof_values_t *t8dg_element_dof_values_clone (const t8dg_element_dof_values_t * element_dof_values);

t8dg_face_dof_values_t *t8dg_face_dof_values_duplicate (const t8dg_face_dof_values_t * face_dof_values);

void                t8dg_element_dof_values_set_zero (t8dg_element_dof_values_t * dof_values);

int                 t8dg_element_dof_values_is_valid (t8dg_element_dof_values_t * element_dof_values);
int                 t8dg_face_dof_values_is_valid (t8dg_face_dof_values_t * face_dof_values);

void                t8dg_element_dof_values_axpy (double a, const t8dg_element_dof_values_t * x, t8dg_element_dof_values_t * y);

T8DG_EXTERN_C_END ();

#endif // SRC_T8DG_DOF_H_
