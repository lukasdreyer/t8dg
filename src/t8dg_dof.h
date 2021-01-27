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

t8dg_element_dof_values_t *t8dg_dof_values_new_element_dof_values_view (t8dg_dof_values_t * dof_values, t8_locidx_t itree,
                                                                        t8_locidx_t ielement);

t8dg_element_dof_values_t *t8dg_dof_values_new_element_dof_values_view_idata_eclass (t8dg_dof_values_t * dof_values, t8_locidx_t idata,
                                                                                     t8_eclass_t eclass);

t8dg_face_dof_values_t *t8dg_dof_values_new_face_dof_values_view_idata_eclass (t8dg_dof_values_t * dof_values, int iface, t8_locidx_t idata,
                                                                               t8_eclass_t element_eclass);

t8dg_face_dof_values_t *t8dg_dof_values_new_face_dof_values_view (t8dg_dof_values_t * dof_values, int iface, t8_locidx_t itree,
                                                                  t8_locidx_t ielement);

void                t8dg_element_dof_values_destroy (t8dg_element_dof_values_t ** dof_values);
void                t8dg_face_dof_values_destroy (t8dg_face_dof_values_t ** dof_values);

void                t8dg_dof_values_copy_from_index_to_index (t8dg_dof_values_t * src_dof, t8_locidx_t src_idata,
                                                              t8dg_dof_values_t * dest_dof, t8_locidx_t dest_idata);

t8dg_dof_values_t  *t8dg_dof_values_new (t8_forest_t forest, t8dg_global_values_t ** global_values_array);

t8dg_dof_values_t  *t8dg_dof_values_new_data_local (t8_forest_t forest, t8dg_global_values_t ** global_values_array, double *array,
                                                    t8dg_dofidx_t num_total_values);

t8dg_dof_values_t  *t8dg_dof_values_duplicate (t8dg_dof_values_t * src_dof_values);

t8dg_dof_values_t  *t8dg_dof_values_clone (t8dg_dof_values_t * src_dof_values);

int                 t8dg_dof_values_equal (t8dg_dof_values_t * dof_values, t8dg_dof_values_t * dof_values_compare);

int                 t8dg_element_dof_values_equal (t8dg_element_dof_values_t * element_dof_values,
                                                   t8dg_element_dof_values_t * element_dof_values_compare);

int                 t8dg_dof_values_is_valid (t8dg_dof_values_t * dof_values);

void                t8dg_dof_values_axpy (double a, const t8dg_dof_values_t * x, t8dg_dof_values_t * y);

void                t8dg_dof_values_axpyz (double a, const t8dg_dof_values_t * x, const t8dg_dof_values_t * y, t8dg_dof_values_t * z);

void                t8dg_dof_values_debug_print (t8dg_dof_values_t * array);

void                t8dg_dof_values_print (t8dg_dof_values_t * array);

void                t8dg_dof_values_square_values (t8dg_dof_values_t * src, t8dg_dof_values_t * dest);

void                t8dg_dof_values_swap (t8dg_dof_values_t ** parray1, t8dg_dof_values_t ** parray2);

void                t8dg_dof_values_copy (t8dg_dof_values_t * src, t8dg_dof_values_t * dest);

void                t8dg_dof_values_set_zero (t8dg_dof_values_t * array);

void                t8dg_dof_values_set_all_values (t8dg_dof_values_t * dof_values, double value);

void                t8dg_dof_values_destroy (t8dg_dof_values_t ** dof_values);

void                t8dg_dof_values_partition (t8dg_dof_values_t * dof_values_old, t8dg_dof_values_t * dof_values_partition);

void                t8dg_dof_values_ghost_exchange (t8dg_dof_values_t * dof_values);

/*replace by set/get value?*/
double             *t8dg_dof_values_get_double_pointer (const t8dg_dof_values_t * dof_values, t8_locidx_t idata);

t8dg_element_dof_values_t *t8dg_element_dof_values_duplicate (t8dg_element_dof_values_t * element_dof_values);

t8dg_element_dof_values_t *t8dg_element_dof_values_clone (t8dg_element_dof_values_t * element_dof_values);

void                t8dg_element_dof_values_copy (t8dg_element_dof_values_t * src, t8dg_element_dof_values_t * dest);

t8dg_face_dof_values_t *t8dg_face_dof_values_duplicate (t8dg_face_dof_values_t * face_dof_values);

double              t8dg_face_dof_values_get_value (t8dg_face_dof_values_t * face_dof_values, t8dg_dofidx_t idof);
void                t8dg_face_dof_values_set_value (t8dg_face_dof_values_t * face_dof_values, t8dg_dofidx_t idof, double value);

void                t8dg_face_dof_values_set_zero (t8dg_face_dof_values_t * face_dof_values);
void                t8dg_face_dof_values_axpy (double a, t8dg_face_dof_values_t * x, t8dg_face_dof_values_t * y);

void                t8dg_element_dof_values_set_zero (t8dg_element_dof_values_t * dof_values);

int                 t8dg_element_dof_values_is_valid (t8dg_element_dof_values_t * element_dof_values);
int                 t8dg_face_dof_values_is_valid (t8dg_face_dof_values_t * face_dof_values);

void                t8dg_element_dof_values_axpy (double a, t8dg_element_dof_values_t * x, t8dg_element_dof_values_t * y);
void                t8dg_element_dof_values_axpyz (double a, t8dg_element_dof_values_t * x, t8dg_element_dof_values_t * y,
                                                   t8dg_element_dof_values_t * z);

double              t8dg_element_dof_values_get_value (t8dg_element_dof_values_t * element_dof_values, t8dg_dofidx_t idof);

void                t8dg_element_dof_values_set_value (t8dg_element_dof_values_t * element_dof_values, t8dg_dofidx_t idof, double value);

void                t8dg_element_dof_values_set_all (t8dg_element_dof_values_t * element_dof_values, double value);
t8dg_element_dof_values_t *t8dg_element_dof_values_new (t8dg_dofidx_t num_dof);

t8dg_dofidx_t       t8dg_element_dof_values_get_num_dof (t8dg_element_dof_values_t * element_dof_values);

double              t8dg_element_dof_values_element_norm_infty (t8dg_element_dof_values_t * element_dof_values);

void                t8dg_element_dof_values_square_values (t8dg_element_dof_values_t * element_dof_values,
                                                           t8dg_element_dof_values_t * element_dof_square_values);

void                t8dg_element_dof_values_debug_print (t8dg_element_dof_values_t * element_dof_values);

t8dg_face_dof_values_t *t8dg_face_dof_values_new (t8dg_dofidx_t num_dof);

void                t8dg_dof_values_add (t8dg_dof_values_t * sum, t8dg_dof_values_t * summand);

void                t8dg_dof_values_subtract (t8dg_dof_values_t * tally, t8dg_dof_values_t * subtrahend);

void                t8dg_dof_values_ax (t8dg_dof_values_t * x, double a);

t8_forest_t         t8dg_dof_values_get_forest (t8dg_dof_values_t * dof_values);

void                t8dg_face_dof_values_orient (t8dg_face_dof_values_t * face_dof_values, t8_eclass_t eclass_face, int orientation);
void                t8dg_face_dof_values_orient_back (t8dg_face_dof_values_t * face_dof_values, t8_eclass_t eclass_face, int orientation);

double              t8dg_dof_values_get_max_value (t8dg_dof_values_t * dof_values, t8_locidx_t itree, t8_locidx_t ielement);
double              t8dg_dof_values_get_min_value (t8dg_dof_values_t * dof_values, t8_locidx_t itree, t8_locidx_t ielement);

int                 t8dg_element_dof_values_is_init (t8dg_element_dof_values_t * element_dof_values);

T8DG_EXTERN_C_END ();

#endif // SRC_T8DG_DOF_H_
