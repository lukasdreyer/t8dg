#ifndef SRC_T8DG_ADAPT_H_
#define SRC_T8DG_ADAPT_H_

#include "t8dg.h"
#include "t8dg_dof.h"
#include "t8dg_values.h"

T8DG_EXTERN_C_BEGIN ();

typedef struct t8dg_adapt_data
{
  int                 maximum_refinement_level;
  int                 uniform_refinement_level;
  t8dg_dof_values_t  *dof_values;
  t8dg_dof_values_t  *dof_values_adapt;
  t8dg_values_t      *dg_values;
} t8dg_adapt_data_t;

void
 
 
 
 
 
 
 
 t8dg_adapt_replace (t8_forest_t forest_old,
                     t8_forest_t forest_new,
                     t8_locidx_t itree,
                     t8_eclass_scheme_c * ts,
                     int num_elems_old, t8_locidx_t first_ielem_old, int num_elems_new, t8_locidx_t first_ielem_new);

int
 
 
 
 
 
 
 
 
t8dg_adapt_gradient (t8_forest_t forest,
                     t8_forest_t forest_from,
                     t8_locidx_t itree, t8_locidx_t lelement_id, t8_eclass_scheme_c * ts, int num_elements, t8_element_t * elements[]);

int
 
 
 
 
 
 
 
 
t8dg_adapt_mass (t8_forest_t forest,
                 t8_forest_t forest_from,
                 t8_locidx_t itree, t8_locidx_t ielement, t8_eclass_scheme_c * ts, int num_elements, t8_element_t * elements[]);

T8DG_EXTERN_C_END ();

#endif
