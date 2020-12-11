#ifndef SRC_T8DG_OUTPUT_H_
#define SRC_T8DG_OUTPUT_H_
#include "t8dg.h"
#include "t8dg_dof.h"
#include <t8_forest.h>

T8DG_EXTERN_C_BEGIN ();

void                t8dg_output_write_vtk (t8dg_dof_values_t * dof_values, t8_forest_t forest, char prefix[BUFSIZ], int *vtk_count);

T8DG_EXTERN_C_END ();

#endif
