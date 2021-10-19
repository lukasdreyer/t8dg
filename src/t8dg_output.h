#ifndef SRC_T8DG_OUTPUT_H_
#define SRC_T8DG_OUTPUT_H_

#include <t8_forest.h>
#include <t8dg.h>
#include <t8dg_dof.h>


T8DG_EXTERN_C_BEGIN ();

typedef struct t8dg_vtk_data
{
  int                 vtk_count;
  const char         *prefix;
  int                 vtk_freq;
} t8dg_vtk_data_t;

t8dg_vtk_data_t    *t8dg_output_vtk_data_new (const char *prefix, int vtk_freq);

void                t8dg_output_vtk_data_destroy (t8dg_vtk_data_t ** p_vtk_data);

/** Write a vtk file with the computed solution.
 * \param [in] dof_values The dof_values to write.
 * \param [in] output_data Meta data for the output.
 * \param [in] write_flow If non-zero will also write the flow_field.
 * \param [in] flow_field Flow field to write if \a write_flow is non-zero.
 * \param [in] time       Time value to be passed to \a flow_field if \a write_flow is non-zero.
 * \param [in] flux_data  Data pointer passer to \a flow_field if \a write_flow is non-zero.
 */
void                t8dg_output_write_vtk (const t8dg_dof_values_t * dof_values, t8dg_vtk_data_t * output_data,
                                           int write_flow, t8dg_linear_flux3D_fn flow_field, const double time, void *flux_data);
T8DG_EXTERN_C_END ();

#endif