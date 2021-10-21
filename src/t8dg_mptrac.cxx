/*
 * t8dg_mptrac.c
 *
 *  Created on: October 15, 2021
 *      Author: Johannes Holke
 */

#include <t8_element_cxx.hxx>
#include <t8dg.h>
#include "t8dg_flux_implementation.h"
#include "t8dg_mptrac.h"


t8_mptrac_context_t* t8dg_mptrac_setup (const char *nc_filename)
{
  const char *mptrac_input = "blub DT_MET 21600 METBASE ei MET_DX 8 MET_DY 8";
  const int dimension = 3;
  const int uniform_level = 3;
  t8_mptrac_context_t *context = t8_mptrac_context_new (0, nc_filename, mptrac_input, dimension, uniform_level);
  double start_six_hours = 0;
  double physical_time;
  time2jsec (2011, 06, 05, start_six_hours, 00, 00, 00, &physical_time);

  t8_mptrac_read_nc (context, 1, physical_time);
  t8dg_debugf ("Initialized mptrac context.\n");

  return context;
}

void
t8dg_mptrac_flow_3D_fn (double x_vec[3], double flux_vec[3], double t, void *flux_data,
                        t8_locidx_t itree, t8_locidx_t ielement)
{
  double lat, lon, pressure;
  double physical_time;
  double start_six_hours = 0;
  time2jsec (2011, 06, 05, start_six_hours, 00, 00, 00, &physical_time);
  t8_mptrac_context_t *mptrac_context = (t8_mptrac_context_t *) flux_data;
  
  /* Handle boundary condition */
  const t8_forest_t forest = mptrac_context->forest;
  const t8_eclass_t tree_class = t8_forest_get_tree_class (forest, itree);
  t8_eclass_scheme_c *scheme = t8_forest_get_eclass_scheme (forest, tree_class);
  const t8_element_t *element = t8_forest_get_element_in_tree (forest, itree, ielement);
  int is_upper_boundary = 0;
  int is_lower_boundary = 0;
  /* Check for upper boundary of unit cube */
  if (scheme->t8_element_is_root_boundary (element, 5)) {
    /* We are at the upper boundary (z = 1) */
    is_upper_boundary = 1;
  }
  /* Check for upper boundary of unit cube */
  if (scheme->t8_element_is_root_boundary (element, 4)) {
    /* We are at the upper boundary (z = 1) */
    is_lower_boundary = 1;
  }
  
  t8_mptrac_coords_to_latlonpressure (mptrac_context, x_vec, &lat, &lon, &pressure);

  int                 ci[3];
  double              cw[3];

  /* If we are at the upper or lower boundary, set the flow in z direction to 0 */
  flux_vec[2] = 0;
  if (!is_lower_boundary) {
    /* Compute interpolation of u */
    intpol_met_time_3d (mptrac_context->mptrac_meteo1, mptrac_context->mptrac_meteo1->u,
                        mptrac_context->mptrac_meteo2, mptrac_context->mptrac_meteo2->u,
                        physical_time, pressure, lon, lat, flux_vec, ci, cw,
                        1);
    /* Compute interpolation of v */
    intpol_met_time_3d (mptrac_context->mptrac_meteo1, mptrac_context->mptrac_meteo1->v, mptrac_context->mptrac_meteo2, mptrac_context->mptrac_meteo2->v, physical_time, pressure, lon, lat, flux_vec + 1, ci, cw, 0);  /* 0 here since we can reuse the interpolation weights */

    if (!is_upper_boundary) {
      /* Compute interpolation of w */
      intpol_met_time_3d (mptrac_context->mptrac_meteo1, mptrac_context->mptrac_meteo1->w, mptrac_context->mptrac_meteo2, mptrac_context->mptrac_meteo2->w, physical_time, pressure, lon, lat, flux_vec + 2, ci, cw, 0);  /* 0 here since we can reuse the interpolation weights */
    }
  }
  else {
    /* We are at the lower boundary. Compute 2D interpolation at ground pressure */
    intpol_met_time_2d (mptrac_context->mptrac_meteo1, mptrac_context->mptrac_meteo1->us,
                        mptrac_context->mptrac_meteo2, mptrac_context->mptrac_meteo2->us,
                        physical_time, lon, lat, flux_vec, ci, cw,
                        1);
    intpol_met_time_2d (mptrac_context->mptrac_meteo1, mptrac_context->mptrac_meteo1->vs,
                        mptrac_context->mptrac_meteo2, mptrac_context->mptrac_meteo2->vs,
                        physical_time, lon, lat, flux_vec + 1, ci, cw,
                        1);
  }
}