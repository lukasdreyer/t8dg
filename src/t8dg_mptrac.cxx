/*
 * t8dg_mptrac.c
 *
 *  Created on: October 15, 2021
 *      Author: Johannes Holke
 */

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
t8dg_mptrac_flow_3D_fn (double x_vec[3], double flux_vec[3], double t, void *flux_data)
{
  double lat, lon, pressure;
  double physical_time;
  double start_six_hours = 0;
  time2jsec (2011, 06, 05, start_six_hours, 00, 00, 00, &physical_time);
  t8_mptrac_context_t *mptrac_context = (t8_mptrac_context_t *) flux_data;
  t8_mptrac_coords_to_latlonpressure (mptrac_context, x_vec, &lat, &lon, &pressure);

  int                 ci[3];
  double              cw[3];
  /* Compute interpolation of u */
  intpol_met_time_3d (mptrac_context->mptrac_meteo1, mptrac_context->mptrac_meteo1->u,
                      mptrac_context->mptrac_meteo2, mptrac_context->mptrac_meteo2->u,
                      physical_time, pressure, lon, lat, flux_vec, ci, cw,
                      1);
  /* Compute interpolation of v */
  intpol_met_time_3d (mptrac_context->mptrac_meteo1, mptrac_context->mptrac_meteo1->v, mptrac_context->mptrac_meteo2, mptrac_context->mptrac_meteo2->v, physical_time, pressure, lon, lat, flux_vec + 1, ci, cw, 0);  /* 0 here since we can reuse the interpolation weights */
  /* Compute interpolation of w */
  intpol_met_time_3d (mptrac_context->mptrac_meteo1, mptrac_context->mptrac_meteo1->w, mptrac_context->mptrac_meteo2, mptrac_context->mptrac_meteo2->w, physical_time, pressure, lon, lat, flux_vec + 2, ci, cw, 0);  /* 0 here since we can reuse the interpolation weights */
}