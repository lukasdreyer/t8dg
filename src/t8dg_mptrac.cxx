/*
 * t8dg_mptrac.c
 *
 *  Created on: October 15, 2021
 *      Author: Johannes Holke
 */

#include <t8_element_cxx.hxx>
#include <t8dg.h>
#include <t8dg_advect_diff_problem.h>
#include <t8dg_flux_implementation.h>
#include <t8dg_timestepping.h>
#include <t8dg_mptrac.h>


t8dg_mptrac_flux_data::t8dg_mptrac_flux_data (const char *nc_filename, int hours_between_file_reads)
: hours_between_file_reads(hours_between_file_reads)
{
  const char *mptrac_input = "blub DT_MET 21600 METBASE ei MET_DX 8 MET_DY 8";
  const int dimension = 3;
  const int uniform_level = 3;
  context = t8_mptrac_context_new (0, nc_filename, mptrac_input, dimension, uniform_level);

  start_six_hours = 0;
  time2jsec (2011, 06, 05, start_six_hours, 00, 00, 00, &physical_time_s);

  t8_mptrac_read_nc (context, 1, physical_time_s);
  hours_since_last_file_read = 0;
  t8dg_debugf ("Initialized mptrac context.\n");
}

t8dg_mptrac_flux_data::~t8dg_mptrac_flux_data ()
{
  t8_mptrac_context_destroy (&context);
}

void t8dg_mptrac_flux_data::initialize(const struct t8dg_linear_advection_diffusion_problem *problem)
{
  context->forest = t8dg_advect_diff_problem_get_forest (problem);
  t8_forest_ref (context->forest);
}


void t8dg_mptrac_flux_data::start_new_time_step (const struct t8dg_linear_advection_diffusion_problem *problem)
{
  const t8dg_timestepping_data_t *time_data = t8dg_advect_diff_problem_get_time_data (problem);
  const double time = t8dg_timestepping_data_get_current_time (time_data);
  const double delta_t_hours = t8dg_timestepping_data_get_time_step (time_data);

  t8dg_debugf ("Mptrac data entering new timestep: %i\t%f\n", t8dg_advect_diff_problem_get_stepnumber (problem), time);

  /* Add current time to internal physical time.
   * and diff since last read. */
  physical_time_s += delta_t_hours * 3600;
  hours_since_last_file_read += delta_t_hours;
  t8dg_debugf ("Mptrac. time = %f  hours_since = %f, diff = %f\n", physical_time_s, hours_since_last_file_read,
   delta_t_hours);
  if (hours_since_last_file_read > hours_between_file_reads) {
    /* Read the nc file. This will only update if physical is advanced to the
    * next file. */
    t8_mptrac_read_nc (context, 1, physical_time_s);
    /* Rescale physical time to interval [0, hours_between_file_reads] */
    hours_since_last_file_read -= (int)(hours_since_last_file_read/hours_between_file_reads) * hours_between_file_reads;
  }
}

inline const double t8dg_mptrac_flux_data::get_time () const
{
  return physical_time_s;
}

inline const t8_mptrac_context_t *t8dg_mptrac_flux_data::get_context() const
{
  return context;
}

void
t8dg_mptrac_flow_3D_fn (double x_vec[3], double flux_vec[3], double t, const t8dg_flux_data_base *flux_data,
                        t8_locidx_t itree, t8_locidx_t ielement)
{
  double lat, lon, pressure;
  T8DG_ASSERT (dynamic_cast<const t8dg_mptrac_flux_data*>(flux_data) != NULL);
  const t8dg_mptrac_flux_data *mptrac_flux_data = static_cast<const t8dg_mptrac_flux_data *> (flux_data);
  const t8_mptrac_context_t *mptrac_context = mptrac_flux_data->get_context();
  const double physical_time_s = mptrac_flux_data->get_time();
  /* Handle boundary condition */
  const t8_forest_t forest = mptrac_context->forest;
  const t8_eclass_t tree_class = t8_forest_get_tree_class (forest, itree);
  t8_eclass_scheme_c *scheme = t8_forest_get_eclass_scheme (forest, tree_class);
  const t8_element_t *element = t8_forest_get_element_in_tree (forest, itree, ielement);
  int is_upper_boundary = 0;
  int is_lower_boundary = 0;
  int is_south_boundary = 0;
  int is_north_boundary = 0;
  /* Check for upper boundary of unit cube */
  if (scheme->t8_element_is_root_boundary (element, 5)) {
    is_upper_boundary = 1;
  }
  /* Check for upper boundary of unit cube */
  if (scheme->t8_element_is_root_boundary (element, 4)) {
    is_lower_boundary = 1;
  }
  /* Check for south boundary of unit cube */
  if (scheme->t8_element_is_root_boundary (element, 2)) {
    is_south_boundary = 1;
  }
  /* Check for north boundary of unit cube */
  if (scheme->t8_element_is_root_boundary (element, 3)) {
    is_north_boundary = 1;
  }
  const int is_at_pole = is_south_boundary || is_north_boundary;
  const int is_at_top_or_bottom = is_upper_boundary || is_lower_boundary;

  /* TODO: Try using coord = 0/1 as indicator...flow in volumen will have
   *        w/v value, flow at boundary not. */
  t8_mptrac_coords_to_latlonpressure (mptrac_context, x_vec, &lat, &lon, &pressure);

  int                 ci[3];
  double              cw[3];

  /* Compute interpolation of u (zonal wind east/west direction, x axis) */
  intpol_met_time_3d (mptrac_context->mptrac_meteo1, mptrac_context->mptrac_meteo1->u,
                      mptrac_context->mptrac_meteo2, mptrac_context->mptrac_meteo2->u,
                      physical_time, pressure, lon, lat, flux_vec, ci, cw,
                      1);
  if (!is_at_pole) {
    /* Compute interpolation of v (meridional wind north/south direction, y axis) */
    intpol_met_time_3d (mptrac_context->mptrac_meteo1, mptrac_context->mptrac_meteo1->v, mptrac_context->mptrac_meteo2, mptrac_context->mptrac_meteo2->v, physical_time, pressure, lon, lat, flux_vec + 1, ci, cw, 0);  /* 0 here since we can reuse the interpolation weights */
  } 
  else {
    flux_vec[1] = 0;
  }
  if (!is_at_top_or_bottom) {
    /* Compute interpolation of w (vertical velocity, z axis) */
    intpol_met_time_3d (mptrac_context->mptrac_meteo1, mptrac_context->mptrac_meteo1->w, mptrac_context->mptrac_meteo2, mptrac_context->mptrac_meteo2->w, physical_time, pressure, lon, lat, flux_vec + 2, ci, cw, 0);  /* 0 here since we can reuse the interpolation weights */
  }
  else {
    flux_vec[2] = 0;
  }
#if 0
/* 2D Interpolation for ground pressure. Currently deactivated, may be useful later. */
  else {
    /* We are at the lower boundary. Compute 2D interpolation at ground pressure */
    intpol_met_time_2d (mptrac_context->mptrac_meteo1, mptrac_context->mptrac_meteo1->us,
                        mptrac_context->mptrac_meteo2, mptrac_context->mptrac_meteo2->us,
                        physical_time_s, lon, lat, flux_vec, ci, cw,
                        1);
    intpol_met_time_2d (mptrac_context->mptrac_meteo1, mptrac_context->mptrac_meteo1->vs,
                        mptrac_context->mptrac_meteo2, mptrac_context->mptrac_meteo2->vs,
                        physical_time_s, lon, lat, flux_vec + 1, ci, cw,
                        1);
    /* TODO: Do we also need flux_vec[2] or is it 0? */
  }
#endif
}

double
t8dg_mptrac_box_source (const double x[3], const double t, void *fn_data)
{
  const double value_inside_box = 10;
  /* Lon/Lat/pressure coordinates of the lower left corner of the box. */
  const double box_lower_left[3] = {7, 50, 2000};
  /* Dimensions of the box in degree x degree x meter */
  const double box_extend[3] = {200, 200, 800};

  double lat, lon, pressure;

  const t8dg_mptrac_flux_data *mptrac_flux_data = static_cast<const t8dg_mptrac_flux_data *> (fn_data);

  t8_mptrac_coords_to_lonlatpressure (mptrac_flux_data->get_context(), x, &lat, &lon, &pressure);

  if (box_lower_left[0] <= lon && lon <= box_lower_left[0] + box_extend[0]
   && box_lower_left[1] <= lat && lat <= box_lower_left[1] + box_extend[1]
   && box_lower_left[2] <= pressure && pressure <= box_lower_left[2] + box_extend[2])
   {
     /* The point x is inside the box. */
     return value_inside_box;
   }
   else return 0;
}