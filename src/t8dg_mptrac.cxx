/*
 * t8dg_mptrac.c
 *
 *  Created on: October 15, 2021
 *      Author: Johannes Holke
 */

#include <t8_element_cxx.hxx>
#include <t8_vec.h>
#include <t8dg.h>
#include <t8dg_advect_diff_problem.h>
#include <t8dg_flux_implementation.h>
#include <t8dg_timestepping.h>
#include <t8dg_common.h>
#include <t8dg_mptrac.h>


t8dg_mptrac_box::t8dg_mptrac_box (double center_lon, double center_lat, double center_p, double extend_lon_deg, double extend_lat_deg, double extend_p_km)
 : box_center{center_lon, center_lat, center_p}, box_extend{extend_lon_deg, extend_lat_deg, extend_p_km}
 {
  pressure_center_in_km = Z(box_center[2]);
  box_minimum[0] = box_center[0] - box_extend[0]/2;
  box_minimum[1] = box_center[1] - box_extend[1]/2;
  box_minimum[2] = pressure_center_in_km - box_extend[2]/2;
  box_maximum[0] = box_center[0] + box_extend[0]/2;
  box_maximum[1] = box_center[1] + box_extend[1]/2;
  box_maximum[2] = pressure_center_in_km + box_extend[2]/2;
 }

int t8dg_mptrac_box::lon_lat_pressure_is_in (const double lon, const double lat, const double pressure) const
{
  const double pressure_in_km = Z(pressure);

  if (box_minimum[0] <= lon && lon <= box_maximum[0]
   && box_minimum[1] <= lat && lat <= box_maximum[1]
   && box_minimum[2] <= pressure_in_km && pressure_in_km <= box_maximum[2])
   {
     /* The point is inside the box. */
     return 1;
   }
   else return 0;
}


double t8dg_mptrac_box::dist_from_center (const double lon, const double lat, const double pressure) const
{
  const double pressure_in_km = Z(pressure);
  double coords_in_km[3], coords_center_in_km[3];


  /* Convert coordinates from long lat to km */
  geo2cart(pressure_in_km, lon, lat, coords_in_km);
  geo2cart(pressure_center_in_km, box_center[0], box_center[1], coords_center_in_km);


  /* Compute the distance and return */
  return t8_vec_dist (coords_in_km, coords_center_in_km);
}

double t8dg_mptrac_box::dist_from_center_maxnorm (const double lon, const double lat, const double pressure) const
{
  const double pressure_in_km = Z(pressure);
  double coords_in_km[3], coords_center_in_km[3];
  double diff[3];

  /* Convert coordinates from long lat to km */
  geo2cart(pressure_in_km, lon, lat, coords_in_km);
  geo2cart(pressure_center_in_km, box_center[0], box_center[1], coords_center_in_km);

  /* Compute | x - y | */
  diff[0] = fabs(coords_in_km[0] - coords_center_in_km[0]);
  diff[1] = fabs(coords_in_km[1] - coords_center_in_km[1]);
  diff[2] = fabs(coords_in_km[2] - coords_center_in_km[2]);

  /* Return the maximum of diff */
  return SC_MAX (diff[0], SC_MAX (diff[1], diff[2]));
}

double t8dg_mptrac_box::horiz_dist_from_center_maxnorm (const double lon, const double lat, const double pressure) const
{
  const double pressure_in_km = Z(pressure);
  double coords_in_km[3], coords_center_in_km[3];
  double diff[2];

  /* Convert coordinates from long lat to km */
  geo2cart(pressure_in_km, lon, lat, coords_in_km);
  geo2cart(pressure_center_in_km, box_center[0], box_center[1], coords_center_in_km);

  /* Compute the first 2 coordinates of | x - y | */
  diff[0] = fabs(coords_in_km[0] - coords_center_in_km[0]);
  diff[1] = fabs(coords_in_km[1] - coords_center_in_km[1]);

  t8dg_debugf ("distance of (%f, %f) from (%f, %f): (%f, %f)\n", lon, lat, box_center[0], box_center[1], diff[0], diff[1]);

  /* Return the maximum of diff */
  return SC_MAX (diff[0], diff[1]);
}

double t8dg_mptrac_box::horiz_dist_from_center_maxnorm_deg (const double lon, const double lat) const
{
  double diff[2];

  /* Compute the first 2 coordinates of | x - y | */
  diff[0] = fabs(box_center[0] - lon);
  diff[1] = fabs(box_center[1] - lat);

  /* We need to calculate the diff for longitude modulo 360 */
  if (diff[0] > 360) {
    diff[0] -= 360; /* Normalize to [0,360) */
  }
  T8DG_ASSERT (0 <= diff[0] && diff[0] <= 360);
  /* If > 180 we need to diff with 360 instead of 0
   * (359 - 0 = 1 and not 359!)
   */
  diff[0] = diff[0] > 180 ? 360 - diff[0] : diff[0];
  T8DG_ASSERT (diff[0] >= 0);

  t8dg_debugf ("distance of (%f, %f) from (%f, %f): (%f, %f)\n", lon, lat, box_center[0], box_center[1], diff[0], diff[1]);

  /* Return the maximum of diff */
  return SC_MAX (diff[0], diff[1]);
}

double t8dg_mptrac_box::vert_dist_from_center_maxnorm (const double pressure) const
{
  const double pressure_in_km = Z(pressure);

  /* Return the maximum of diff */
  return fabs (pressure_in_km - pressure_center_in_km);
}


t8dg_mptrac_flux_data::t8dg_mptrac_flux_data (const char *nc_filename, int hours_between_file_reads, sc_MPI_Comm comm)
: comm (comm), hours_between_file_reads(hours_between_file_reads), point_source(0, 0, 500, 10, 10, 1)
{
  const char *mptrac_input = "blub DT_MET 21600 METBASE ei MET_DX 1 MET_DY 1";
  const int dimension = 3;
  const int uniform_level = 3;
  context = t8_mptrac_context_new (0, nc_filename, mptrac_input, dimension, uniform_level, comm);

  start_six_hours = 0;
  time2jsec (2017, 01, 01, start_six_hours, 00, 00, 00, &physical_time_s);

  t8_mptrac_read_nc (context, 1, physical_time_s, comm);
  hours_since_last_file_read = 0;
  t8dg_global_productionf ("Initialized mptrac context.\n");
}

t8dg_mptrac_flux_data::~t8dg_mptrac_flux_data ()
{
  t8_mptrac_context_destroy (&context, comm);
}

void t8dg_mptrac_flux_data::initialize(const struct t8dg_linear_advection_diffusion_problem *problem)
{
  context->forest = t8dg_advect_diff_problem_get_forest (problem);
  t8_forest_ref (context->forest);
}

void t8dg_mptrac_flux_data::before_first_call_on_element (t8_forest_t forest, t8_locidx_t itree, t8_locidx_t element_in_tree)
{  
  /* Determine whether this element is at a pole or at the top/bottom. */
  int is_upper_boundary = 0;
  int is_lower_boundary = 0;
  int is_south_boundary = 0;
  int is_north_boundary = 0;
  const t8_eclass_t tree_class = t8_forest_get_tree_class (forest, itree);
  t8_eclass_scheme_c *scheme = t8_forest_get_eclass_scheme (forest, tree_class);
  const t8_element_t *element = t8_forest_get_element_in_tree (forest, itree, element_in_tree);

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
  current_element_is_at_pole = is_south_boundary || is_north_boundary;
  current_element_is_at_top_or_bottom = is_upper_boundary || is_lower_boundary;
}

bool t8dg_mptrac_flux_data::is_current_element_at_pole () const
{
  return current_element_is_at_pole;
}
    
bool t8dg_mptrac_flux_data::is_current_element_at_top_or_bottom () const
{
  return current_element_is_at_top_or_bottom;
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
    t8_mptrac_read_nc (context, 1, physical_time_s, comm);
    /* Rescale physical time to interval [0, hours_between_file_reads] */
    hours_since_last_file_read -= (int)(hours_since_last_file_read/hours_between_file_reads) * hours_between_file_reads;
  }
}

const double t8dg_mptrac_flux_data::get_time () const
{
  return physical_time_s;
}

const t8_mptrac_context_t *t8dg_mptrac_flux_data::get_context() const
{
  return context;
}

int t8dg_mptrac_flux_data::point_is_in_box (const double x[3]) const
{
  double lon, lat, pressure;

  t8_mptrac_coords_to_lonlatpressure (context, x, &lon, &lat, &pressure);

  return point_source.lon_lat_pressure_is_in (lon, lat, pressure);
}

double t8dg_mptrac_flux_data::point_distance_from_box_center (const double x[3]) const
{
  double lon, lat, pressure;

  t8_mptrac_coords_to_lonlatpressure (context, x, &lon, &lat, &pressure);

  return point_source.dist_from_center (lon, lat, pressure);
}

double t8dg_mptrac_flux_data::point_distance_from_box_center_maxnorm (const double x[3]) const
{
  double lon, lat, pressure;

  t8_mptrac_coords_to_lonlatpressure (context, x, &lon, &lat, &pressure);

  return point_source.dist_from_center_maxnorm (lon, lat, pressure);
}

void t8dg_mptrac_flux_data::point_max_vert_and_horiz_distance_from_box (const double x[3], double *vert_dist, double *horiz_dist) const
{
  double lon, lat, pressure;

  t8_mptrac_coords_to_lonlatpressure (context, x, &lon, &lat, &pressure);

  *vert_dist = point_source.vert_dist_from_center_maxnorm (pressure);
  *horiz_dist = point_source.horiz_dist_from_center_maxnorm_deg (lon, lat);
}

void
t8dg_mptrac_flow_3D_fn (double x_vec[3], double flux_vec[3], double t, const t8dg_flux_data_base *flux_data)
{
  double lat, lon, pressure;
  T8DG_ASSERT (dynamic_cast<const t8dg_mptrac_flux_data*>(flux_data) != NULL);
  const t8dg_mptrac_flux_data *mptrac_flux_data = static_cast<const t8dg_mptrac_flux_data *> (flux_data);
  const t8_mptrac_context_t *mptrac_context = mptrac_flux_data->get_context();
  const double physical_time_s = mptrac_flux_data->get_time();

  /* TODO: Try using coord = 0/1 as indicator...flow in volumen will have
   *        w/v value, flow at boundary not. */
  t8_mptrac_coords_to_lonlatpressure (mptrac_context, x_vec, &lon, &lat, &pressure);

  int                 ci[3];
  double              cw[3];
  const met_t        *meteo1 =
    (const met_t *) t8_shmem_array_index (mptrac_context->mptrac_meteo, 0);
  const met_t        *meteo2 =
    (const met_t *) t8_shmem_array_index (mptrac_context->mptrac_meteo, 1);
  /* intpol_met_time_3d does not modify meteo1 and meteo2, but also does
   * not declare them const. Hence we need to cast away the constness manually. */
  met_t *meteo1_noconst = (met_t *)     meteo1;
  met_t *meteo2_noconst = (met_t *)     meteo2;


  /* Compute interpolation of u (zonal wind east/west direction, x axis) */
  intpol_met_time_3d (meteo1_noconst, meteo1_noconst->u,
                      meteo2_noconst, meteo2_noconst->u,
                      physical_time_s, pressure, lon, lat, flux_vec, ci, cw,
                      1);
  if (!mptrac_flux_data->is_current_element_at_pole()) {
    /* Compute interpolation of v (meridional wind north/south direction, y axis) */
    intpol_met_time_3d (meteo1_noconst, meteo1_noconst->v, 
                        meteo2_noconst, meteo2_noconst->v, 
                        physical_time_s, pressure, lon, lat, flux_vec + 1, ci, cw, 0);  /* 0 here since we can reuse the interpolation weights */
  } 
  else {
    flux_vec[1] = 0;
  }
  if (!mptrac_flux_data->is_current_element_at_top_or_bottom()) {
    /* Compute interpolation of w (vertical velocity, z axis) */
    intpol_met_time_3d (meteo1_noconst, meteo1_noconst->w,
                        meteo2_noconst, meteo2_noconst->w, 
                        physical_time_s, pressure, lon, lat, flux_vec + 2, ci, cw, 0);  /* 0 here since we can reuse the interpolation weights */
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
  /* Currently deactivated to properly test indicator function. */
  return 0;
  /* We only spit out material for the first minute. */
  if (t > 1./60) {
    return 0;
  }
  /* The flux value inside the box */
  const double value_inside_box = 1;

  const t8dg_mptrac_flux_data *mptrac_flux_data = static_cast<const t8dg_mptrac_flux_data *> (fn_data);

  return mptrac_flux_data->point_is_in_box (x) * value_inside_box;
}

double
t8dg_mptrac_box_indicator_fn (const double x[3], const double t, void *fn_data)
{
  const t8dg_mptrac_flux_data *mptrac_flux_data = static_cast<const t8dg_mptrac_flux_data *> (fn_data);
  double vert_dist, horiz_dist;
  if (mptrac_flux_data->point_is_in_box (x)) {
    return 1;
  }
  mptrac_flux_data->point_max_vert_and_horiz_distance_from_box (x, &vert_dist, &horiz_dist);
  double              radius_vert = 0.1;
  double              radius_horiz = 0.05;
  double              smoothing_factor = 0.5;

  if (vert_dist > (1 + smoothing_factor) * radius_vert 
  || horiz_dist > (1 + smoothing_factor) * radius_horiz) {
    return 0;
  }
  /* Smooth horiz and vert dist */
  if (horiz_dist > radius_horiz);
  horiz_dist = (horiz_dist - radius_horiz) / (radius_horiz * smoothing_factor); /* transform to [0,1] */
  return t8dg_smooth_g (horiz_dist);
}