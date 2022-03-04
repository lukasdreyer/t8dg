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
: comm (comm), hours_between_file_reads(hours_between_file_reads), point_source(0, 45, 100, 10, 10, 1)
{
  //const char *mptrac_input = "blub DT_MET 21600 METBASE ei MET_DX 1 MET_DY 1";
  const char *mptrac_input = "blub DT_MET 21600 METBASE wind MET_DX 1 MET_DY 1";
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
  {
    /* TODO: Testing only. Use this, when you need to know the result of the
     *        flux computation (checking conversion formula etc. 
     *        The final result should be (1,1,p) 
     */
    flux_vec[0] = 11111.1111111 * cos (lat/360 * 2*M_PI); /* circumference km/h in [m/s] */
    flux_vec[1] = 5555.55555556; /* 20k */
    flux_vec[2] = 0.277777777778; /* 1000 hPa/h in [hPa/s]  */
  }
#endif
  /* Convert units. u and v are in m/s, which we convert to km/h.
   * w is in hPa/s which we convert to hPa/h */
   /* X m/s = X * 0.001 km/s = X * 3600 * 0.001 km/h */
   static const double sec_per_hour = 3600;
   static const double m_per_km = 1000;
   flux_vec[0] = flux_vec[0] * sec_per_hour / m_per_km;
   flux_vec[1] = flux_vec[1] * sec_per_hour / m_per_km;
   /* X hPa/s is X * 3600 hPa/h */
   flux_vec[2] = flux_vec[2] * sec_per_hour;
   t8dg_debugf ("Flow at %f %f %f = %f %f %f\n", x_vec[0], x_vec[1], x_vec[2], flux_vec[0], flux_vec[1], flux_vec[2]);
   /* Scale down to [0,1]^3 */
   /* So 1 km/h currently actually corresponds to 1 unit cube/hour.
    * We thus need to scale with length of a latitude in u.
    * The length of a longitude in v and the height of the atmosphere in w */
   static const double circumference_equator_km = RE * 2 * M_PI; /* U = 2Pi * radius */
   /* Circumference of latitude is circumference * cos(lat) */
   double circumference_lat_km = circumference_equator_km * cos (lat/360. * 2.* M_PI);
   if (circumference_lat_km <= 0) {
     /* Prevent division by 0. 1m is minimum circumference */
     circumference_lat_km = 1e-3;
   }
   /* Zonal wind (west-east) is scaled with the circumference of the latitude. */
   /* U [Cubelength/h] = U [km/h] / cubelength [km] */
   flux_vec[0] = flux_vec[0] / circumference_lat_km;
   /* Meridional wind (north-south) is scaled with half the circumference. */
   /* V [Cubelength/h] = V [km/h] / cubelength [km] */
   flux_vec[1] = flux_vec[1] / (circumference_equator_km / 2);

  /* To compute the speed in z direction we need to convert the pressure
   * gradient in hPa/h to km/h.
   * This is not the same as converting hPa to km.
   * 
   * We have: P[hPa] + v[hPa/h] = P_1[hPa]
   *  for the pressure P_1 in 1 hour in hPa.
   * we need v_x such that
   *          P_x[] + v_x[] = P_1[]
   * for the pressure in [0,1] coordinates.
   * With F    = t8_mptrac_coords_to_lonlatpressure
   * and  F^-1 = t8_mptrac_pressure_to_coord
   * we have P[hPa] = F(P_x[]) and get
   *  v_x[] = F^-1(P[hPa] + v[hPa/h]) - P_x[]
   */
  /* TODO: Use DP2DZ */
  double pressure_temp;
  /* Compute F^-1(P[hPa] + v[hPa/h]) */
  t8_mptrac_pressure_to_coord (mptrac_context, pressure + flux_vec[2], &pressure_temp);
  /* Subtract P_x */
  flux_vec[2] = pressure_temp - x_vec[2];
#if T8_ENABLE_DEBUG
  /* Check that F(F^-1) = id */
  double pressure_inv;
  t8_mptrac_pressure_to_coord (mptrac_context, pressure, &pressure_inv);
  T8_ASSERT (fabs (pressure_inv - x_vec[2]) < 1e-5);
  double x_vec_inv[3] = {x_vec[0], x_vec[1], pressure_inv};
  double pressure_inv_orig;
  t8_mptrac_coords_to_lonlatpressure (mptrac_context, x_vec_inv, &lon, &lat, &pressure_inv_orig);
  T8_ASSERT (fabs (pressure_inv_orig - pressure) < 1e-5);
#endif

#ifdef T8_ENABLE_DEBUG
   static const int           max_p_idx = meteo1->np;
   static const double        pressure_min_in_km = Z (meteo1->p[0]);
   static const double        pressure_max_in_km = Z (meteo1->p[max_p_idx - 1]);
   static const double        height_of_atmosphere = pressure_max_in_km - pressure_min_in_km;
#if 0
   flux_vec[2] = flux_vec[2] / height_of_atmosphere;
#endif
   t8dg_debugf ("Flow at %f %f %f = %f %f %f\n", x_vec[0], x_vec[1], x_vec[2], flux_vec[0], flux_vec[1], flux_vec[2]);
   t8dg_debugf ("height of atmosphere [km] = %f\n", height_of_atmosphere);
   t8dg_debugf ("circumference of lat [km] = %f\n", circumference_lat_km);
   t8dg_debugf ("lon/lat coords %f %f\n", lon, lat);
#endif
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
 
  double vert_dist_from_center, horiz_dist_from_center;
  mptrac_flux_data->point_max_vert_and_horiz_distance_from_box (x, &vert_dist_from_center, &horiz_dist_from_center);
  const double              radius_horiz_deg = 1;
  const double              radius_vert_km = 1;
  const double              smoothing_factor_horiz = 0;//0.2; //0.5; TODO: values != 0 result in kind of weird shape
  const double              smoothing_factor_vert = 0;

  if (vert_dist_from_center < radius_vert_km
   && horiz_dist_from_center < radius_horiz_deg) {
    return 1;
  }

  if (horiz_dist_from_center > (1 + smoothing_factor_horiz) * radius_horiz_deg ||
    vert_dist_from_center > (1 + smoothing_factor_vert) * radius_vert_km) {
    return 0;
  }
  /* transform to [0,1] */
  horiz_dist_from_center = smoothing_factor_horiz == 0 ? 0 : (horiz_dist_from_center - radius_horiz_deg) / (radius_horiz_deg * smoothing_factor_horiz); 
  vert_dist_from_center = smoothing_factor_vert == 0 ? 0 : (vert_dist_from_center - radius_vert_km) / (radius_vert_km * smoothing_factor_vert);
  const double dist = sqrt(horiz_dist_from_center*horiz_dist_from_center + vert_dist_from_center*vert_dist_from_center);
  /* Smooth horiz and vert dist */
  return t8dg_smooth_g (dist);
}
