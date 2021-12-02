/*
 * t8dg_mptrac.h
 *
 *  Created on: October 15, 2021
 *      Author: Johannes Holke
 */

/** @file t8dg_mptrac.h */
#ifndef SRC_T8DG_MPTRAC_H_
#define SRC_T8DG_MPTRAC_H_

#include <t8/example/mptrac/t8_mptrac_interpolate.h>
#include <t8dg.h>
#include <t8dg_flux_data_base.h>

/** Defines a box via a midpoint in lon,lat,pressure coordinates
 *  and an extend given in km.
 */
class t8dg_mptrac_box
{ 
  public:
  /*< Construct a new box with given center and extend. */
    t8dg_mptrac_box (double center_lon, double center_lat, double center_p, double extend_lon_deg, double extend_lat_deg, double extend_p_km);

    /*< Query whether a given lon/lat/pressure point is inside the box. */
    int lon_lat_pressure_is_in (const double lon, const double lat, const double pressure) const;

    /*< Determine the horizontal distance of a point from the center in maximum norm in degrees. */
    double max_horiz_distance_from_center (const double lon, const double lat) const;

    /*< Determine the vertical distance of a point from the center in maximum norm in km. */
    double vert_distance_from_center (const double pressure_km) const;

  protected:
    const double box_center[3]; /* units are deg, deg, hpa */
    const double box_extend[3]; /* units are deg, deg, km */
    double box_minimum[3];
    double box_maximum[3];
    double pressure_center_in_km;
};

/* Source function for a box source */
double t8dg_mptrac_box_source (const double x[3], const double t, void *fn_data);

/* Indicator function for a box source. */
double t8dg_mptrac_box_indicator_fn (const double x[3], const double t, void *fn_data);

class t8dg_mptrac_flux_data : public t8dg_flux_data_base
{
    public:

    t8dg_mptrac_flux_data (const char *nc_filename, int hours_between_file_reads, sc_MPI_Comm comm);

    ~t8dg_mptrac_flux_data ();

    void initialize (const struct t8dg_linear_advection_diffusion_problem *problem);

    void start_new_time_step (const struct t8dg_linear_advection_diffusion_problem *problem);

    /* Before we enter an element for the first time, we check on which boundary
     * it lies so that we can adapt the flux later. */
    void before_first_call_on_element (t8_forest_t forest, t8_locidx_t itree, t8_locidx_t element_in_tree);


    /** Get the current physical time. */
    inline const double get_time () const;

    inline const t8_mptrac_context_t *get_context() const;

    inline bool is_current_element_at_pole () const;
    
    inline bool is_current_element_at_top_or_bottom () const;

    inline int point_is_in_box (const double x[3]) const;

    /** Determin the maximum vertical distance in km and horizontal distance in degrees
     *  of a point x in [0,1]^3 from the box point source. */
    void point_max_vert_and_horiz_distance_from_box (const double x[3], double *vert_dist, double *horiz_dist) const;

    protected:
    t8_mptrac_context_t * context;
    sc_MPI_Comm comm; /*< MPI Communicator used for \a context. */

    /* How many simulation hours after 05.06.2011 00:00 to start. */
    double start_six_hours;
    /* Current simulation time in seconds */
    double physical_time_s;

    /* How many simulation hours have passed since last file read.
    * Every 6 hours we read a new file. */
    double hours_since_last_file_read;

    /* How many hours to pass before a new input file is read. */
    double hours_between_file_reads;

    bool current_element_is_at_pole;
    bool current_element_is_at_top_or_bottom;

    t8dg_mptrac_box point_source;
};

#endif /* SRC_T8DG_MPTRAC_H_ */
