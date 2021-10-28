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


/* Source function for a box source */
double t8dg_mptrac_box_source (const double x[3], const double t, void *fn_data);

class t8dg_mptrac_flux_data : public t8dg_flux_data_base
{
    public:

    t8dg_mptrac_flux_data (const char *nc_filename, int hours_between_file_reads);

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

    protected:
    t8_mptrac_context_t * context;

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
};

#endif /* SRC_T8DG_MPTRAC_H_ */
