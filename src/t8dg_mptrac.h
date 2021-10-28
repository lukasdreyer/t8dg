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

class t8dg_mptrac_flux_data : public t8dg_flux_data_base
{
    public:

    t8dg_mptrac_flux_data (const char *nc_filename);

    ~t8dg_mptrac_flux_data ();

    void initialize (const struct t8dg_linear_advection_diffusion_problem *problem);

    void start_new_time_step (const struct t8dg_linear_advection_diffusion_problem *problem);


    const t8_mptrac_context_t *get_context() const;

    protected:
        t8_mptrac_context_t * context;
};

#endif /* SRC_T8DG_MPTRAC_H_ */
