#ifndef SRC_T8DG_FLUX_DATA_BASE_H_
#define SRC_T8DG_FLUX_DATA_BASE_H_

#include <t8dg.h>

/* Forward declaration to prevent circle includes */
struct t8dg_linear_advection_diffusion_problem;


class t8dg_flux_data_base
{
public:
  t8dg_flux_data_base () {}
   /** The destructor. It does nothing but has to be defined since
   * we may want to delete geometry that is actually inherited
   * and providing an implementation
   * for the destructor ensures that the
   * destructor of the child class will be executed. */
  virtual ~t8dg_flux_data_base () {
  }

  /* The following functions can be used to update the flux data during the 
   * computation.
   * In the base class these are all empty. */

  /** Initialize the data. This function is called after the advection diffusion problem
   * was setup.
   * \param [in] problem The advection diffusion problem.
   */
  virtual void initialize (const struct t8dg_linear_advection_diffusion_problem *problem) {};

  /** This function is called at the start of a new time step.
   * \param [in] problem The advection diffusion problem.
   */
  virtual void start_new_time_step (const struct t8dg_linear_advection_diffusion_problem *problem) {};


  virtual void before_first_call_on_element (t8_forest_t forest, t8_locidx_t itree, t8_locidx_t element_in_tree) {};
  virtual void after_last_call_on_element (t8_forest_t forest, t8_locidx_t itree, t8_locidx_t element_in_tree) {};
protected:
};

class t8dg_linear_flux3D_constant_flux_data : public t8dg_flux_data_base
{
  public:

  t8dg_linear_flux3D_constant_flux_data (double flow_x, double flow_y, double flow_z, double flow_velocity_in)
  : flow_direction{flow_x, flow_y, flow_z}, flow_velocity(flow_velocity_in)
  {
  }

  const double* get_flow_direction () const {
    return flow_direction;
  }
  const double get_flow_velocity () const {
    return flow_velocity;
  }
  protected:
  double              flow_direction[3];
  double              flow_velocity;
};

#endif
