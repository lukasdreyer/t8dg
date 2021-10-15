#ifndef SRC_T8DG_FLUX_IMPLEMENTATION_H_
#define SRC_T8DG_FLUX_IMPLEMENTATION_H_

T8DG_EXTERN_C_BEGIN ();

typedef struct t8dg_linear_flux3D_constant_flux_data
{
  double              flow_direction[3];
  double              flow_velocity;
} t8dg_linear_flux3D_constant_flux_data_t;

void                t8dg_linear_flux3D_constant_flux_fn (double x_vec[3], double flux_vec[3], double t, void *flux_data);

void                t8dg_rotating_flux_2D_fn (double x_vec[3], double flux_vec[3], double t, void *flux_data);

void                t8dg_spiral_flux_3D_fn (double x_vec[3], double flux_vec[3], double t, void *flux_data);

void                t8dg_mptrac_flow_3D_fn (double x_vec[3], double flux_vec[3], double t, void *flux_data);

typedef struct t8dg_linear_flux_1D_from_3D_data
{
  t8dg_linear_flux3D_fn flux_3D_fn;
  int                 component;
  void               *flux_3D_data;
} t8dg_linear_flux_1D_from_3D_data_t;

void                t8dg_linear_flux_1D_from_3D_fn (double x_vec[3], double *flux_value, double t, void *flux_data);

double              t8dg_linear_numerical_flux3D_lax_friedrich_fn (const double u_minus, const double u_plus, const double flow_vector[3],
                                                                   const double normal_vector[3], const void *numerical_flux_data);

double              t8dg_numerical_flux1D_central (const double u_minus, const double u_plus,
                                                   const double normal_component, const void *numerical_flux_data);
double
 
 
 
 
     t8dg_numerical_flux1D_left (const double u_minus, const double u_plus, const double normal_component, const void *numerical_flux_data);

double
 
 
 
 
 
   t8dg_numerical_flux1D_right (const double u_minus, const double u_plus, const double normal_component, const void *numerical_flux_data);

T8DG_EXTERN_C_END ();

#endif
