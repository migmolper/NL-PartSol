// Global libs
#include "Macros.h"
#include "Types.h"
#include "Globals.h"
#include "Matlib.h"
#include "Particles.h"

typedef struct {

#ifdef USE_PETSC
  Vec value;
  Vec d_value_dt;
  Vec d2_value_dt2;
#else
  double *value;
  double *d_value_dt;
  double *d2_value_dt2;
#endif

} Nodal_Field;

typedef struct {
  double alpha_1;
  double alpha_2;
  double alpha_3;
  double alpha_4;
  double alpha_5;
  double alpha_6;
} Newmark_parameters;



/*!
  \brief Finite strains Newmark-beta

  \param[in] beta: First Newmark-beta parameter
  \param[in] gamma: Second Newmark-beta parameter
  \param[in] DeltaTimeStep: Timestep increment
*/
static Newmark_parameters __compute_Newmark_parameters(double beta,
                                                       double gamma,
                                                       double DeltaTimeStep);
/**************************************************************/

static Matrix __compute_nodal_lumped_mass(Particle, Mesh, Mask, double);



static void compute_Gravity_field(Mask, Particle, int, int);
static Nodal_Field __get_nodal_field_tn(Particle, Mesh, Mask);
static Nodal_Field __initialise_nodal_incrementss(Nodal_Field, Mask, Mesh,
                                               Newmark_parameters, int, int);
static void update_Local_State(Nodal_Field, Mask, Particle, Mesh);
static Matrix __assemble_residual(Nodal_Field, Nodal_Field, Mask, Particle,
                               Mesh, Newmark_parameters, int, int);
static void compute_Inertial_Forces_Mixture(Nodal_Field, Nodal_Field, Matrix,
                                            Matrix, Mask, Particle,
                                            Newmark_parameters);
static void compute_Internal_Forces_Mixture(Matrix, Mask, Particle, Mesh);
static Tensor compute_total_first_Piola_Kirchhoff_stress(Tensor, double,
                                                         Tensor);
static void compute_Rate_Mass_Fluid(Matrix, Mask, Particle, Mesh);
static void compute_Flow_contribution_Fluid(Nodal_Field, Nodal_Field, Matrix,
                                            Mask, Particle, Mesh);
static Tensor compute_Kirchoff_Pore_water_pressure_gradient_n1(Matrix, Matrix,
                                                               Matrix);
static void compute_nominal_traction_and_fluid_flux(Matrix, Mask, Particle,
                                                    Mesh, int, int);
static Matrix compute_Nodal_Reactions(Mesh, Matrix, Mask, int, int);
static bool check_convergence(Matrix, double, int, int, int);

static Matrix __assemble_tangent_stiffness(Nodal_Field, Nodal_Field, Mask,
                                         Particle, Mesh, double,
                                         Newmark_parameters);

static Matrix compute_mixture_stiffness_density(Tensor, Tensor, Tensor, double,
                                                double);
static Matrix compute_mixture_inertial_density(Tensor, Tensor, double, double,
                                               double, double, double, double,
                                               double, double, double, double,
                                               double);
static Matrix compute_water_flux_density(Tensor, Tensor, Tensor, Tensor, Tensor,
                                         double, double, double, double, double,
                                         double);
static Matrix compute_water_inertial_density(Tensor, Tensor, double, double,
                                             double, double, double, double,
                                             double, double, double, double,
                                             double);

static void system_reduction(Matrix, Matrix, Mask, Mesh, int, int);
static void solve_system(Nodal_Field, Matrix, Matrix);

static void __update_Nodal_Increments(Nodal_Field, Nodal_Field,
                                            Newmark_parameters);
static void __update_Particles(Nodal_Field, Particle, Mesh, Mask);
static void particle_results_vtk__InOutFun__(Nodal_Field, Nodal_Field, Particle, Mesh, Mask, int,
                            int);
static char Error_message[MAXW];
static void standard_error();