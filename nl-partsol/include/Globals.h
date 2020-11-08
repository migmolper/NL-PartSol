
#define MAXW 100
#define MAXC 1000
#define NumberDimensions 2
#define TOL_InOut 10E-23
#define TOL_NR 10E-6
#define TOL_zero 10E-10

extern char * SimulationFile;
extern char * RestartFile;
extern char * FEM_MeshFileName;
extern char * MPM_MeshFileName;
extern char * TimeIntegrationScheme;
extern char * ShapeFunctionGP;
extern char * Formulation;
extern char   OutputParticlesFile[MAXC];
extern char   OutputNodesFile[MAXC];
extern char   OutputDir[MAXC];
extern char   Field[10];

extern double gamma_LME;
extern double TOL_lambda;
extern double CFL; /* Courant number (0-1) */
extern double DeltaTimeStep;
extern double Error0;
extern double SpectralRadius;
extern double epsilon_Mass_Matrix;
extern double beta_Newmark_beta;
extern double gamma_Newmark_beta;
extern double TOL_Newmark_beta;

extern int NumTimeStep;
extern int ResultsTimeStep;
extern int NumberDOF;

extern bool Out_global_coordinates;
extern bool Out_element_coordinates;
extern bool Out_mass;
extern bool Out_density;
extern bool Out_damage;
extern bool Out_nodal_idx;
extern bool Out_material_idx;
extern bool Out_velocity;
extern bool Out_acceleration;
extern bool Out_displacement;
extern bool Out_stress;
extern bool Out_eigenvalues_stress;
extern bool Out_volumetric_stress;
extern bool Out_strain;
extern bool Out_eigenvalues_strain;
extern bool Out_deformation_gradient;
extern bool Out_energy;