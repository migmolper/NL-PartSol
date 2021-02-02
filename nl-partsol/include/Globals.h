
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


/* Convergence parameters for radia returning algorithm */
extern double TOL_Radial_Returning;
extern int Max_Iterations_Radial_Returning;


extern int NumTimeStep;
extern int ResultsTimeStep;
extern int NumberDOF;


/*
	Soil-water mixture
*/
extern int Number_Soil_Water_Mixtures;
extern Mixture * Soil_Water_Mixtures;

/*
	Variables for the outputs 
*/

extern Event * Out_nodal_path_csv;
extern int Number_Out_nodal_path_csv;

extern Event * Out_particles_path_csv;
extern int Number_Out_particles_path_csv;

extern Event * Out_Gauss_Point_evolution_csv;
extern int Number_Out_Gauss_Point_evolution_csv;

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
extern bool Out_Pw;
extern bool Out_dPw_dt;
extern bool Out_strain;
extern bool Out_eigenvalues_strain;
extern bool Out_deformation_gradient;
extern bool Out_energy;
extern bool Out_Von_Mises;
extern bool Out_EPS;

