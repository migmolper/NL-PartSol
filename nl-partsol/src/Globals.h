/**
 * @file Globals.h
 * @author Miguel Molinos (@migmolper)
 * @brief 
 * @version 0.1
 * @date 2022-05-25
 * 
 * @copyright Copyright (c) 2022
 * 
 */

#ifndef _GLOBALS_H_
#define _GLOBALS_H_

extern char   SimulationFile[MAXC];
extern char   Static_conditons[MAXC];
extern char   Formulation[MAXC];
extern char * FEM_MeshFileName;
extern char * MPM_MeshFileName;
extern char * TimeIntegrationScheme;
extern char   ShapeFunctionGP[MAXC];
extern char   OutputParticlesFile[MAXC];
extern char   OutputNodesFile[MAXC];
extern char   OutputDir[MAXC];
extern char   Field[10];

//! Select the type of solver 
extern bool Petsc_Direct_solver;
extern bool Petsc_Iterative_solver;

//	Global variables for simulation purposes
extern double DeltaTimeStep;
extern int NumberDOF;

// Gravity field
extern Load gravity_field;

// Use Fracture modul
extern bool Driver_EigenErosion;
extern bool Driver_EigenSoftening;

// Material parameters
extern int Number_Soil_Water_Mixtures;
extern Mixture * Soil_Water_Mixtures;

// Global variables for the LME shape functions
extern int max_iter_LME;
extern double TOL_zero_LME;
extern double TOL_wrapper_LME;
extern double gamma_LME;
extern char   wrapper_LME[MAXC];

// Parameter for plain stress simulations
extern double Thickness_Plain_Stress;

/* Convergence parameters for radial returning algorithm */
extern double TOL_Radial_Returning;
extern int Max_Iterations_Radial_Returning;

/*
	Variables for the outputs 
*/
extern bool Flag_Print_Convergence;

extern int ResultsTimeStep;

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
extern bool Out_water_pressure;
extern bool Out_Pw;
extern bool Out_dPw_dt;
extern bool Out_strain;
extern bool Out_eigenvalues_strain;
extern bool Out_deformation_gradient;
extern bool Out_green_lagrange;
extern bool Out_plastic_deformation_gradient;
extern bool Out_plastic_jacobian;
extern bool Out_Metric;
extern bool Out_energy;
extern bool Out_Von_Mises;
extern bool Out_EPS;
extern bool Out_Partition_Unity;


/* Variables for the backups */
extern char Output_Backup_File[MAXC];
extern bool Backup_damage;
extern bool Backup_plastic_deformation_gradient;
extern bool Backup_EPS;

#endif