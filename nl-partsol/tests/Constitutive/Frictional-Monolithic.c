/*
Unitary test for the smooth Mohr-Coulomb model.
One single point is initially confined with an
hydrostatic stress state of -200 kPa. Before that,
the stress in the II direction is incresed using
strain control.
-------------------------------------------------------------------
Written by Miguel Molinos, 2021
Universidad Politécnica de Madrid
-------------------------------------------------------------------
Useful information is in Borja et al. (2003):
'On the numerical integration of three-invariant elastoplastic constitutive
models' https://doi.org/10.1016/S0045-7825(02)00620-5
-------------------------------------------------------------------
Requires:
  * accelerate library or LAPACK (solver)
  * gnuplot (plots)
-------------------------------------------------------------------
Compilation recipy (MacOSX):
gcc -framework Accelerate Frictional-Monolithic.c -o Frictional-Monolithic
-------------------------------------------------------------------
Compilation recipy (Linux)
gcc Frictional-Monolithic.c -o Frictional-Monolithic -llapack -lm
-------------------------------------------------------------------
*/

#include <check.h>

/**************************************************************/
/******************** Solver Parameters ***********************/
/**************************************************************/
#define TOL_Radial_Returning 10E-12
#define Max_Iterations_Radial_Returning 10

/**************************************************************/
/******************* Material Parameters **********************/
/**************************************************************/
#define YoungMouduls 100E3
#define PoissonRatio 0.2
#define AtmosphericPressure -100
#define m_Parameter 0.0
#define c0_Parameter 9.0
#define a1_Parameter 20000
#define a2_Parameter 0.005
#define a3_Parameter 35
#define alpha_Parameter -0.71 // -0.71
#define NumberSteps 5000      // 3500 // 2500 //
#define Delta_strain_II -0.00001
#define Yield_Function "Matsuoka-Nakai"

/**************************************************************/
/********************* Required Libraries *********************/
/**************************************************************/

#ifdef __linux__
#include <lapacke.h>
#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <string.h>
#elif __APPLE__
#include <Accelerate/Accelerate.h>
#include <stdbool.h>
#include <stdio.h>
#endif

//#include "CUnit/Basic.h"

#include "../../include/Macros.h"
//#include "../../include/Globals.h"
#include "../../include/Types.h"
#include "../../include/Constitutive.h"

/**************************************************************/
/******************** Auxiliar variables **********************/
/**************************************************************/

typedef struct {
  double alpha;
  double m;
  double pa;
  double c0;
  double a1;
  double a2;
  double a3;
  double CC[9];
} Model_Parameters;

/*
  Define local global variables
*/
int Particle_Idx;
double Error0;
bool Is_Matsuoka_Nakai;
bool Is_Lade_Duncan;
bool Is_Modified_Lade_Duncan;


//START_TEST(test_frictional_monolitic) {
//  State_Parameters Input;
//  State_Parameters Output;
//  Material MatProp;

void main(){
  State_Parameters Input;
  State_Parameters Output;
  Material MatProp;

  // Define material variables
  MatProp.E = YoungMouduls;
  MatProp.nu = PoissonRatio;
  MatProp.atmospheric_pressure = AtmosphericPressure;
  strcpy(MatProp.Yield_Function_Frictional, Yield_Function);
  MatProp.m_Frictional = m_Parameter;
  MatProp.c0_Frictional = c0_Parameter;
  MatProp.a_Hardening_Borja[0] = a1_Parameter;
  MatProp.a_Hardening_Borja[1] = a2_Parameter;
  MatProp.a_Hardening_Borja[2] = a3_Parameter;
  MatProp.alpha_Hardening_Borja = alpha_Parameter;

  // Initialize state variables
  double *stress = (double *)calloc(3 * NumberSteps, sizeof(double));
  double *strain = (double *)calloc(3 * NumberSteps, sizeof(double));
  double *kappa1 = (double *)calloc(NumberSteps, sizeof(double));
  double *Equiv_Plast_Str = (double *)calloc(NumberSteps, sizeof(double));
  double Increment_E_plastic[3] = {0.0, 0.0, 0.0};

  // Set initial stress and strain
  double stress_tr[3] = {0.0, 0.0, 0.0};
  for (int i = 0; i < 3; i++) {
    stress[i] = -200;
    strain[i] = 0.0;
  }

  // Start time integration
  for (int i = 1; i < NumberSteps; i++) {
    // Trial strain
    strain[i * 3 + 1] = strain[(i - 1) * 3 + 1] + Delta_strain_II;

    // Trial stress
    stress[i * 3 + 0] = -200;
    stress[i * 3 + 1] =
        stress[(i - 1) * 3 + 1] + YoungMouduls * Delta_strain_II;
    stress[i * 3 + 2] = -200;

    // Asign variables to the solver
    Input.Particle_Idx = 0;
    Input.Stress = &stress[i * 3];
    Input.Strain = &strain[i * 3];
    Input.Increment_E_plastic = Increment_E_plastic;
    Input.Equiv_Plast_Str = &Equiv_Plast_Str[i - 1];
    Input.Kappa = &kappa1[i - 1];

    // Run solver
    Output = Frictional_Monolithic(Input, MatProp);

    // Update scalar state variables
    Equiv_Plast_Str[i] = *Output.Equiv_Plast_Str;
    kappa1[i] = *Output.Kappa;
  }

  // Save data in a csv file
  FILE *MN_output = fopen("MN_output.csv", "w");
  for (int i = 0; i < NumberSteps; i++) {
    fprintf(MN_output, "%e, %e\n", -strain[i * 3 + 1],
            fabs(stress[i * 3 + 1] - stress[i * 3 + 0]));
  }
  fclose(MN_output);

  // Print data with gnuplot
  FILE *gnuplot = popen("gnuplot -persistent", "w");
  fprintf(gnuplot, "set termoption enhanced \n");
  fprintf(gnuplot, "set datafile separator ',' \n");
  fprintf(gnuplot,
          "set xlabel '- {/Symbol e}_{II}' font 'Times,20' enhanced \n");
  fprintf(gnuplot, "set ylabel 'abs({/Symbol s}_{II} - {/Symbol s}_{I})' font "
                   "'Times,20' enhanced \n");
  fprintf(gnuplot, "plot %s, %s \n",
          "'Borja_compression.csv' title 'Borja et al. (2003)' with line lt "
          "rgb 'blue' lw 2",
          "'MN_output.csv' title 'Us' with line lt rgb 'red' lw 2");
  fflush(gnuplot);

  // Free memory
  free(stress);
  free(strain);
  free(kappa1);
  free(Equiv_Plast_Str);
}
//END_TEST

/**************************************************************/