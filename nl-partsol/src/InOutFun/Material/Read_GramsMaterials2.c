#include <string.h>
#include "nl-partsol.h"

/*
  Call global variables
*/
double TOL_Radial_Returning;
int Max_Iterations_Radial_Returning;

/*
  Local structures
*/
typedef struct {
  int Idx;
  char Model[100];

} Param_Index_and_Model;

/*
  Auxiliar functions and variables
*/

#ifdef _WIN32
static char *delimiters_1 = " ,()\r\n\t";
static char *delimiters_2 = " =\t\r\n";
#else
static char *delimiters_1 = " ,()\n\t";
static char *delimiters_2 = " =\t\n";
#endif
static char *delimiters_3 = "=";

static char Error_message[MAXW];

static Param_Index_and_Model Read_Index_and_Model(char *, char *);
static bool Activate_Options(char *, char *);
static void standard_error();
static void standard_output(char *);
static FILE *Open_and_Check_simulation_file(char *);

/**********************************************************************/

Material *Read_Materials__InOutFun__(char *SimulationFile, int NumberMaterials)
/*

Define-Material(idx=0,Model=Drucker-Prager-Plane-Strain)
{
        E=1
        nu=0.3
        Reference-Plastic-Strain=
        Yield_stress=
        Hardening_modulus=
        Hardening_exponent=
        Cohesion=
        Friction_angle=
        Dilatancy_angle=
}

*/
{
  /* Simulation file */
  FILE *Sim_dat;

  /* Variables for reading purposes */
  char line[MAXC] = {0};
  char *kwords[MAXW] = {NULL};
  int nkwords;

  /* Index for the materials */
  int idx = 0;

  /* Auxiliar parameter */
  Param_Index_and_Model Index_and_Model;

  /* Allocate table with the material */
  Material *List_Materials =
      (Material *)malloc(NumberMaterials * sizeof(Material));
  if (List_Materials == NULL) {
    sprintf(Error_message, "%s", "Memory error for table of material");
    standard_error();
  }

  /* Open and check file */
  Sim_dat = Open_and_Check_simulation_file(SimulationFile);

  while (fgets(line, sizeof line, Sim_dat) != NULL) {

    /* Read the line with the delimiter_1 */
    nkwords = parse(kwords, line, delimiters_1);
    if (nkwords < 0) {
      sprintf(Error_message, "%s", "Parser failed");
      standard_error();
    }

    /* Read Initial-nodal-values */
    if ((nkwords > 0) && (strcmp(kwords[0], "Define-Material") == 0)) {

      /* Read index and model */
      Index_and_Model = Read_Index_and_Model(kwords[1], kwords[2]);

      /* Read material properties  and asign to the material */
      if (strcmp(Index_and_Model.Model, "Solid-Rigid") == 0) {
        List_Materials[idx] = Define_Solid_Rigid(Sim_dat, Index_and_Model.Model,
                                                 Index_and_Model.Idx);
      } else if (strcmp(Index_and_Model.Model, "LE") == 0) {
        List_Materials[idx] = Define_Linear_Elastic(
            Sim_dat, Index_and_Model.Model, Index_and_Model.Idx);
      } else if (strcmp(Index_and_Model.Model, "Saint-Venant-Kirchhoff") == 0) {
        List_Materials[idx] = Define_Saint_Venant_Kirchhoff(
            Sim_dat, Index_and_Model.Model, Index_and_Model.Idx);
      } else if (strcmp(Index_and_Model.Model, "Neo-Hookean-Wriggers") == 0) {
        List_Materials[idx] = Define_Neo_Hookean_Wriggers(
            Sim_dat, Index_and_Model.Model, Index_and_Model.Idx);
      } else if (strcmp(Index_and_Model.Model, "Von-Mises") == 0) {
        List_Materials[idx] = Define_Von_Mises(Sim_dat, Index_and_Model.Model,
                                               Index_and_Model.Idx);
      } else if (strcmp(Index_and_Model.Model, "Granular") == 0) {
        List_Materials[idx] = Define_Frictional(Sim_dat, Index_and_Model.Model,
                                                Index_and_Model.Idx);
      } else if (strcmp(Index_and_Model.Model,
                        "Newtonian-Fluid-Compressible") == 0) {
        List_Materials[idx] = Define_Compressible_Newtonian_Fluid(
            Sim_dat, Index_and_Model.Model, Index_and_Model.Idx);
      } else if (strcmp(Index_and_Model.Model,
                        "Newtonian-Fluid-Incompressible") == 0) {
        List_Materials[idx] = Define_Incompressible_Newtonian_Fluid(
            Sim_dat, Index_and_Model.Model, Index_and_Model.Idx);
      } else {
        sprintf(Error_message, "%s", "Unrecognized kind of material");
        standard_error();
      }

      idx++;
    }
  }

  /* Close  file */
  fclose(Sim_dat);

  return List_Materials;
}

/**********************************************************************/

static Param_Index_and_Model Read_Index_and_Model(char *String_Index,
                                                  char *String_Model) {

  int Num_Idx;
  char *Idx_pars[MAXW] = {NULL};

  int Num_Model;
  char *Model_pars[MAXW] = {NULL};

  /* Define outputs */
  Param_Index_and_Model Parameters;

  Num_Idx = parse(Idx_pars, String_Index, delimiters_3);
  Num_Model = parse(Model_pars, String_Model, delimiters_3);

  /* Fill index of the material */
  if (Num_Idx == 2) {
    Parameters.Idx = atoi(Idx_pars[1]);
  } else if (Num_Idx == 1) {
    Parameters.Idx = atoi(Idx_pars[0]);
  }

  /* Fill model information */
  if (Num_Model == 2) {
    strcpy(Parameters.Model, Model_pars[1]);
  } else if (Num_Model == 1) {
    strcpy(Parameters.Model, Model_pars[0]);
  }

  /* Return outputs */
  return Parameters;
}

/**********************************************************************/

static void standard_error() {
  fprintf(stderr, "%s : %s !!! \n", "Error in Define-Material()",
          Error_message);
  exit(EXIT_FAILURE);
}

/***************************************************************************/

static void standard_output(char *Status_message) {
  fprintf(stdout, "%s \n", Status_message);
}

/**********************************************************************/

static FILE *Open_and_Check_simulation_file(char *Name_File) {
  FILE *Simulation_file = fopen(Name_File, "r");

  if (Simulation_file == NULL) {
    sprintf(Error_message, "%s %s", "Incorrect lecture of", Name_File);
    standard_error();
  }

  return Simulation_file;
}

/***************************************************************************/
