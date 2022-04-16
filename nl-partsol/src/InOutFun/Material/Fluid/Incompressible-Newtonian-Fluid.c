#include <string.h>
#include "nl-partsol.h"

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

typedef struct {
  bool Is_rho;
  bool Is_Viscosity;

} Check_Material;

static bool Activate_Options(char *);
static void standard_error();
static Check_Material Initialise_Check_Material();
static void check_Compressible_Newtonian_Fluid_Material(Material,
                                                        Check_Material, int);

/***************************************************************************/

Material Define_Incompressible_Newtonian_Fluid(FILE *Simulation_file,
                                               char *Material_Model,
                                               int Material_Idx) {

  int Ndim = NumberDimensions;
  /* Define outputs */
  Material New_Material;

  /* Variables for reading purposes */
  char Parameter_line[MAXC] = {0};
  char *Parameter_pars[MAXW] = {NULL};
  int Parser_status;

  /* Check variables for sintax */
  bool Is_Open = false;
  bool Is_Close = false;
  Check_Material ChkMat = Initialise_Check_Material();

  while (fgets(Parameter_line, sizeof(Parameter_line), Simulation_file) !=
         NULL) {
    /* Parse line */
    Parser_status = parse(Parameter_pars, Parameter_line, delimiters_2);

    if ((strcmp(Parameter_pars[0], "{") == 0) && (Parser_status == 1)) {
      Is_Open = true;
    }
    /**************************************************/
    else if (strcmp(Parameter_pars[0], "rho") == 0) {
      ChkMat.Is_rho = true;
      New_Material.rho = atof(Parameter_pars[1]);
    }
    /**************************************************/
    else if (strcmp(Parameter_pars[0], "Viscosity") == 0) {
      ChkMat.Is_Viscosity = true;
      New_Material.Viscosity = atof(Parameter_pars[1]);
    }
    /**************************************************/
    else if ((strcmp(Parameter_pars[0], "}") == 0) && (Parser_status == 1)) {
      Is_Close = true;
      break;
    }
    /**************************************************/
    else if (Parser_status > 0) {
      sprintf(Error_message, "%s %s", "Undefined", Parameter_pars[0]);
      standard_error();
    }
  }

  strcpy(New_Material.Type, Material_Model);

  check_Compressible_Newtonian_Fluid_Material(New_Material, ChkMat,
                                              Material_Idx);

  /* Return outputs */
  return New_Material;
}

/***************************************************************************/

static Check_Material Initialise_Check_Material() {
  Check_Material ChkMat;

  ChkMat.Is_rho = false; // Reference fensity
  ChkMat.Is_Viscosity = false;

  return ChkMat;
}

/**********************************************************************/

static void check_Compressible_Newtonian_Fluid_Material(Material Mat_particle,
                                                        Check_Material ChkMat,
                                                        int Idx) {
  if (ChkMat.Is_rho && ChkMat.Is_Viscosity) {
    printf("\t -> %s \n", "Incompressible Newtonian Fluid material");
    printf("\t \t -> %s : %i \n", "Idx", Idx);
    printf("\t \t -> %s : %f \n", "Density", Mat_particle.rho);
    printf("\t \t -> %s : %f \n", "Viscosity", Mat_particle.Viscosity);

  } else {
    fprintf(
        stderr, "%s : %s \n", "Error in GramsMaterials()",
        "Some parameter is missed for Incompressible Newtonian Fluid material");
    fputs(ChkMat.Is_rho ? "Density : true \n" : "Density : false \n", stdout);
    fputs(ChkMat.Is_Viscosity ? "Viscosity : true \n" : "Viscosity : false \n",
          stdout);
    exit(EXIT_FAILURE);
  }
}

/**********************************************************************/

static void standard_error() {
  fprintf(stderr, "%s : %s !!! \n", "Error in Define-Material()",
          Error_message);
  exit(EXIT_FAILURE);
}

/***************************************************************************/