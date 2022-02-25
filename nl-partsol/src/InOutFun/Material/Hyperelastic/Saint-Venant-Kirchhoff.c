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

  bool Is_rho;                  // Reference fensity
  bool Is_E;                    // Young modulus
  bool Is_nu;                   // Poisson cefficient
  bool Is_Locking_Control_Fbar; // For incompressible materials
  bool Is_alpha_Fbar;

} Check_Material;

static void standard_error();
static Check_Material Initialise_Check_Material();
static void check_Saint_Venant_Kirchhoff_Material(Material, Check_Material,
                                                  int);
static bool Activate_Options(char *);

/***************************************************************************/

Material Define_Saint_Venant_Kirchhoff(FILE *Simulation_file,
                                       char *Material_Model, int Material_Idx) {

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
    else if (strcmp(Parameter_pars[0], "E") == 0) {
      ChkMat.Is_E = true;
      New_Material.E = atof(Parameter_pars[1]);
    }
    /**************************************************/
    else if (strcmp(Parameter_pars[0], "nu") == 0) {
      ChkMat.Is_nu = true;
      New_Material.nu = atof(Parameter_pars[1]);
    }
    /**************************************************/
    else if (strcmp(Parameter_pars[0], "Fbar") == 0) {
      ChkMat.Is_Locking_Control_Fbar = true;
      New_Material.Locking_Control_Fbar = Activate_Options(Parameter_pars[1]);
    }
    /**************************************************/
    else if (strcmp(Parameter_pars[0], "Fbar-alpha") == 0) {
      ChkMat.Is_alpha_Fbar = true;
      New_Material.alpha_Fbar = atof(Parameter_pars[1]);

      if ((New_Material.alpha_Fbar < 0.0) || (New_Material.alpha_Fbar > 1.0)) {
        sprintf(Error_message, "The range for Fbar-alpha is [0,1]");
        standard_error(Error_message);
      }
    }
    /**************************************************/
    else if ((strcmp(Parameter_pars[0], "}") == 0) && (Parser_status == 1)) {
      Is_Close = true;
      break;
    } else if (Parser_status > 0) {
      sprintf(Error_message, "%s %s", "Undefined", Parameter_pars[0]);
      standard_error();
    }
  }

  strcpy(New_Material.Type, Material_Model);

  check_Saint_Venant_Kirchhoff_Material(New_Material, ChkMat, Material_Idx);

  /* Return outputs */
  return New_Material;
}

/***************************************************************************/

static Check_Material Initialise_Check_Material() {
  Check_Material ChkMat;

  ChkMat.Is_rho = false; // Reference fensity
  ChkMat.Is_E = false;   // Young modulus
  ChkMat.Is_nu = false;  // Poisson cefficient
  ChkMat.Is_Locking_Control_Fbar = false;
  ChkMat.Is_alpha_Fbar = false;

  return ChkMat;
}

/**********************************************************************/

static void check_Saint_Venant_Kirchhoff_Material(Material Mat_particle,
                                                  Check_Material ChkMat,
                                                  int Idx) {
  if (ChkMat.Is_rho && ChkMat.Is_E && ChkMat.Is_nu) {
    printf("\t -> %s \n", "Saint-Venant-Kirchhoff material");
    printf("\t \t -> %s : %i \n", "Idx", Idx);
    printf("\t \t -> %s : %f \n", "Density", Mat_particle.rho);
    printf("\t \t -> %s : %f \n", "Elastic modulus", Mat_particle.E);
    printf("\t \t -> %s : %f \n", "Poisson modulus", Mat_particle.nu);

    if (ChkMat.Is_Locking_Control_Fbar) {
      printf("\t \t -> %s : %s \n", "F-bar", "Enabled");

      if (ChkMat.Is_alpha_Fbar) {
        printf("\t \t -> %s : %f \n", "alpha F-bar", Mat_particle.alpha_Fbar);
      }
    } else {
      printf("\t \t -> %s : %s \n", "F-bar", "Disabled");
    }

  } else {
    fprintf(stderr, "%s : %s \n", "Error in GramsMaterials()",
            "Some parameter is missed for Saint-Venant-Kirchhoff material");
    fputs(ChkMat.Is_rho ? "Density : true \n" : "Density : false \n", stdout);
    fputs(ChkMat.Is_E ? "Elastic modulus : true \n"
                      : "Elastic modulus : false \n",
          stdout);
    fputs(ChkMat.Is_nu ? "Poisson modulus : true \n"
                       : "Poisson modulus : false \n",
          stdout);
    exit(EXIT_FAILURE);
  }
}

/***************************************************************************/

static bool Activate_Options(char *status_text) {
  bool status;

  if (strcmp(status_text, "true") == 0) {
    return true;
  } else if (strcmp(status_text, "false") == 0) {
    return false;
  } else {
    sprintf(Error_message, "The status was %s. Please, use : true/false",
            status_text);
    standard_error();
  }

  return status;
}

/**********************************************************************/

static void standard_error() {
  fprintf(stderr, "%s : %s !!! \n", "Error in Define-Material()",
          Error_message);
  exit(EXIT_FAILURE);
}

/***************************************************************************/