// clang-format off
#include <stdlib.h>
#include <stdbool.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "Macros.h"
#include "Types.h"
#include "Globals.h"
#include "Matlib.h"
#include "Particles.h"
#include "Constitutive/Constitutive.h"
#include "InOutFun.h"
// clang-format on

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

  bool Is_Ceps; // Normalizing constant (Eigenerosion)
  bool Is_Gf;  // Failure energy (Eigenerosion)

} Check_Material;

static void standard_error();
static Check_Material __Initialise_Check_Material();
static int __check_material(Material *NH_Material, Check_Material ChkMat,
                            int Idx);
static bool Activate_Options(char *);

/***************************************************************************/

int Define_Neo_Hookean_Wriggers(Material *NH_Material, FILE *Simulation_file,
                                     char *Material_Model, int Material_Idx) {

  int STATUS = EXIT_SUCCESS;
  int Ndim = NumberDimensions;

  /* Variables for reading purposes */
  char Parameter_line[MAXC] = {0};
  char *Parameter_pars[MAXW] = {NULL};
  int Parser_status;

  /* Check variables for sintax */
  bool Is_Open = false;
  bool Is_Close = false;
  Check_Material ChkMat = __Initialise_Check_Material();

  /* Default parameters */
  (*NH_Material).Locking_Control_Fbar = false;
  (*NH_Material).alpha_Fbar = 0.0;

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
      (*NH_Material).rho = atof(Parameter_pars[1]);
    }
    /**************************************************/
    else if (strcmp(Parameter_pars[0], "E") == 0) {
      ChkMat.Is_E = true;
      (*NH_Material).E = atof(Parameter_pars[1]);
    }
    /**************************************************/
    else if (strcmp(Parameter_pars[0], "nu") == 0) {
      ChkMat.Is_nu = true;
      (*NH_Material).nu = atof(Parameter_pars[1]);
    }
    /**************************************************/
    else if (strcmp(Parameter_pars[0], "Fbar") == 0) {
      ChkMat.Is_Locking_Control_Fbar = true;
      (*NH_Material).Locking_Control_Fbar = Activate_Options(Parameter_pars[1]);
    }
    /**************************************************/
    else if (strcmp(Parameter_pars[0], "Fbar-alpha") == 0) {
      ChkMat.Is_alpha_Fbar = true;
      (*NH_Material).alpha_Fbar = atof(Parameter_pars[1]);

      if (((*NH_Material).alpha_Fbar < 0.0) || ((*NH_Material).alpha_Fbar > 1.0)) {
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

  strcpy((*NH_Material).Type, Material_Model);

  STATUS = __check_material(NH_Material, ChkMat, Material_Idx);
  if (STATUS == EXIT_FAILURE) {
    fprintf(stderr, "" RED " Error in __check_material() " RESET " \n");
    return EXIT_FAILURE;
  }

  /* Return outputs */
  return STATUS;
}

/***************************************************************************/

static Check_Material __Initialise_Check_Material() {
  Check_Material ChkMat;

  ChkMat.Is_rho = false; 
  ChkMat.Is_E = false;
  ChkMat.Is_nu = false; 
  ChkMat.Is_Locking_Control_Fbar = false;
  ChkMat.Is_alpha_Fbar = false;
  ChkMat.Is_Ceps = true;
  ChkMat.Is_Gf = true;

  return ChkMat;
}

/**********************************************************************/

static int __check_material(Material * NH_Material,
                                                Check_Material ChkMat,
                                                int Idx) {

  int STATUS = EXIT_SUCCESS;

  if (ChkMat.Is_rho && ChkMat.Is_E && ChkMat.Is_nu) {

    printf("\t -> %s \n", "Neo-Hookean material");

    printf("\t \t -> %s : %f \n", "" MAGENTA "[rho]" RESET "",
           (*NH_Material).rho);

    printf("\t \t -> %s : %f \n", "" MAGENTA "[E]" RESET "", (*NH_Material).E);

    printf("\t \t -> %s : %f \n", "" MAGENTA "[nu]" RESET "", (*NH_Material).nu);

    if ((Driver_EigenErosion == true) || (Driver_EigenSoftening == true)) {
      printf("\t \t -> %s : %f \n", "" MAGENTA "[Ceps]" RESET "",
             (*NH_Material).Ceps);

      printf("\t \t -> %s : %f \n", "" MAGENTA "[Gf]" RESET "",
             (*NH_Material).Gf);
    }

    if (ChkMat.Is_Locking_Control_Fbar) {
      printf("\t \t -> %s : %s \n", "F-bar", "Enabled");

      if (ChkMat.Is_alpha_Fbar) {
        printf("\t \t -> %s : %f \n", "alpha F-bar", (*NH_Material).alpha_Fbar);
      }
    } else {
      printf("\t \t -> %s : %s \n", "F-bar", "Disabled");
    }

  } else {

    fprintf(stderr, "" RED " %s : %s " RESET " \n", "Error in GramsMaterials()",
            "Some parameter is missed for Neo-Hookean");

    fputs(ChkMat.Is_rho
              ? "" MAGENTA "[rho]" RESET " : " GREEN "true" RESET " \n"
              : "" MAGENTA "[rho]" RESET " : " RED "false" RESET " \n",
          stderr);

    fputs(ChkMat.Is_E ? "" MAGENTA "[E]" RESET " : " GREEN "true" RESET " \n"
                      : "" MAGENTA "[E]" RESET " : " RED "false" RESET " \n",
          stderr);

    fputs(ChkMat.Is_nu ? "" MAGENTA "[nu]" RESET " : " GREEN "true" RESET " \n"
                       : "" MAGENTA "[nu]" RESET " : " RED "false" RESET " \n",
          stderr);

    if ((Driver_EigenErosion == true) || (Driver_EigenSoftening == true)) {
      fputs(ChkMat.Is_Ceps
                ? "" MAGENTA "[Ceps]" RESET " : " GREEN "true" RESET " \n"
                : "" MAGENTA "[Ceps]" RESET " : " RED "false" RESET " \n",
            stderr);

      fputs(ChkMat.Is_Gf
                ? "" MAGENTA "[Gf]" RESET " : " GREEN "true" RESET " \n"
                : "" MAGENTA "[Gf]" RESET " : " RED "false" RESET " \n",
            stderr);
    }

    STATUS = EXIT_FAILURE;
  }

  return STATUS;

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