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

typedef struct {

  bool Is_rho;                      // Reference fensity
  bool Is_E;                        // Young modulus
  bool Is_nu;                       // Poisson cefficient
  bool Is_Hardening_Exponent;       // Hardening exponent
  bool Is_Hardening_modulus;        // Hardening modulus
  bool Is_Reference_pressure;       // Reference pressure
  bool Is_Reference_plastic_strain; // Reference Plastic Strain
  bool Is_kappa_0;                  // Initial yield
  bool Is_friction_angle;           // Friction angle
  bool Is_dilatancy_angle;          // Dilatancy angle
  bool Is_J2_degradated;            // Degradation limit
  bool Is_Ceps;                     // Normalizing constant (Eigenerosion)
  bool Is_Gf;                       // Failure energy (Eigenerosion)
  bool Is_ft;                       // Tensile strengt of the material
  bool Is_heps;                     // Bandwidth of the cohesive fracture
  bool Is_wcrit;                    // Critical opening displacement
} Check_Material;

static Check_Material __Initialise_Check_Material();
static int __check_material(Material *, Check_Material, int);
static void __print_miss_variables(Check_Material ChkMat);
static int __read_boolean(bool *Value, char *Value_char);

/**********************************************************************/

int Define_Drucker_Prager(Material *DP_Material, FILE *Simulation_file,
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
  (*DP_Material).kappa_0 = 0.0;
  (*DP_Material).ReferencePressure = 0.0;
  (*DP_Material).J2_degradated = 0.0;
  TOL_Radial_Returning = 1E-14;
  Max_Iterations_Radial_Returning = 10;

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
      (*DP_Material).rho = atof(Parameter_pars[1]);
    }
    /**************************************************/
    else if (strcmp(Parameter_pars[0], "E") == 0) {
      ChkMat.Is_E = true;
      (*DP_Material).E = atof(Parameter_pars[1]);
    }
    /**************************************************/
    else if (strcmp(Parameter_pars[0], "nu") == 0) {
      ChkMat.Is_nu = true;
      (*DP_Material).nu = atof(Parameter_pars[1]);
    }
    /**************************************************/
    else if (strcmp(Parameter_pars[0], "m") == 0) {
      ChkMat.Is_Hardening_Exponent = true;
      (*DP_Material).Exponent_Hardening_Ortiz = atof(Parameter_pars[1]);
    }
    /**************************************************/
    else if (strcmp(Parameter_pars[0], "Hardening-modulus") == 0) {
      ChkMat.Is_Hardening_modulus = true;
      (*DP_Material).Hardening_modulus = atof(Parameter_pars[1]);
    }
    /**************************************************/
    else if (strcmp(Parameter_pars[0], "Reference-pressure") == 0) {
      ChkMat.Is_Reference_pressure = true;
      (*DP_Material).ReferencePressure = atof(Parameter_pars[1]);
    }
    /**************************************************/
    else if (strcmp(Parameter_pars[0], "Reference-plastic-strain") == 0) {
      ChkMat.Is_Reference_plastic_strain = true;
      (*DP_Material).Plastic_Strain_0 = atof(Parameter_pars[1]);
    }
    /**************************************************/
    else if (strcmp(Parameter_pars[0], "kappa-0") == 0) {
      ChkMat.Is_kappa_0 = true;
      (*DP_Material).kappa_0 = atof(Parameter_pars[1]);
    }
    /**************************************************/
    else if (strcmp(Parameter_pars[0], "Friction-angle") == 0) {
      ChkMat.Is_friction_angle = true;
      (*DP_Material).phi_Frictional = atof(Parameter_pars[1]);

      if (((*DP_Material).phi_Frictional > 90.0) &&
          ((*DP_Material).phi_Frictional < 0.0)) {
        fprintf(stderr,
                "" RED " Invalid value of the [Friction-angle] " RESET " \n");
        return EXIT_FAILURE;
      }
    }
    /**************************************************/
    else if (strcmp(Parameter_pars[0], "Dilatancy-angle") == 0) {
      ChkMat.Is_dilatancy_angle = true;
      (*DP_Material).psi_Frictional = atof(Parameter_pars[1]);
    }
    /**************************************************/
    else if (strcmp(Parameter_pars[0], "J2-degradated") == 0) {
      ChkMat.Is_J2_degradated = true;
      (*DP_Material).J2_degradated = atof(Parameter_pars[1]);

      if (((*DP_Material).psi_Frictional > 90.0) &&
          ((*DP_Material).psi_Frictional < 0.0)) {
        fprintf(stderr,
                "" RED " Invalid value of the [Dilatancy-angle] " RESET " \n");
        return EXIT_FAILURE;
      }
    }
    /**************************************************/
    else if (strcmp(Parameter_pars[0], "Ceps") == 0) {
      if ((Driver_EigenErosion == true) || (Driver_EigenSoftening == true)) {
        ChkMat.Is_Ceps = true;
        (*DP_Material).Ceps = atof(Parameter_pars[1]);
      }
    }
    /**************************************************/
    else if (strcmp(Parameter_pars[0], "Gf") == 0) {
      if (Driver_EigenErosion == true) {
        ChkMat.Is_Gf = true;
        (*DP_Material).Gf = atof(Parameter_pars[1]);
      }
    }
    /**************************************************/
    else if (strcmp(Parameter_pars[0], "ft") == 0) {
      if (Driver_EigenSoftening == true) {
        ChkMat.Is_ft = true;
        (*DP_Material).ft = atof(Parameter_pars[1]);
      }
    }
    /**************************************************/
    else if (strcmp(Parameter_pars[0], "heps") == 0) {
      if (Driver_EigenSoftening == true) {
        ChkMat.Is_heps = true;
        (*DP_Material).heps = atof(Parameter_pars[1]);
      }
    }
    /**************************************************/
    else if (strcmp(Parameter_pars[0], "wcrit") == 0) {
      if (Driver_EigenSoftening == true) {
        ChkMat.Is_wcrit = true;
        (*DP_Material).wcrit = atof(Parameter_pars[1]);
      }
    }
    /**************************************************/
    else if ((strcmp(Parameter_pars[0], "}") == 0) && (Parser_status == 1)) {
      Is_Close = true;
      break;
    } else if (Parser_status > 0) {
      fprintf(stderr, RED "Undefined %s \n" RESET, Parameter_pars[0]);
      return EXIT_FAILURE;
    }
  }

  strcpy((*DP_Material).Type, Material_Model);

  double kappa0 = (*DP_Material).kappa_0;
  double m = (*DP_Material).Exponent_Hardening_Ortiz;
  double H = (*DP_Material).Hardening_modulus;

  if (((m > 0) && (H > 0)) || ((m < 0) && (H < 0))) {
    if (ChkMat.Is_Reference_plastic_strain == false) {
      (*DP_Material).Plastic_Strain_0 =
          (kappa0 / (m * H)) * pow(1, (1.0 / m - 1.0));
    }
  } else {
    fprintf(stderr,
            "" RED " [Hardening-modulus] and [m] must have the same sign " RESET
            " \n");
    return EXIT_FAILURE;
  }

  STATUS = __check_material(DP_Material, ChkMat, Material_Idx);
  if (STATUS == EXIT_FAILURE) {
    fprintf(stderr, "" RED " Error in " RESET "" BOLDRED
                    "__check_material() " RESET " \n");
    __print_miss_variables(ChkMat);
    return EXIT_FAILURE;
  }

  return STATUS;
}

/***************************************************************************/

static Check_Material __Initialise_Check_Material() {
  Check_Material ChkMat;

  ChkMat.Is_rho = false;
  ChkMat.Is_E = false;
  ChkMat.Is_nu = false;
  ChkMat.Is_Hardening_Exponent = false;
  ChkMat.Is_Hardening_modulus = false;
  ChkMat.Is_Reference_pressure = false;
  ChkMat.Is_Reference_plastic_strain = false;
  ChkMat.Is_kappa_0 = false;
  ChkMat.Is_friction_angle = false;
  ChkMat.Is_dilatancy_angle = false;
  ChkMat.Is_J2_degradated = false;
  ChkMat.Is_Ceps = false;
  ChkMat.Is_Gf = false;
  ChkMat.Is_ft = false;
  ChkMat.Is_heps = false;
  ChkMat.Is_wcrit = false;

  return ChkMat;
}

/**********************************************************************/

static int __read_boolean(bool *Value, char *Value_char) {
  if ((strcmp(Value_char, "True") == 0) ||
      (strcmp(Value_char, "TRUE") == 0 || (strcmp(Value_char, "true") == 0)) ||
      (strcmp(Value_char, "1") == 0)) {
    *Value = true;
    return EXIT_SUCCESS;
  } else if ((strcmp(Value_char, "False") == 0) ||
             (strcmp(Value_char, "FALSE") == 0 ||
              (strcmp(Value_char, "false") == 0) ||
              (strcmp(Value_char, "0") == 0))) {
    *Value = false;
    return EXIT_SUCCESS;
  } else {
    fprintf(stderr, RED "Undefined %s \n" RESET, Value_char);
    return EXIT_FAILURE;
  }
}

/**********************************************************************/

static int __check_material(Material *DP_Material, Check_Material ChkMat,
                            int Idx) {

  int STATUS = EXIT_SUCCESS;

  printf("\t -> %s \n", "Drucker-Prager material");

  if (ChkMat.Is_rho == true) {
    printf("\t \t -> %s : %f \n", "" MAGENTA "[rho]" RESET "",
           (*DP_Material).rho);
  } else {
    return EXIT_FAILURE;
  }

  if (ChkMat.Is_E == true) {
    printf("\t \t -> %s : %f \n", "" MAGENTA "[E]" RESET "", (*DP_Material).E);
  } else {
    return EXIT_FAILURE;
  }

  if (ChkMat.Is_nu == true) {
    printf("\t \t -> %s : %f \n", "" MAGENTA "[nu]" RESET "",
           (*DP_Material).nu);
  } else {
    return EXIT_FAILURE;
  }

  if (ChkMat.Is_friction_angle == true) {
    printf("\t \t -> %s : %f \n", "" MAGENTA "[Friction-angle]" RESET "",
           (*DP_Material).phi_Frictional);
  } else {
    return EXIT_FAILURE;
  }

  if (ChkMat.Is_dilatancy_angle) {
    printf("\t \t -> %s : %f \n", "" MAGENTA "[Dilatancy-angle]" RESET "",
           (*DP_Material).psi_Frictional);
  } else {
    return EXIT_FAILURE;
  }

  printf("\t \t -> %s : %f \n", "" MAGENTA "[Kappa-0]" RESET "",
         (*DP_Material).kappa_0);

  printf("\t \t -> %s : %f \n",
         "" MAGENTA "[Reference-plastic-strain]" RESET "",
         (*DP_Material).Plastic_Strain_0);

  printf("\t \t -> %s : %f \n", "" MAGENTA "[Hardening-modulus]" RESET "",
         (*DP_Material).Hardening_modulus);

  printf("\t \t -> %s : %f \n", "" MAGENTA "[m]" RESET "",
         (*DP_Material).Exponent_Hardening_Ortiz);

  printf("\t \t -> %s : %f \n", "" MAGENTA "[J2-degradated]" RESET "",
         (*DP_Material).J2_degradated);

  printf("\t \t -> %s : %f \n", "" MAGENTA "[Reference-pressure]" RESET "",
         (*DP_Material).ReferencePressure);

  if ((Driver_EigenErosion == true) || (Driver_EigenSoftening == true)) {
    if (ChkMat.Is_Ceps == true) {
      printf("\t \t -> %s : %f \n", "" MAGENTA "[Ceps]" RESET "",
             (*DP_Material).Ceps);
    } else {
      return EXIT_FAILURE;
    }
  }

  if (Driver_EigenErosion == true) {
    if (ChkMat.Is_Gf == true) {
      printf("\t \t -> %s : %f \n", "" MAGENTA "[Gf]" RESET "",
             (*DP_Material).Gf);
    } else {
      return EXIT_FAILURE;
    }
  }

  if (Driver_EigenSoftening == true) {
    if (ChkMat.Is_ft == true) {
      printf("\t \t -> %s : %f \n", "" MAGENTA "[ft]" RESET "",
             (*DP_Material).ft);
    } else {
      return EXIT_FAILURE;
    }

    if (ChkMat.Is_heps == true) {
      printf("\t \t -> %s : %f \n", "" MAGENTA "[heps]" RESET "",
             (*DP_Material).heps);
    } else {
      return EXIT_FAILURE;
    }

    if (ChkMat.Is_wcrit == true) {
      printf("\t \t -> %s : %f \n", "" MAGENTA "[wcrit]" RESET "",
             (*DP_Material).wcrit);
    } else {
      return EXIT_FAILURE;
    }
  }

  return STATUS;
}

/**********************************************************************/

static void __print_miss_variables(Check_Material ChkMat) {

  if (ChkMat.Is_rho == false) {
    fputs("\t \t -> " MAGENTA "[rho]" RESET " : " RED "false" RESET " \n",
          stderr);
  }

  if (ChkMat.Is_E == false) {
    fputs("\t \t -> " MAGENTA "[E]" RESET " : " RED "false" RESET " \n",
          stderr);
  }

  if (ChkMat.Is_nu == false) {
    fputs("\t \t -> " MAGENTA "[nu]" RESET " : " RED "false" RESET " \n",
          stderr);
  }

  if (ChkMat.Is_friction_angle == false) {
    fputs("\t \t -> " MAGENTA "[Friction-angle]" RESET " : " RED "false" RESET
          " \n",
          stderr);
  }

  if (ChkMat.Is_dilatancy_angle == false) {
    fputs("\t \t -> " MAGENTA "[Dilatancy-angle]" RESET " : " RED "false" RESET
          " \n",
          stderr);
  }

  if (ChkMat.Is_kappa_0 == false) {
    fputs("\t \t -> " MAGENTA "[kappa-0]" RESET " : " CYAN "optional" RESET
          " \n",
          stderr);
  }

  if (ChkMat.Is_Hardening_modulus == false) {
    fputs("\t \t -> " MAGENTA "[Hardening-modulus]" RESET " : " CYAN
          "optional" RESET " \n",
          stderr);
  }

  if (ChkMat.Is_Hardening_Exponent == false) {
    fputs("\t \t -> " MAGENTA "[m]" RESET " : " CYAN "optional" RESET " \n",
          stderr);
  }

  if (ChkMat.Is_J2_degradated == false) {
    fputs("\t \t -> " MAGENTA "[J2-degradated]" RESET " : " CYAN
          "optional" RESET " \n",
          stderr);
  }

  if (ChkMat.Is_Reference_pressure == false) {
    fputs("\t \t -> " MAGENTA "[Reference-pressure]" RESET " : " CYAN
          "optional" RESET " \n",
          stderr);
  }

  if ((Driver_EigenErosion == true) || (Driver_EigenSoftening == true)) {

    if (ChkMat.Is_Ceps == false) {
      fputs("\t \t -> " MAGENTA "[Ceps]" RESET " : " RED "false" RESET " \n",
            stderr);
    }
  }

  if (Driver_EigenErosion == true) {
    if (ChkMat.Is_Gf == false) {
      fputs("\t \t -> " MAGENTA "[Gf]" RESET " : " RED "false" RESET " \n",
            stderr);
    }
  }

  if (Driver_EigenSoftening == true) {

    if (ChkMat.Is_ft == false) {
      fputs("\t \t -> " MAGENTA "[ft]" RESET " : " RED "false" RESET " \n",
            stderr);
    }

    if (ChkMat.Is_heps == false) {
      fputs("\t \t -> " MAGENTA "[heps]" RESET " : " RED "false" RESET " \n",
            stderr);
    }

    if (ChkMat.Is_wcrit == false) {
      fputs("\t \t -> " MAGENTA "[wcrit]" RESET " : " RED "false" RESET " \n",
            stderr);
    }
  }
}

/**********************************************************************/
