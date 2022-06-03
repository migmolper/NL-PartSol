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

  bool Is_rho;               // Reference fensity
  bool Is_E;                 // Young modulus
  bool Is_nu;                // Poisson cefficient
  bool Is_Yield;             // Initial Yield stress
  bool Is_Hardening_modulus; // Hardening modulus
  bool Is_theta;
  bool Is_K_0;
  bool Is_K_inf;
  bool Is_delta;
  bool Is_Ceps;  // Normalizing constant (Eigenerosion)
  bool Is_Gf;    // Failure energy (Eigenerosion)
  bool Is_ft;    // Tensile strengt of the material
  bool Is_heps;  // Bandwidth of the cohesive fracture
  bool Is_wcrit; // Critical opening displacement

} Check_Material;

static Check_Material __Initialise_Check_Material();
static int __check_material(Material *, Check_Material, int);

/***************************************************************************/

int Define_Von_Mises(Material *VM_Material, FILE *Simulation_file,
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

  /* Default options */
  (*VM_Material).Hardening_modulus = 0.0;
  (*VM_Material).theta_Hardening_Voce = 1.0;
  (*VM_Material).K_0_Hardening_Voce = 0.0;
  (*VM_Material).K_inf_Hardening_Voce = 0.0;
  (*VM_Material).delta_Hardening_Voce = 0.0;
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
      (*VM_Material).rho = atof(Parameter_pars[1]);

      if ((*VM_Material).rho <= 0.0) {
        fprintf(stderr, "" RED " Invalid value of [rho] (0,+inf)" RESET " \n");
        return EXIT_FAILURE;
      }

    }
    /**************************************************/
    else if (strcmp(Parameter_pars[0], "E") == 0) {
      ChkMat.Is_E = true;
      (*VM_Material).E = atof(Parameter_pars[1]);

      if ((*VM_Material).E <= 0.0) {
        fprintf(stderr, "" RED " Invalid value of [E] (0,+inf) " RESET " \n");
        return EXIT_FAILURE;
      }

    }
    /**************************************************/
    else if (strcmp(Parameter_pars[0], "nu") == 0) {
      ChkMat.Is_nu = true;
      (*VM_Material).nu = atof(Parameter_pars[1]);

      if ((*VM_Material).nu < 0.0) {
        fprintf(stderr, "" RED " Invalid value of [nu] [0,+inf)" RESET " \n");
        return EXIT_FAILURE;
      }

    }
    /**************************************************/
    else if (strcmp(Parameter_pars[0], "Yield-stress") == 0) {
      ChkMat.Is_Yield = true;
      (*VM_Material).kappa_0 = atof(Parameter_pars[1]);

      if ((*VM_Material).kappa_0 <= 0.0) {
        fprintf(stderr,
                "" RED " Invalid value of [Yield-stress] (0,+inf) " RESET
                " \n");
        return EXIT_FAILURE;
      }

    }
    /**************************************************/
    else if (strcmp(Parameter_pars[0], "Hardening-Modulus") == 0) {
      ChkMat.Is_Hardening_modulus = true;
      (*VM_Material).Hardening_modulus = atof(Parameter_pars[1]);

      if ((*VM_Material).Hardening_modulus <= 0.0) {
        fprintf(stderr,
                "" RED " Invalid value of [Hardening-Modulus] (0,+inf)] " RESET
                " \n");
        return EXIT_FAILURE;
      }

    }
    /**************************************************/
    else if (strcmp(Parameter_pars[0], "theta") == 0) {
      ChkMat.Is_theta = true;
      (*VM_Material).theta_Hardening_Voce = atof(Parameter_pars[1]);

      if (((*VM_Material).theta_Hardening_Voce < 0.0) ||
          ((*VM_Material).theta_Hardening_Voce > 1.0)) {
        fprintf(stderr, "" RED " Invalid value of [theta] [0,1]" RESET " \n");
        return EXIT_FAILURE;
      }

    }
    /**************************************************/
    else if (strcmp(Parameter_pars[0], "K-0") == 0) {
      ChkMat.Is_K_0 = true;
      (*VM_Material).K_0_Hardening_Voce = atof(Parameter_pars[1]);

      if ((*VM_Material).K_0_Hardening_Voce <= 0.0) {
        fprintf(stderr, "" RED " Invalid value of [K-0] (0,+inf)" RESET " \n");
        return EXIT_FAILURE;
      }

    }
    /**************************************************/
    else if (strcmp(Parameter_pars[0], "K-inf") == 0) {
      ChkMat.Is_K_inf = true;
      (*VM_Material).K_inf_Hardening_Voce = atof(Parameter_pars[1]);

      if ((*VM_Material).K_inf_Hardening_Voce <= 0.0) {
        fprintf(stderr,
                "" RED " Invalid value of [K-inf] (0,+inf)" RESET " \n");
        return EXIT_FAILURE;
      }

    }
    /**************************************************/
    else if (strcmp(Parameter_pars[0], "delta") == 0) {
      ChkMat.Is_delta = true;
      (*VM_Material).delta_Hardening_Voce = atof(Parameter_pars[1]);

      if ((*VM_Material).delta_Hardening_Voce < 0.0) {
        fprintf(stderr,
                "" RED " Invalid value of [delta] [0,+inf)" RESET " \n");
        return EXIT_FAILURE;
      }

    }
    /**************************************************/
    else if (strcmp(Parameter_pars[0], "Ceps") == 0) {
      if ((Driver_EigenErosion == true) || (Driver_EigenSoftening == true)) {
        ChkMat.Is_Ceps = true;
        (*VM_Material).Ceps = atof(Parameter_pars[1]);
      }
    }
    /**************************************************/
    else if (strcmp(Parameter_pars[0], "Gf") == 0) {
      if (Driver_EigenErosion == true) {
        ChkMat.Is_Gf = true;
        (*VM_Material).Gf = atof(Parameter_pars[1]);
      }
    }
    /**************************************************/
    else if (strcmp(Parameter_pars[0], "ft") == 0) {
      if (Driver_EigenSoftening == true) {
        ChkMat.Is_ft = true;
        (*VM_Material).ft = atof(Parameter_pars[1]);
      }
    }
    /**************************************************/
    else if (strcmp(Parameter_pars[0], "heps") == 0) {
      if (Driver_EigenSoftening == true) {
        ChkMat.Is_heps = true;
        (*VM_Material).heps = atof(Parameter_pars[1]);
      }
    }
    /**************************************************/
    else if (strcmp(Parameter_pars[0], "wcrit") == 0) {
      if (Driver_EigenSoftening == true) {
        ChkMat.Is_wcrit = true;
        (*VM_Material).wcrit = atof(Parameter_pars[1]);
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

  if ((*VM_Material).K_inf_Hardening_Voce < (*VM_Material).K_0_Hardening_Voce) {
    fprintf(stderr,
            "" RED " [K-inf] should be larger or equal than [K-0] " RESET
            " \n");
    return EXIT_FAILURE;
  }

  strcpy((*VM_Material).Type, Material_Model);

  __check_material(VM_Material, ChkMat, Material_Idx);
  if (STATUS == EXIT_FAILURE) {
    fprintf(stderr, "" RED " Error in __check_material() " RESET " \n");
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
  ChkMat.Is_Hardening_modulus = false;
  ChkMat.Is_Yield = false;
  ChkMat.Is_theta = false;
  ChkMat.Is_K_0 = false;
  ChkMat.Is_K_inf = false;
  ChkMat.Is_delta = false;
  ChkMat.Is_Ceps = false;
  ChkMat.Is_Gf = false;
  ChkMat.Is_ft = false;
  ChkMat.Is_heps = false;
  ChkMat.Is_wcrit = false;

  return ChkMat;
}

/***************************************************************************/

static int __check_material(Material *VM_Material, Check_Material ChkMat,
                            int Idx) {

  int STATUS = EXIT_SUCCESS;

  if (ChkMat.Is_rho && ChkMat.Is_E && ChkMat.Is_nu && ChkMat.Is_Yield &&
      ChkMat.Is_Yield) {

    printf("\t -> %s \n", "Von-Mises material");

    printf("\t \t -> %s : %f \n", "" MAGENTA "[rho]" RESET "",
           (*VM_Material).rho);

    printf("\t \t -> %s : %f \n", "" MAGENTA "[E]" RESET "", (*VM_Material).E);

    printf("\t \t -> %s : %f \n", "" MAGENTA "[nu]" RESET "",
           (*VM_Material).nu);

    printf("\t \t -> %s : %f \n", "" MAGENTA "[Yield-stress]" RESET "",
           (*VM_Material).kappa_0);

    printf("\t \t -> %s : %f \n", "" MAGENTA "[Hardening-Modulus]" RESET "",
           (*VM_Material).Hardening_modulus);

    printf("\t \t -> %s : %f \n", "" MAGENTA "[theta]" RESET "",
           (*VM_Material).theta_Hardening_Voce);

    printf("\t \t -> %s : %f \n", "" MAGENTA "[K-0]" RESET "",
           (*VM_Material).K_0_Hardening_Voce);

    printf("\t \t -> %s : %f \n", "" MAGENTA "[K-inf]" RESET "",
           (*VM_Material).K_inf_Hardening_Voce);

    printf("\t \t -> %s : %f \n", "" MAGENTA "[delta]" RESET "",
           (*VM_Material).delta_Hardening_Voce);

    if ((Driver_EigenErosion == true) || (Driver_EigenSoftening == true)) {
      printf("\t \t -> %s : %f \n", "" MAGENTA "[Ceps]" RESET "",
             (*VM_Material).Ceps);
    }

    if (Driver_EigenErosion == true) {
      printf("\t \t -> %s : %f \n", "" MAGENTA "[Gf]" RESET "",
             (*VM_Material).Gf);
    }

    if (Driver_EigenSoftening == true) {
      printf("\t \t -> %s : %f \n", "" MAGENTA "[ft]" RESET "",
             (*VM_Material).ft);

      printf("\t \t -> %s : %f \n", "" MAGENTA "[heps]" RESET "",
             (*VM_Material).heps);

      printf("\t \t -> %s : %f \n", "" MAGENTA "[wcrit]" RESET "",
             (*VM_Material).wcrit);
    }

  } else {

    fprintf(stderr, "" RED " %s : %s " RESET " \n", "Error in GramsMaterials()",
            "Some parameter is missed for Von-Mises");

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

    fputs(ChkMat.Is_Yield
              ? "" MAGENTA "[Yield-stress]" RESET " : " GREEN "true" RESET " \n"
              : "" MAGENTA "[Yield-stress]" RESET " : " RED "false" RESET " \n",
          stderr);

    fputs(ChkMat.Is_Hardening_modulus ? "" MAGENTA "[Hardening-Modulus]" RESET
                                        " : " GREEN "true" RESET " \n"
                                      : "" MAGENTA "[Hardening-Modulus]" RESET
                                        " : " RED "false" RESET " \n",
          stderr);

    fputs(ChkMat.Is_theta
              ? "" MAGENTA "[theta]" RESET " : " GREEN "true" RESET " \n"
              : "" MAGENTA "[theta]" RESET " : " RED "false" RESET " \n",
          stderr);

    fputs(ChkMat.Is_K_0
              ? "" MAGENTA "[K-0]" RESET " : " GREEN "true" RESET " \n"
              : "" MAGENTA "[K-0]" RESET " : " RED "false" RESET " \n",
          stderr);

    fputs(ChkMat.Is_K_inf
              ? "" MAGENTA "[K-inf]" RESET " : " GREEN "true" RESET " \n"
              : "" MAGENTA "[K-inf]" RESET " : " RED "false" RESET " \n",
          stderr);

    fputs(ChkMat.Is_delta
              ? "" MAGENTA "[delta]" RESET " : " GREEN "true" RESET " \n"
              : "" MAGENTA "[delta]" RESET " : " RED "false" RESET " \n",
          stderr);

    if ((Driver_EigenErosion == true) || (Driver_EigenSoftening == true)) {
      fputs(ChkMat.Is_Ceps
                ? "" MAGENTA "[Ceps]" RESET " : " GREEN "true" RESET " \n"
                : "" MAGENTA "[Ceps]" RESET " : " RED "false" RESET " \n",
            stderr);
    }

    if (Driver_EigenErosion == true) {
      fputs(ChkMat.Is_Gf
                ? "" MAGENTA "[Gf]" RESET " : " GREEN "true" RESET " \n"
                : "" MAGENTA "[Gf]" RESET " : " RED "false" RESET " \n",
            stderr);
    }

    if (Driver_EigenSoftening == true) {
      fputs(ChkMat.Is_ft
                ? "" MAGENTA "[ft]" RESET " : " GREEN "true" RESET " \n"
                : "" MAGENTA "[ft]" RESET " : " RED "false" RESET " \n",
            stderr);

      fputs(ChkMat.Is_heps
                ? "" MAGENTA "[heps]" RESET " : " GREEN "true" RESET " \n"
                : "" MAGENTA "[heps]" RESET " : " RED "false" RESET " \n",
            stderr);

      fputs(ChkMat.Is_wcrit
                ? "" MAGENTA "[wcrit]" RESET " : " GREEN "true" RESET " \n"
                : "" MAGENTA "[wcrit]" RESET " : " RED "false" RESET " \n",
            stderr);
    }

    STATUS = EXIT_FAILURE;
  }

  return STATUS;
}

/***************************************************************************/
