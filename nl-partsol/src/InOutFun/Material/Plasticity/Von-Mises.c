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
  bool Is_rho;            // Reference fensity
  bool Is_E;              // Young modulus
  bool Is_nu;             // Poisson cefficient
  bool Is_Plastic_solver; // FE or BE solver
  bool Is_yield_stress;   // Initial Yield stress

  bool Is_Hardening_modulus;
  bool Is_Hardening;
  bool Is_Hardening_Hughes;
  bool Is_Parameter_Hardening_Hughes;
  bool Is_Hardening_Cervera;
  bool Is_Hardening_Ortiz;
  bool Is_Exponent_Hardening_Ortiz;
  bool Is_Reference_Plastic_Strain_Ortiz;

  bool Is_Viscous_regularization;
  bool Is_fluidity_param;

  bool Is_Locking_Control_Fbar; // Locking control
  bool Is_alpha_Fbar;           // Tunning paramer for the F-bar

} Check_Material;

static bool Activate_Options(char *);
static void standard_error();
static Check_Material Initialise_Check_Material();
static void check_Von_Mises_Material(Material, Check_Material, int);

/***************************************************************************/

Material Define_Von_Mises(FILE *Simulation_file, char *Material_Model,
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

  /* Default options */
  New_Material.Hardening_Ortiz = false;
  New_Material.Eigenerosion = false;
  New_Material.Eigensoftening = false;
  New_Material.Hardening_Hughes = false;
  New_Material.Hardening_Cervera = false;
  New_Material.Exponent_Hardening_Ortiz = false;
  New_Material.Viscous_regularization = false;
  New_Material.Locking_Control_Fbar = false;
  New_Material.alpha_Fbar = 0.0;

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
    else if (strcmp(Parameter_pars[0], "Plastic-Solver") == 0) {
      ChkMat.Is_Plastic_solver = true;
      strcpy(New_Material.Plastic_Solver, Parameter_pars[1]);
    }
    /**************************************************/
    else if (strcmp(Parameter_pars[0], "Yield-stress") == 0) {
      ChkMat.Is_yield_stress = true;
      New_Material.yield_stress_0 = atof(Parameter_pars[1]);
    }
    /**************************************************/
    else if (strcmp(Parameter_pars[0], "Hardening-Criteria") == 0) {
      ChkMat.Is_Hardening = true;

      if (strcmp(Parameter_pars[1], "Hughes") == 0) {
        ChkMat.Is_Hardening_Hughes = true;
        New_Material.Hardening_Hughes = true;
      } else if (strcmp(Parameter_pars[1], "Cervera") == 0) {
        ChkMat.Is_Hardening_Cervera = true;
        New_Material.Hardening_Cervera = true;
      } else if (strcmp(Parameter_pars[1], "Ortiz") == 0) {
        ChkMat.Is_Hardening_Ortiz = true;
        New_Material.Hardening_Ortiz = true;
      } else {
        sprintf(Error_message, "%s",
                "Options for Hardening-Criteria -> Hughes/Cervera/Ortiz");
        standard_error(Error_message);
      }
    }
    /**************************************************/
    else if (strcmp(Parameter_pars[0], "Hardening-Modulus") == 0) {
      ChkMat.Is_Hardening_modulus = true;
      New_Material.Hardening_modulus = atof(Parameter_pars[1]);
    }
    /**************************************************/

    else if (strcmp(Parameter_pars[0], "Parameter-Hardening-Hughes") == 0) {
      ChkMat.Is_Parameter_Hardening_Hughes = true;
      New_Material.Parameter_Hardening_Hughes = atof(Parameter_pars[1]);
    }
    /**************************************************/
    else if (strcmp(Parameter_pars[0], "Exponent-Hardening-Ortiz") == 0) {
      ChkMat.Is_Exponent_Hardening_Ortiz = true;
      New_Material.Exponent_Hardening_Ortiz = atof(Parameter_pars[1]);
    }
    /**************************************************/
    else if (strcmp(Parameter_pars[0], "Reference-Plastic-Strain-Ortiz") == 0) {
      ChkMat.Is_Reference_Plastic_Strain_Ortiz = true;
      New_Material.Reference_Plastic_Strain_Ortiz = atof(Parameter_pars[1]);
    }
    /**************************************************/
    else if (strcmp(Parameter_pars[0], "Viscous-regularization") == 0) {
      ChkMat.Is_Viscous_regularization = true;

      if (strcmp(Parameter_pars[1], "true") == 0) {
        New_Material.Viscous_regularization = true;
      } else if (strcmp(Parameter_pars[1], "false") == 0) {
        New_Material.Viscous_regularization = false;
      } else {
        sprintf(Error_message, "%s",
                "Options for Viscous-regularization -> true/false");
        standard_error(Error_message);
      }
    }
    /**************************************************/
    else if (strcmp(Parameter_pars[0], "Fluidity-Parameter") == 0) {
      ChkMat.Is_fluidity_param = true;
      New_Material.fluidity_param = atof(Parameter_pars[1]);
    }
    /**************************************************/
    else if (strcmp(Parameter_pars[0], "Hardening-Ortiz") == 0) {
      New_Material.Hardening_Ortiz = Activate_Options(Parameter_pars[1]);
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

  check_Von_Mises_Material(New_Material, ChkMat, Material_Idx);

  TOL_Radial_Returning = 1E-10;
  Max_Iterations_Radial_Returning = 300;

  /* Return outputs */
  return New_Material;
}
/***************************************************************************/

static Check_Material Initialise_Check_Material() {
  Check_Material ChkMat;

  ChkMat.Is_rho = false;            // Reference fensity
  ChkMat.Is_E = false;              // Young modulus
  ChkMat.Is_nu = false;             // Poisson cefficient
  ChkMat.Is_Plastic_solver = false; // FE or BE solver
  ChkMat.Is_yield_stress = false;   // Initial Yield stress
  ChkMat.Is_Hardening_modulus = false;
  ChkMat.Is_Hardening = false;
  ChkMat.Is_Hardening_Hughes = false;
  ChkMat.Is_Parameter_Hardening_Hughes = false;
  ChkMat.Is_Hardening_Cervera = false;
  ChkMat.Is_Hardening_Ortiz = false;
  ChkMat.Is_Exponent_Hardening_Ortiz = false;
  ChkMat.Is_Reference_Plastic_Strain_Ortiz = false;
  ChkMat.Is_Viscous_regularization = false; // Viscous regularization
  ChkMat.Is_fluidity_param = false;         // Viscoplasticity parameter
  ChkMat.Is_Locking_Control_Fbar = false;   // Locking control
  ChkMat.Is_alpha_Fbar = false;             // Tunning parameter for the F-bar

  return ChkMat;
}

/***************************************************************************/

static void check_Von_Mises_Material(Material Mat_particle,
                                     Check_Material ChkMat, int Idx) {
  if (ChkMat.Is_rho && ChkMat.Is_E && ChkMat.Is_nu && ChkMat.Is_yield_stress &&
      ChkMat.Is_Plastic_solver) {
    printf("\t -> %s \n", "Von-Mises material");
    printf("\t \t -> %s : %f \n", "Density", Mat_particle.rho);
    printf("\t \t -> %s : %f \n", "Elastic modulus", Mat_particle.E);
    printf("\t \t -> %s : %f \n", "Poisson modulus", Mat_particle.nu);
    printf("\t \t -> %s : %f \n", "Yield stress", Mat_particle.yield_stress_0);
    printf("\t \t -> %s : %s \n", "Plastic solver",
           Mat_particle.Plastic_Solver);

    if (Mat_particle.Locking_Control_Fbar) {
      printf("\t \t -> %s : %s \n", "F-bar", "Enabled");

      if (ChkMat.Is_alpha_Fbar) {
        printf("\t \t -> %s : %f \n", "alpha F-bar", Mat_particle.alpha_Fbar);
      }
    } else {
      printf("\t \t -> %s : %s \n", "F-bar", "Disabled");
    }

    if (ChkMat.Is_Hardening) {
      if (ChkMat.Is_Hardening_Hughes) {
        if (ChkMat.Is_Hardening_modulus &&
            ChkMat.Is_Parameter_Hardening_Hughes) {
          printf("\t \t -> %s : %f \n", "Hardening-Modulus",
                 Mat_particle.Hardening_modulus);
          printf("\t \t -> %s : %f \n", "Parameter-Hardening-Hughes",
                 Mat_particle.Parameter_Hardening_Hughes);
        } else {
          fprintf(stderr, "%s : %s \n", "Error in GramsMaterials()",
                  "Some parameter is missed for Von-Mises material (Hughes "
                  "Hardening)");
          fputs(ChkMat.Is_Hardening_modulus ? "Hardening-Modulus : true \n"
                                            : "Hardening-Modulus : false \n",
                stdout);
          fputs(ChkMat.Is_Parameter_Hardening_Hughes
                    ? "Parameter-Hardening-Hughes : true \n"
                    : "Parameter-Hardening-Hughes : false \n",
                stdout);
          exit(EXIT_FAILURE);
        }

      } else if (ChkMat.Is_Hardening_Cervera) {

        if (strcmp(Mat_particle.Plastic_Solver, "Forward-Euler") == 0) {
          fprintf(stderr, "%s : %s \n", "Error in GramsMaterials()",
                  "Switch to Backward-Euler for non-linear Hardening laws)");
          exit(EXIT_FAILURE);
        }

        if (ChkMat.Is_Hardening_modulus) {
          printf("\t \t -> %s : %f \n", "Hardening-Modulus",
                 Mat_particle.Hardening_modulus);
        } else {
          fprintf(stderr, "%s : %s \n", "Error in GramsMaterials()",
                  "Some parameter is missed for Von-Mises material (Cervera "
                  "Hardening)");
          fputs(ChkMat.Is_Hardening_modulus ? "Hardening-Modulus : true \n"
                                            : "Hardening-Modulus : false \n",
                stdout);
          exit(EXIT_FAILURE);
        }
      } else if (ChkMat.Is_Hardening_Ortiz) {

        if (strcmp(Mat_particle.Plastic_Solver, "Forward-Euler") == 0) {
          fprintf(stderr, "%s : %s \n", "Error in GramsMaterials()",
                  "Switch to Backward-Euler for non-linear Hardening laws)");
          exit(EXIT_FAILURE);
        }

        if (ChkMat.Is_Exponent_Hardening_Ortiz &&
            ChkMat.Is_Reference_Plastic_Strain_Ortiz) {
          printf("\t \t -> %s : %f \n", "Exponent-Hardening-Ortiz",
                 Mat_particle.Exponent_Hardening_Ortiz);
          printf("\t \t -> %s : %f \n", "Reference-Plastic-Strain_Ortiz",
                 Mat_particle.Reference_Plastic_Strain_Ortiz);
        } else {
          fprintf(stderr, "%s : %s \n", "Error in GramsMaterials()",
                  "Some parameter is missed for Von-Mises material (Ortiz "
                  "Hardening)");
          fputs(ChkMat.Is_Exponent_Hardening_Ortiz
                    ? "Exponent-Hardening-Ortiz : true \n"
                    : "Exponent-Hardening-Ortiz : false \n",
                stdout);
          fputs(ChkMat.Is_Reference_Plastic_Strain_Ortiz
                    ? "Reference-Plastic-Strain_Ortiz : true \n"
                    : "Reference-Plastic-Strain_Ortiz : false \n",
                stdout);
          exit(EXIT_FAILURE);
        }
      }
    }

    if (ChkMat.Is_Viscous_regularization) {
      if (ChkMat.Is_fluidity_param) {
        printf("\t \t -> %s : %f \n", "Fluidity parameter",
               Mat_particle.fluidity_param);
      } else {
        fprintf(stderr, "%s : %s \n", "Error in GramsMaterials()",
                "Some parameter is missed for Von-Mises material (Viscoplastic "
                "regularization)");
        fputs(ChkMat.Is_fluidity_param ? "Fluidity parameter : true \n"
                                       : "Fluidity parameter : false \n",
              stdout);
        exit(EXIT_FAILURE);
      }
    }

  } else {
    fprintf(stderr, "%s : %s \n", "Error in GramsMaterials()",
            "Some parameter is missed for Von-Mises material");
    fputs(ChkMat.Is_rho ? "Density : true \n" : "Density : false \n", stdout);
    fputs(ChkMat.Is_E ? "Elastic modulus : true \n"
                      : "Elastic modulus : false \n",
          stdout);
    fputs(ChkMat.Is_nu ? "Poisson modulus : true \n"
                       : "Poisson modulus : false \n",
          stdout);
    fputs(ChkMat.Is_yield_stress ? "Yield stress : true \n"
                                 : "Yield stress : false \n",
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