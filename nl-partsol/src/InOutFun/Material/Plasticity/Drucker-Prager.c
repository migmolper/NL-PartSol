#include <string.h>
#include <math.h>
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

typedef struct {

  bool Is_rho;                       // Reference fensity
  bool Is_E;                         // Young modulus
  bool Is_nu;                        // Poisson cefficient
  bool Is_Hardening_Exponent;        // Hardening exponent
  bool Is_Hardening_modulus;         // Hardening modulus
  bool Is_Reference_pressure;      // Reference pressure
  bool Is_Reference_plastic_strain; // Reference Plastic Strain
  bool Is_kappa_0;                   // Initial yield
  bool Is_friction_angle;            // Friction angle
  bool Is_dilatancy_angle;           // Dilatancy angle
  bool Is_J2_degradated; // Degradation limit

} Check_Material;

static Check_Material __Initialise_Check_Material();
static int __check_DP_Material(Material *, Check_Material, int);

/**********************************************************************/

int Define_Drucker_Prager(Material *DP_Material,FILE *Simulation_file, char *Material_Model, int Material_Idx) {

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
  (*DP_Material).yield_stress_0 = 0.0;
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
      (*DP_Material).Reference_Plastic_Strain_Ortiz = atof(Parameter_pars[1]);
    }
    /**************************************************/
    else if (strcmp(Parameter_pars[0], "kappa-0") == 0) {
      ChkMat.Is_kappa_0 = true;
      (*DP_Material).yield_stress_0 = atof(Parameter_pars[1]);
    }
    /**************************************************/
    else if (strcmp(Parameter_pars[0], "Friction-angle") == 0) {
      ChkMat.Is_friction_angle = true;
      (*DP_Material).phi_Frictional = atof(Parameter_pars[1]);

      if(((*DP_Material).phi_Frictional > 90.0) && ((*DP_Material).phi_Frictional < 0.0))
      {
        fprintf(stderr, ""RED" Invalid value of the [Friction-angle] "RESET" \n");
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

      if(((*DP_Material).psi_Frictional > 90.0) && ((*DP_Material).psi_Frictional < 0.0))
      {
        fprintf(stderr, ""RED" Invalid value of the [Dilatancy-angle] "RESET" \n");
        return EXIT_FAILURE;
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

  double kappa0 = (*DP_Material).yield_stress_0;
  double m = (*DP_Material).Exponent_Hardening_Ortiz;
  double H = (*DP_Material).Hardening_modulus;

  if(((m>0) && (H>0)) || ((m<0) && (H<0)))
  {
    if(ChkMat.Is_Reference_plastic_strain == false)
    {
      (*DP_Material).Reference_Plastic_Strain_Ortiz = (kappa0 / (m * H)) *
      pow(1, (1.0 / m - 1.0));
    }
  }
  else
  {
    fprintf(stderr, ""RED" [Hardening-modulus] and [m] must have the same sign "RESET" \n");
    return EXIT_FAILURE;
  }

  STATUS = __check_DP_Material(DP_Material, ChkMat, Material_Idx);
  if(STATUS == EXIT_FAILURE)
  {
        fprintf(stderr, ""RED" Error in __check_DP_Material() "RESET" \n");
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

  return ChkMat;
}

/**********************************************************************/

static int __check_DP_Material(Material * DP_Material,
                                      Check_Material ChkMat, int Idx) {

int STATUS = EXIT_SUCCESS;

  if (ChkMat.Is_rho && ChkMat.Is_E && ChkMat.Is_nu && ChkMat.Is_Reference_pressure &&
      ChkMat.Is_friction_angle && ChkMat.Is_dilatancy_angle) {


    printf("\t -> %s \n", "Drucker-Prager material");

    printf("\t \t -> %s : %f \n", ""MAGENTA"[rho]"RESET"", (*DP_Material).rho);

    printf("\t \t -> %s : %f \n", ""MAGENTA"[E]"RESET"", (*DP_Material).E);

    printf("\t \t -> %s : %f \n", ""MAGENTA"[nu]"RESET"", (*DP_Material).nu);

    printf("\t \t -> %s : %f \n", ""MAGENTA"[Friction-angle]"RESET"",
           (*DP_Material).phi_Frictional);

    printf("\t \t -> %s : %f \n", ""MAGENTA"[Dilatancy-angle]"RESET"",
           (*DP_Material).psi_Frictional);

    printf("\t \t -> %s : %f \n", ""MAGENTA"[Kappa-0]"RESET"", (*DP_Material).yield_stress_0);

    printf("\t \t -> %s : %f \n", ""MAGENTA"[Reference-plastic-strain]"RESET"", (*DP_Material).Reference_Plastic_Strain_Ortiz);

    printf("\t \t -> %s : %f \n", ""MAGENTA"[Hardening-modulus]"RESET"", (*DP_Material).Hardening_modulus);

    printf("\t \t -> %s : %f \n", ""MAGENTA"[m]"RESET"", (*DP_Material).Exponent_Hardening_Ortiz);

    printf("\t \t -> %s : %f \n", ""MAGENTA"[J2-degradated]"RESET"", (*DP_Material).J2_degradated);

    printf("\t \t -> %s : %f \n", ""MAGENTA"[Reference-pressure]"RESET"", (*DP_Material).ReferencePressure);

  } else {
    
    
    fprintf(stderr,""RED" %s : %s "RESET" \n", "Error in GramsMaterials()",
            "Some parameter is missed for Drucker-Prager");
    
    fputs(ChkMat.Is_rho ? 
    ""MAGENTA"[rho]"RESET" : "GREEN"true"RESET" \n" : 
    ""MAGENTA"[rho]"RESET" : "RED"false"RESET" \n", stderr);
    
    fputs(ChkMat.Is_E ? 
    ""MAGENTA"[E]"RESET" : "GREEN"true"RESET" \n" :
    ""MAGENTA"[E]"RESET" : "RED"false"RESET" \n",stderr);

    fputs(ChkMat.Is_nu ? 
    ""MAGENTA"[nu]"RESET" : "GREEN"true"RESET" \n"  : 
    ""MAGENTA"[nu]"RESET" : "RED"false"RESET" \n",stderr);

    fputs(ChkMat.Is_Hardening_Exponent ? 
    ""MAGENTA"[m]"RESET" : "GREEN"true"RESET" \n" : 
    ""MAGENTA"[m]"RESET" : "RED"false"RESET" \n", stderr);

    fputs(ChkMat.Is_Hardening_modulus ? 
    ""MAGENTA"[Hardening-modulus]"RESET" : "GREEN"true"RESET" \n" : 
    ""MAGENTA"[Hardening-modulus]"RESET" : "RED"false"RESET" \n", stderr);

    fputs(ChkMat.Is_Reference_pressure ? 
    ""MAGENTA"[Reference-pressure]"RESET" : "GREEN"true"RESET" \n" :
    ""MAGENTA"[Reference-pressure]"RESET" : "RED"false"RESET" \n",
          stderr);

    fputs(ChkMat.Is_friction_angle ? ""MAGENTA"[Friction-angle]"RESET" : "GREEN"true"RESET" \n"
                                   : ""MAGENTA"[Friction-angle]"RESET" : "RED"false"RESET" \n",
          stderr);

    fputs(ChkMat.Is_dilatancy_angle ? ""MAGENTA"[Dilatancy-angle]"RESET" : "GREEN"true"RESET" \n"
                                    : ""MAGENTA"[Dilatancy-angle]"RESET" : "RED"false"RESET" \n",
          stderr);

    fputs(ChkMat.Is_J2_degradated ? 
    ""MAGENTA"[J2-degradated]"RESET" : "GREEN"true"RESET" \n" :
    ""MAGENTA"[J2-degradated]"RESET" : "RED"false"RESET" \n",
          stderr);

    fputs(ChkMat.Is_kappa_0 ? 
    ""MAGENTA"[kappa-0]"RESET" : "GREEN"true"RESET" \n":
    ""MAGENTA"[kappa-0]"RESET" : "RED"false"RESET" \n",
          stderr);

    STATUS = EXIT_SUCCESS;
  }

  return STATUS;
}

/**********************************************************************/