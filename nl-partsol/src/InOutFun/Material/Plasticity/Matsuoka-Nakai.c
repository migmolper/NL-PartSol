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

static char Error_message[MAXW];

typedef struct {

  bool Is_rho;                       // Reference fensity
  bool Is_E;                         // Young modulus
  bool Is_nu;                        // Poisson cefficient
  bool Is_Reference_pressure;      // Reference pressure
  bool Is_a1_Hardening_Borja;        // Hardening parameter
  bool Is_a2_Hardening_Borja;        // Hardening parameter
  bool Is_a3_Hardening_Borja;        // Hardening parameter
  bool Is_atmospheric_pressure;      // Reference pressure
  bool Is_friction_angle;           // Friction angle
  bool Is_dilatancy_parameter;        // dilatancy parameter
  bool Is_Reference_plastic_strain; // Reference Plastic Strain
  bool Is_kappa_0;                   // Initial yield
  bool Is_J2_degradated; // Degradation limit
  bool Is_cohesion;

} Check_Material;

static void standard_error();
static Check_Material Initialise_Check_Material();
static int __check_material(Material *, Check_Material, int);

/**********************************************************************/

int Define_Matsuoka_Nakai(Material * MN_Material,FILE *Simulation_file, char *Material_Model,
                           int Material_Idx) {

  int STATUS = EXIT_SUCCESS;
  int Ndim = NumberDimensions;

  /* Variables for reading purposes */  
  char Parameter_line[MAXC] = {0};
  char *Parameter_pars[MAXW] = {NULL};
  int Parser_status;

  /* Check variables for sintax */
  bool Is_Open = false;
  bool Is_Close = false;
  Check_Material ChkMat = Initialise_Check_Material();

  /* Default parameters */
  (*MN_Material).kappa_0 = 0.0; // Cohesionless
  (*MN_Material).ReferencePressure = 0.0; 
  (*MN_Material).J2_degradated = 0.0;
  (*MN_Material).Cohesion = 0.0;
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
      (*MN_Material).rho = atof(Parameter_pars[1]);
    }
    /**************************************************/
    else if (strcmp(Parameter_pars[0], "E") == 0) {
      ChkMat.Is_E = true;
      (*MN_Material).E = atof(Parameter_pars[1]);
    }
    /**************************************************/
    else if (strcmp(Parameter_pars[0], "nu") == 0) {
      ChkMat.Is_nu = true;
      (*MN_Material).nu = atof(Parameter_pars[1]);
    }
    /**************************************************/
    else if (strcmp(Parameter_pars[0], "alpha") == 0) {
      ChkMat.Is_dilatancy_parameter = true;
      (*MN_Material).alpha_Hardening_Borja = atof(Parameter_pars[1]);
    }
    /**************************************************/
    else if (strcmp(Parameter_pars[0], "a1") == 0) {
      ChkMat.Is_a1_Hardening_Borja = true;
      (*MN_Material).a_Hardening_Borja[0] = atof(Parameter_pars[1]);
    }
    /**************************************************/
    else if (strcmp(Parameter_pars[0], "a2") == 0) {
      ChkMat.Is_a2_Hardening_Borja = true;
      (*MN_Material).a_Hardening_Borja[1] = atof(Parameter_pars[1]);
    }
    /**************************************************/
    else if (strcmp(Parameter_pars[0], "a3") == 0) {
      ChkMat.Is_a3_Hardening_Borja = true;
      (*MN_Material).a_Hardening_Borja[2] = atof(Parameter_pars[1]);
    }
    /**************************************************/
    else if (strcmp(Parameter_pars[0], "Atmospheric-pressure") == 0) {
      ChkMat.Is_atmospheric_pressure = true;
      (*MN_Material).atmospheric_pressure = atof(Parameter_pars[1]);
    }
    /**************************************************/
    else if (strcmp(Parameter_pars[0], "Reference-pressure") == 0) {
      ChkMat.Is_Reference_pressure = true;
      (*MN_Material).ReferencePressure = atof(Parameter_pars[1]);
    }
    /**************************************************/
    else if (strcmp(Parameter_pars[0], "Friction-angle") == 0) {
      ChkMat.Is_friction_angle = true;
      (*MN_Material).phi_Frictional = atof(Parameter_pars[1]);

      if(((*MN_Material).phi_Frictional > 90.0) && ((*MN_Material).phi_Frictional < 0.0))
      {
        fprintf(stderr, ""RED" Invalid value of the [Friction-angle] "RESET" \n");
        return EXIT_FAILURE;
      }
    }
    /**************************************************/
    else if (strcmp(Parameter_pars[0], "EPS-0") == 0) {
      ChkMat.Is_Reference_plastic_strain = true;
      (*MN_Material).Plastic_Strain_0 = atof(Parameter_pars[1]);
    }    
    /**************************************************/
    else if (strcmp(Parameter_pars[0], "kappa-0") == 0) {
      ChkMat.Is_kappa_0 = true;
      (*MN_Material).kappa_0 = atof(Parameter_pars[1]);
    }    
    /**************************************************/
    else if (strcmp(Parameter_pars[0], "Cohesion") == 0) {
      ChkMat.Is_cohesion = true;
      (*MN_Material).Cohesion = atof(Parameter_pars[1]);
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

  strcpy((*MN_Material).Type, Material_Model);

  if(ChkMat.Is_friction_angle == true)
  {
    double rad_friction_angle = (PI__MatrixLib__ / 180.0) * (*MN_Material).phi_Frictional;  
    double c_cotphi = (*MN_Material).Cohesion/tan(rad_friction_angle);
    double a1 = (*MN_Material).a_Hardening_Borja[0];
    double a2 = (*MN_Material).a_Hardening_Borja[1];
    double a3 = (*MN_Material).a_Hardening_Borja[2];
    double f, df;
    double I1 = 3*((*MN_Material).ReferencePressure - c_cotphi);
    double EPS_0 = 0.0;
    double kappa_0 = 8.0 * sin(rad_friction_angle)*sin(rad_friction_angle) / (1.0 - sin(rad_friction_angle)*sin(rad_friction_angle));
    int iter = 0;

    f = kappa_0 - a1 * EPS_0 * exp(a2*I1) * exp(-a3 * EPS_0);
    while (fabs(f) > TOL_Radial_Returning) {
      iter++;
      df = (a3 * EPS_0 - 1) * a1 * exp(a2*I1) * exp(-a3 * EPS_0);
      EPS_0 += - f / df;
      f = kappa_0 - a1 * EPS_0 * exp(a2*I1) * exp(-a3 * EPS_0);
      if (iter > 10) {
        fprintf(stderr, ""RED" Iter > 10 "RESET" \n");
        return EXIT_FAILURE;
      }
    }

    (*MN_Material).kappa_0 = kappa_0;
    (*MN_Material).Plastic_Strain_0 = EPS_0;  
  }
  else if((ChkMat.Is_Reference_plastic_strain = true) 
  && (ChkMat.Is_kappa_0 == true))
  {
    (*MN_Material).phi_Frictional = asin(sqrt((*MN_Material).kappa_0/((*MN_Material).kappa_0+8)));
  }
  else {
    fprintf(stderr, "%s : %s \n", "Error in GramsMaterials()",
            "Some parameters are missed for Matsuoka-Nakai initialization");

    fputs(ChkMat.Is_Reference_plastic_strain ? 
    ""MAGENTA"[EPS-0]"RESET" : "GREEN"true"RESET" \n":
    ""MAGENTA"[EPS-0]"RESET" : "RED"false"RESET" \n",stderr);

    fputs(ChkMat.Is_friction_angle ? 
    ""MAGENTA"[Friction-angle]"RESET" : "GREEN"true"RESET" \n":
    ""MAGENTA"[Friction-angle]"RESET" : "RED"false"RESET" \n",stderr);

    fputs(ChkMat.Is_kappa_0 ? 
    ""MAGENTA"[kappa-0]"RESET" : "GREEN"true"RESET" \n":
    ""MAGENTA"[kappa-0]"RESET" : "RED"false"RESET" \n",stderr);

    return EXIT_FAILURE;
  }
  

  STATUS = __check_material(MN_Material, ChkMat, Material_Idx);
  if(STATUS == EXIT_FAILURE)
  {
        fprintf(stderr, ""RED" Error in __check_DP_Material() "RESET" \n");
    return EXIT_FAILURE;
  }

  /* Return outputs */
  return EXIT_SUCCESS;
}

/***************************************************************************/

static Check_Material Initialise_Check_Material() {
  Check_Material ChkMat;

  ChkMat.Is_rho = false;
  ChkMat.Is_E = false;
  ChkMat.Is_nu = false;
  ChkMat.Is_Reference_pressure = false;
  ChkMat.Is_a1_Hardening_Borja = false;
  ChkMat.Is_a2_Hardening_Borja = false;
  ChkMat.Is_a3_Hardening_Borja = false;
  ChkMat.Is_atmospheric_pressure = false;
  ChkMat.Is_kappa_0 = false;
  ChkMat.Is_dilatancy_parameter = false;
  ChkMat.Is_J2_degradated = false;
  ChkMat.Is_kappa_0 = false;
  ChkMat.Is_Reference_plastic_strain = false;
  ChkMat.Is_cohesion = false;

  return ChkMat;
}

/**********************************************************************/

static int __check_material(Material * MN_Material, Check_Material ChkMat, int Idx) {

  int STATUS = EXIT_SUCCESS;

  if (ChkMat.Is_rho && ChkMat.Is_E && ChkMat.Is_nu && ChkMat.Is_atmospheric_pressure &&
      ChkMat.Is_dilatancy_parameter && ChkMat.Is_a1_Hardening_Borja &&
      ChkMat.Is_a2_Hardening_Borja && ChkMat.Is_a3_Hardening_Borja) {

    printf("\t -> %s \n", "Matsuoka-Nakai material");

    printf("\t \t -> %s : %f \n", ""MAGENTA"[rho]"RESET"", (*MN_Material).rho);

    printf("\t \t -> %s : %f \n", ""MAGENTA"[E]"RESET"", (*MN_Material).E);

    printf("\t \t -> %s : %f \n", ""MAGENTA"[nu]"RESET"", (*MN_Material).nu);
 
    printf("\t \t -> %s : %f \n", ""MAGENTA"[Atmospheric-pressure]"RESET"", (*MN_Material).atmospheric_pressure);

    printf("\t \t -> %s : %f \n", ""MAGENTA"[Friction-angle]"RESET"", (*MN_Material).phi_Frictional);

    printf("\t \t -> %s : %f \n", ""MAGENTA"[Kappa-0]"RESET"", (*MN_Material).kappa_0);

    printf("\t \t -> %s : %f \n", ""MAGENTA"[EPS-0]"RESET"", (*MN_Material).Plastic_Strain_0);

    printf("\t \t -> %s : %f \n", ""MAGENTA"[alpha]"RESET"",(*MN_Material).alpha_Hardening_Borja);

    printf("\t \t -> %s : %f \n", ""MAGENTA"[a1]"RESET"",(*MN_Material).a_Hardening_Borja[0]);

    printf("\t \t -> %s : %f \n", ""MAGENTA"[a2]"RESET"",(*MN_Material).a_Hardening_Borja[1]);

    printf("\t \t -> %s : %f \n", ""MAGENTA"[a3]"RESET"",(*MN_Material).a_Hardening_Borja[2]);

    printf("\t \t -> %s : %f \n", ""MAGENTA"[J2-degradated]"RESET"", (*MN_Material).J2_degradated);

    printf("\t \t -> %s : %f \n", ""MAGENTA"[Reference-pressure]"RESET"", (*MN_Material).ReferencePressure);

    printf("\t \t -> %s : %f \n", ""MAGENTA"[Cohesion]"RESET"", (*MN_Material).Cohesion);


  } else {
    fprintf(stderr, "%s : %s \n", "Error in GramsMaterials()",
            "Some parameter is missed for Matsuoka-Nakai");

    fputs(ChkMat.Is_rho ? 
    ""MAGENTA"[rho]"RESET" : "GREEN"true"RESET" \n" : 
    ""MAGENTA"[rho]"RESET" : "RED"false"RESET" \n", stderr);
    
    fputs(ChkMat.Is_E ? 
    ""MAGENTA"[E]"RESET" : "GREEN"true"RESET" \n" :
    ""MAGENTA"[E]"RESET" : "RED"false"RESET" \n",stderr);

    fputs(ChkMat.Is_nu ? 
    ""MAGENTA"[nu]"RESET" : "GREEN"true"RESET" \n"  : 
    ""MAGENTA"[nu]"RESET" : "RED"false"RESET" \n",stderr);
    
    fputs(ChkMat.Is_atmospheric_pressure ? 
    ""MAGENTA"[Atmospheric-pressure]"RESET" : "GREEN"true \n":
    ""MAGENTA"[Atmospheric-pressure]"RESET" : "RED"false \n",stdout);

    fputs(ChkMat.Is_dilatancy_parameter ? 
    ""MAGENTA"[alpha]"RESET" : "GREEN"true"RESET" \n":
    ""MAGENTA"[alpha]"RESET" : "RED"false"RESET" \n",stdout);

    fputs(ChkMat.Is_a1_Hardening_Borja ? 
    ""MAGENTA"[a1]"RESET" : "GREEN"true"RESET" \n" :
    ""MAGENTA"[a1]"RESET" : "RED"false"RESET" \n",stdout);

    fputs(ChkMat.Is_a2_Hardening_Borja ? 
    ""MAGENTA"[a2]"RESET" : "GREEN"true"RESET" \n": 
    ""MAGENTA"[a2]"RESET" : "RED"false"RESET" \n",stdout);

    fputs(ChkMat.Is_a3_Hardening_Borja ? 
    ""MAGENTA"[a3]"RESET" : "GREEN"true"RESET" \n": 
    ""MAGENTA"[a3]"RESET" : "RED"false"RESET" \n",stdout);

    fputs(ChkMat.Is_J2_degradated ? 
    ""MAGENTA"[J2-degradated]"RESET" : "GREEN"true"RESET" \n" :
    ""MAGENTA"[J2-degradated]"RESET" : "RED"false"RESET" \n",stderr);
    
    return EXIT_FAILURE;
  }

  return EXIT_SUCCESS;
}

/**********************************************************************/

static void standard_error() {
  fprintf(stderr, "%s : %s !!! \n", "Error in Define-Material()",
          Error_message);
  exit(EXIT_FAILURE);
}

/***************************************************************************/
