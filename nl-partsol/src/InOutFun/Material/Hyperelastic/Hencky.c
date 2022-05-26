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

  bool Is_rho;                       // Reference fensity
  bool Is_E;                         // Young modulus
  bool Is_nu;                        // Poisson cefficient

} Check_Material;

static Check_Material __Initialise_Check_Material();
static int __check_material(Material *, Check_Material, int);

/**********************************************************************/

int Define_Hencky(Material *H_Material,FILE *Simulation_file, char *Material_Model, int Material_Idx) {

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
      (*H_Material).rho = atof(Parameter_pars[1]);
    }
    /**************************************************/
    else if (strcmp(Parameter_pars[0], "E") == 0) {
      ChkMat.Is_E = true;
      (*H_Material).E = atof(Parameter_pars[1]);
    }
    /**************************************************/
    else if (strcmp(Parameter_pars[0], "nu") == 0) {
      ChkMat.Is_nu = true;
      (*H_Material).nu = atof(Parameter_pars[1]);
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

  strcpy((*H_Material).Type, Material_Model);


  STATUS = __check_material(H_Material, ChkMat, Material_Idx);
  if(STATUS == EXIT_FAILURE)
  {
        fprintf(stderr, ""RED" Error in __check_material() "RESET" \n");
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

  return ChkMat;
}

/**********************************************************************/

static int __check_material(Material * H_Material,
                                      Check_Material ChkMat, int Idx) {

int STATUS = EXIT_SUCCESS;

  if (ChkMat.Is_rho && ChkMat.Is_E && ChkMat.Is_nu) {


    printf("\t -> %s \n", "Drucker-Prager material");

    printf("\t \t -> %s : %f \n", ""MAGENTA"[rho]"RESET"", (*H_Material).rho);

    printf("\t \t -> %s : %f \n", ""MAGENTA"[E]"RESET"", (*H_Material).E);

    printf("\t \t -> %s : %f \n", ""MAGENTA"[nu]"RESET"", (*H_Material).nu);

  } else {
    
    
    fprintf(stderr,""RED" %s : %s "RESET" \n", "Error in GramsMaterials()",
            "Some parameter is missed for Hencky");
    
    fputs(ChkMat.Is_rho ? 
    ""MAGENTA"[rho]"RESET" : "GREEN"true"RESET" \n" : 
    ""MAGENTA"[rho]"RESET" : "RED"false"RESET" \n", stderr);
    
    fputs(ChkMat.Is_E ? 
    ""MAGENTA"[E]"RESET" : "GREEN"true"RESET" \n" :
    ""MAGENTA"[E]"RESET" : "RED"false"RESET" \n",stderr);

    fputs(ChkMat.Is_nu ? 
    ""MAGENTA"[nu]"RESET" : "GREEN"true"RESET" \n"  : 
    ""MAGENTA"[nu]"RESET" : "RED"false"RESET" \n",stderr);

    STATUS = EXIT_FAILURE;
  }

  return STATUS;
}

/**********************************************************************/