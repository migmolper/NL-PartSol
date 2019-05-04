#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "ToolsLib/TypeDefinitions.h"
#include "ToolsLib/Utils.h"
#include "InOutFun/InOutFun.h"
#include "Constitutive/Constitutive.h"
#include "ElementsFunctions/ElementTools.h"

#include "GaussPointsFunctions/GaussPointsTools.h"
#include "Solvers/Solvers.h"

void main(int argc, char *argv[])
/*
  Inputs parameters :
  * Mesh file
  * Data file
*/
{
  /* Check command-line arguments */
  if(argc == 1){
    perror("Error in main(), insuficient number of input files !");
    exit(0);
  }

  /* Read mesh data */
  ReadGidMesh(argv[1]);
  
  /* Initialize the element mesh */
  char * TypeElement = "L2";
  Element ElementMesh = Initialize_Element(0,TypeElement);
  
  /* Read the input fields as a CSV */
  Matrix InputFields = Read_CSV(argv[2], 5);

  /* Read the .dat file */
  ReadDatFile(argv[3]);
    
  /* Read the data file */  /* Physical parameters */
  double YoungModulus = 1;
  
  /* Get material Constitutive matrix */
  Matrix D = LinearElastic1D(YoungModulus);

  /* Define Gauss point mesh */
  GaussPoint GP_Mesh  = Initialize_GP_Mesh(InputFields,D);

  /* Localize all the Gauss-Points */
  UpdateElementLocationGP(GP_Mesh);
    
  exit(EXIT_SUCCESS);
  
}
