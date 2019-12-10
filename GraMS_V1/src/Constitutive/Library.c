#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "../GRAMS/TypeDefinitions.h"

/*********************************************************************/

ConstLib Contitutive(void){
 
  /* Define variable with the functions */
  ConstLib CL;

  CL.LE2D = LinearElastic2D;
  
  return CL;
}

/*********************************************************************/
