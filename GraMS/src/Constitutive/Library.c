#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "../GRAMS/grams.h"

/*********************************************************************/

ConstLib Contitutive(void){
 
  /* Define variable with the functions */
  ConstLib CL;

  CL.LE = LinearElastic;
  
  return CL;
}

/*********************************************************************/
