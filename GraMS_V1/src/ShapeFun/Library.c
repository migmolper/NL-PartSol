#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "../GRAMS/grams.h"

/*********************************************************************/

SHPF ShapeFunLib(void){
   
  /* Define variable with the functions */
  SHPF shpf;

  /* Asign functions to the library */
  shpf.N = LME_pa;
  shpf.dN = LME_dpa;
  
  /* Return the library */
  return shpf;
}

/*********************************************************************/
