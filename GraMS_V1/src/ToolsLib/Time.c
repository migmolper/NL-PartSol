#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include "TypeDefinitions.h"

/*********************************************************************/

/* Time measure function declaration */

void TIC_toc(clock_t start_t){

  /* Start of tic-toc */
  start_t = clock();
  printf("Starting of the program, start_t = %ld\n", start_t);
  
}

/*********************************************************************/

void tic_TOC(clock_t start_t){

  /* Define internal parameters */
  clock_t end_t, total_t;

  /* End of tic-toc */
  end_t = clock();
  printf("End of the big loop, end_t = %ld\n", end_t);
  total_t = (double)(end_t - start_t) / CLOCKS_PER_SEC;
  printf("Total time taken by CPU: %ld\n", total_t  );
  
}

/*********************************************************************/
