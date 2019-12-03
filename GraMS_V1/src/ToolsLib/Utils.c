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

void tic_TOC(clock_t start_t){

  /* Define internal parameters */
  clock_t end_t, total_t;

  /* End of tic-toc */
  end_t = clock();
  printf("End of the big loop, end_t = %ld\n", end_t);
  total_t = (double)(end_t - start_t) / CLOCKS_PER_SEC;
  printf("Total time taken by CPU: %ld\n", total_t  );
  
}


MatLib MatrixOperators(void){
  
#include "../MathTools/MathTools.h"
  
  /* Define variable with the functions */
  MatLib LibraryMatOp;

  /* Asign functions to the library */
  LibraryMatOp.Alloc = MatAlloc;
  LibraryMatOp.AllocZ = MatAllocZ;
  LibraryMatOp.Assign = MatAssign;
  LibraryMatOp.FreeMat= FreeMat;
  LibraryMatOp.Print= PrintMatrix;
  LibraryMatOp.Copy = CopyMat;
  LibraryMatOp.Norm = Norm_Mat;
  LibraryMatOp.Cond = Cond_Mat;
  LibraryMatOp.Det= Get_Determinant;
  LibraryMatOp.Inv = Get_Inverse;
  LibraryMatOp.Trans = Transpose_Mat;
  LibraryMatOp.Sprod = Scalar_prod;
  LibraryMatOp.Vprod = Vectorial_prod;
  LibraryMatOp.Tprod = Tensorial_prod;
  LibraryMatOp.Incr = Incr_Mat;
  LibraryMatOp.Add = Add_Mat;
  LibraryMatOp.Sub = Sub_Mat;
  LibraryMatOp.Lumped = Get_Lumped_Matrix;

  /* Return the library */
  return LibraryMatOp;
}

ConstLib Contitutive(void){

#include "../Constitutive/Constitutive.h"
  
  /* Define variable with the functions */
  ConstLib CL;

  CL.LE2D = LinearElastic2D;
  
  return CL;
}
