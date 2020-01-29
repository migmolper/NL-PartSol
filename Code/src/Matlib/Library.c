#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "grams.h"

MatLib MatrixOperators(void){
    
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
