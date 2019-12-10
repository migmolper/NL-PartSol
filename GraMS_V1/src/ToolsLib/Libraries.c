#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include "TypeDefinitions.h"

/*********************************************************************/

MatLib MatrixOperators(void){
  
#include "../Matlib/Matlib.h"
  
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

/*********************************************************************/

ConstLib Contitutive(void){

#include "../Constitutive/Constitutive.h"
  
  /* Define variable with the functions */
  ConstLib CL;

  CL.LE2D = LinearElastic2D;
  
  return CL;
}

/*********************************************************************/

LME Load_LME(void){

#include "../ShapeFun/ShapeFun.h"
   
  /* Define variable with the functions */
  LME LME_shpf;

  /* Asign functions to the library */
  LME_shpf.lambda = LME_lambda;
  LME_shpf.fa = LME_fa;
  LME_shpf.pa = LME_pa;
  LME_shpf.r= LME_r;
  LME_shpf.J= LME_J;
  LME_shpf.dpa = LME_dpa;
  
  /* Return the library */
  return LME_shpf;
}

/*********************************************************************/

