#include <stdio.h>
#include <stdlib.h>
#include "ElementTools.h"

int main(void){

  Element * E2D = (Element *)Allocate_Array(1,sizeof(Element));
  GaussPoint * GP_e = (GaussPoint *)Allocate_Array(1,sizeof(GaussPoint));

  
  double ** X_g = (double **)Allocate_Matrix(4,2,sizeof(double));
  Vector * X_e = (Vector *)Allocate_Array(1,sizeof(Vector));
  
  int * I_g = (int *)Allocate_Array(4,sizeof(int));
  /* int * N_id = (int *)Allocate_Array(4,sizeof(int)); */  
  double ** Inv_F_e = (double **)Allocate_Matrix(2,2,sizeof(double));

  char ElementType [2] = "Q4";

  X_e[0].x = 0.5;
  X_e[0].y = 0.5;
  
  X_g[0][0] = 0;
  X_g[0][1] = 0;
  X_g[1][0] = 1;
  X_g[1][1] = 0;
  X_g[2][0] = 4;
  X_g[2][1] = 4;
  X_g[3][0] = 0.5;
  X_g[3][1] = 1;

  I_g[0] = 0;
  I_g[1] = 1;
  I_g[2] = 2;
  I_g[2] = 3;
  
  Initialize_Element(E2D,ElementType,X_g,I_g);
  Initialize_GP(GP_e,X_e);

  double ** B_ref = (double **)Allocate_Matrix(E2D->NumberNodes,
					       E2D->NumberDOF,
					       sizeof(double));
  
  double ** F_e = (double **)Allocate_Matrix(E2D->NumberDOF,
					   E2D->NumberDOF,
					   sizeof(double));

  double ** C_e = (double **)Allocate_Matrix(E2D->NumberDOF,
					     E2D->NumberDOF,
					     sizeof(double));
  
  Get_Deformation_Gradient(GP_e,E2D,B_ref,F_e);
  printf("%f ; %f ; %f ; %f \n",
  	 F_e[0][0],
  	 F_e[0][1],
  	 F_e[1][0],
  	 F_e[1][1]);

  double J_e = Get_Determinant(F_e,2);
  printf("%f \n",J_e);

  Get_Lagrangian_CG_Tensor(GP_e,E2D,C_e);
  printf("%f ; %f ; %f ; %f \n",
  	 C_e[0][0],
  	 C_e[0][1],
  	 C_e[1][0],
  	 C_e[1][1]);
  
  Get_Inverse(F_e,Inv_F_e,2);
  printf("%f ; %f ; %f ; %f \n",
  	 Inv_F_e[0][0],
  	 Inv_F_e[0][1],
  	 Inv_F_e[1][0],
  	 Inv_F_e[1][1]);
  
  
  return 0;
}
