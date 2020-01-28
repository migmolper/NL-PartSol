#include <stdbool.h> 

#ifndef TypeDefinitions
#define TypeDefinitions
#endif

/*******************************************************/

/* 1D of two nodes shape functions */
Matrix L2(Matrix);
Matrix dL2(Matrix);
Matrix Get_F_Ref_L2(Matrix,Matrix);
Matrix Get_X_GC_L2(Matrix,Matrix);

/*******************************************************/

/* Triangle of three nodes shape functions */
Matrix T3(Matrix);
Matrix dT3(Matrix);
Matrix Get_F_Ref_T3(Matrix,Matrix);
Matrix Get_dNdX_T3(Matrix,Matrix);
Matrix Get_X_GC_T3(Matrix,Matrix);
void Get_X_EC_T3(Matrix,Matrix,Matrix);

/*******************************************************/

/* Quadrilateral of four nodes shape functions */
void Q4_Initialize(GaussPoint, Mesh);
Matrix Q4(Matrix);
Matrix dQ4(Matrix);
Matrix Get_F_Ref_Q4(Matrix,Matrix);
Matrix Get_dNdX_Q4(Matrix,Matrix);
Matrix Get_X_GC_Q4(Matrix,Matrix);
void Get_X_EC_Q4(Matrix,Matrix,Matrix);

/*******************************************************/

/* GIMP shape functions */
void GIMP_Initialize(GaussPoint, Mesh);
double uGIMP(double, double, double);
double d_uGIMP(double, double, double);
Matrix GIMP_2D(Matrix, Matrix, double);
Matrix dGIMP_2D(Matrix, Matrix, double);
ChainPtr Tributary_Nodes_GIMP(Matrix, int,
			      Matrix, Mesh);

/*******************************************************/

/* LME shape functions */
void LME_Initialize(GaussPoint, Mesh);
Matrix LME_Beta(Matrix, Matrix, double);
Matrix LME_lambda_NR(Matrix, Matrix, Matrix);
double LME_fa(Matrix, Matrix, Matrix);
Matrix LME_p(Matrix, Matrix, Matrix);
Matrix LME_r(Matrix, Matrix);
Matrix LME_J(Matrix, Matrix, Matrix);
Matrix LME_dp(Matrix, Matrix);
ChainPtr LME_Tributary_Nodes(Matrix, Matrix, int, Mesh);

/*******************************************************/

/* Operators */
Matrix Get_Operator(char *, Element, GaussPoint, Mesh); 
Matrix Get_B_GP(Matrix);

/*******************************************************/
