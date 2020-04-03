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
Matrix Q4_N(Matrix);
Matrix Q4_dN_Ref(Matrix);
Matrix Q4_F_Ref(Matrix,Matrix);
Matrix Q4_dN(Matrix,Matrix);
Matrix Q4_Xi_to_X(Matrix,Matrix);
void Q4_X_to_Xi(Matrix,Matrix,Matrix);

/*******************************************************/

/* GIMP shape functions */
void uGIMP_Initialize(GaussPoint, Mesh);
double uGIMP_Sip(double, double, double);
double uGIMP_dSip(double, double, double);
Matrix uGIMP_N(Matrix, Matrix, double);
Matrix uGIMP_dN(Matrix, Matrix, double);
ChainPtr uGIMP_Tributary_Nodes(Matrix, int,
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
Matrix compute_ShapeFunction(Element, GaussPoint, Mesh );
Matrix compute_ShapeFunction_Gradient(Element, GaussPoint, Mesh);

/*******************************************************/
