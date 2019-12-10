#include <stdbool.h> 

#ifndef TypeDefinitions
#define TypeDefinitions
#endif

/* 1D of two nodes shape functions */
Matrix L2(Matrix);
Matrix dL2(Matrix);
Matrix Get_F_Ref_L2(Matrix,Matrix);
Matrix Get_X_GC_L2(Matrix,Matrix);

/* Triangle of three nodes shape functions */
Matrix T3(Matrix);
Matrix dT3(Matrix);
Matrix Get_F_Ref_T3(Matrix,Matrix);
Matrix Get_dNdX_T3(Matrix,Matrix);
Matrix Get_X_GC_T3(Matrix,Matrix);
void Get_X_EC_T3(Matrix,Matrix,Matrix);

/* Quadrilateral of four nodes shape functions */
Matrix Q4(Matrix);
Matrix dQ4(Matrix);
Matrix Get_F_Ref_Q4(Matrix,Matrix);
Matrix Get_dNdX_Q4(Matrix,Matrix);
Matrix Get_X_GC_Q4(Matrix,Matrix);
void Get_X_EC_Q4(Matrix,Matrix,Matrix);

/* GIMP shape functions */
double uGIMP(double, double, double);
double d_uGIMP(double, double, double);
Matrix GIMP_2D(Matrix, Matrix, double);
Matrix dGIMP_2D(Matrix, Matrix, double);
ChainPtr Tributary_Nodes_GIMP(Matrix, int,
			      Matrix, Mesh);

/* LME shape functions */
Matrix LME_lambda(Matrix, Matrix,
		  double, double);
double LME_fa(Matrix, Matrix,
	      double, double);
Matrix LME_pa(Matrix, Matrix,
	      double, double);
Matrix LME_r(Matrix, Matrix);
Matrix LME_J(Matrix, Matrix, Matrix);
Matrix LME_dpa(Matrix, Matrix);
ChainPtr LME_Tributary_Nodes(Matrix, int,
			     Mesh, double);

/*******************************************************/

/* Library */

typedef struct {

  Matrix (*N)(Matrix, Matrix, double);
  Matrix (*dN)(Matrix, Matrix);
  
} SHPF;

SHPF ShapeFunLib(void);


