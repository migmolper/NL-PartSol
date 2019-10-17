#include "../MathTools/MathTools.h"

#ifndef TypeDefinitions
#define TypeDefinitions
#endif

/* Initialize the backgraund mesh */
Mesh InitializeMesh(char * GDF);

/* Functions definitions :*/
int ** GetNodalConnectivity(Mesh FEM_Mesh);
Matrix Get_B_GP(Matrix);
void GetNaturalCoordinates(Matrix,Matrix,Matrix);
void BCC_Nod_VALUE(Mesh, Matrix, int);


/* 1D of two nodes functions */
Matrix L2(Matrix);
Matrix dL2(Matrix);
Matrix Get_GlobalCoordinates_L2(Matrix,Matrix);
Matrix Get_RefDeformation_Gradient_L2(Matrix,Matrix);

/* Triangle of three nodes functions */
Matrix T3(Matrix);
Matrix dT3(Matrix);

/* Quadrilateral of four nodes functions */
Matrix Q4(Matrix);
Matrix dQ4(Matrix);
Matrix Get_F_Ref_Q4(Matrix,Matrix);
Matrix Get_dNdX_Q4(Matrix,Matrix);
Matrix Get_X_GC_Q4(Matrix,Matrix);

/* GIMP shape functions */
double uGIMP(double, double, double, double);
double d_uGIMP(double, double, double, double);
Matrix GIMP_2D(Matrix, Matrix, Matrix, double);
Matrix dGIMP_2D(Matrix, Matrix, Matrix, double);
