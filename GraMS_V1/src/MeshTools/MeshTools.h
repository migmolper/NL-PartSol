#include <stdbool.h> 
#include "../MathTools/MathTools.h"

#ifndef TypeDefinitions
#define TypeDefinitions
#endif

/* Initialize the backgraund mesh */
Mesh InitializeMesh(char * GDF);

/* Functions definitions :*/
void GetNodalConnectivity(Mesh FEM_Mesh);
void GlobalSearchGaussPoints(GaussPoint, Mesh);
void LocalSearchGaussPoints(GaussPoint, Mesh);
Matrix Get_B_GP(Matrix);

/* Boundary conditions functions */
void BCC_Nod_VALUE(Mesh, Matrix, int);

/* Nodes operations */
ChainPtr ArrayToChain(int *, int);
int * ChainToArray(ChainPtr, int);
void FreeChain(ChainPtr);
int LenghtChain(ChainPtr);
bool IsPresentNode (ChainPtr, int);
void PushNodeTop (ChainPtr *, int);
void PopNode (ChainPtr *, int);
ChainPtr CopyChain(ChainPtr);
ChainPtr ChainUnion(ChainPtr, ChainPtr);
ChainPtr ChainIntersection(ChainPtr, ChainPtr);
void printList (ChainPtr);

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
