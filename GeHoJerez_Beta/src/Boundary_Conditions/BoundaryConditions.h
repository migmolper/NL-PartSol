
#ifndef TypeDefinitions
#define TypeDefinitions
#endif

void Read_FEM_BCC(char *, Mesh *);

BoundaryConditions SetBCC(int *, int, int *, char *, char *);

void BCC_Nod_Momentum(Mesh, Matrix, int);

void BCC_GP_Forces(GaussPoint, BoundaryConditions, int);

