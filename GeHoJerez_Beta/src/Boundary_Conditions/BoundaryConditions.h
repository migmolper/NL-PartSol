
#ifndef TypeDefinitions
#define TypeDefinitions
#endif

BoundaryConditions SetBCC(int, int *, char *, char *);

void BCC_Nod_Momentum(Mesh, Matrix, int);

void BCC_GP_Forces(GaussPoint, BoundaryConditions, int);

