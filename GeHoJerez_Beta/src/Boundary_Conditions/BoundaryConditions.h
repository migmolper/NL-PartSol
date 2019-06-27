
#ifndef TypeDefinitions
#define TypeDefinitions
#endif

BoundaryConditions SetBoundaryConditions(int *, int, int *, char *, char *);

void BCC_Nod_Momentum(Mesh, BoundaryConditions, Matrix, int);

void BCC_GP_Forces(GaussPoint, BoundaryConditions, int);

