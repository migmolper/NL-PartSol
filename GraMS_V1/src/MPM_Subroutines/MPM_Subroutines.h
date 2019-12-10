
#ifndef GlobalVariables
#define GlobalVariables
#endif

#ifndef TypeDefinitions
#define TypeDefinitions
#endif


/* Define and initialize both mesh */
GaussPoint Define_GP_Mesh(char *, double);
GaussPoint InitializeGP(char *, Mesh);
Mesh InitializeMesh(char *);

/* int SearchGaussPoint(int, Matrix, Matrix, Mesh); */
Matrix GetNodalMassMomentum(GaussPoint, Mesh);
Matrix GetNodalVelocity(Mesh, Matrix, Matrix);

void BCC_Nod_VALUE(Mesh, Matrix, int);
void UpdateGaussPointStrain(GaussPoint, Mesh, Matrix);
double UpdateGaussPointDensity(double, double);
void UpdateGaussPointStress(GaussPoint);
Matrix GetNodalForces(GaussPoint, Mesh, int);
void UpdateGridNodalMomentum(Mesh, Matrix, Matrix);
void UpdateVelocityAndPositionGP(GaussPoint, Mesh,
				 Matrix, Matrix, Matrix);

void GetNodalConnectivity(Mesh);
double GetMinElementSize(Mesh);
void GlobalSearchGaussPoints(GaussPoint, Mesh);
void LocalSearchGaussPoints(GaussPoint, Mesh);
Matrix Get_B_GP(Matrix);
