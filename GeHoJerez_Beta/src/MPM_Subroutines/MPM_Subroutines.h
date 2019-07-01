
#ifndef GlobalVariables
#define GlobalVariables
#endif

#ifndef TypeDefinitions
#define TypeDefinitions
#endif

/* Mesh RectangularMesh(double, double, */
/* 		     double, double, */
/* 		     double, double, char *); */


void GlobalSearchGaussPoints(GaussPoint, Mesh);

void LocalSearchGaussPoints(GaussPoint, Mesh);

Matrix GetNodalMassMomentum(GaussPoint, Mesh);

Matrix GetNodalVelocity(Mesh, Matrix, Matrix);

void UpdateGaussPointStrain(GaussPoint, Mesh, Matrix);

double UpdateGaussPointDensity(double, double);

void UpdateGaussPointStress(GaussPoint);

Matrix GetNodalForces(GaussPoint, Mesh, int);

void UpdateGridNodalMomentum(Mesh, Matrix, Matrix);

void UpdateVelocityAndPositionGP(GaussPoint, Mesh,
				 Matrix, Matrix, Matrix);
