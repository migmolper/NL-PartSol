
#ifndef GlobalVariables
#define GlobalVariables
#endif

#ifndef TypeDefinitions
#define TypeDefinitions
#endif

Mesh RectangularMesh(double, double,
		     double, double,
		     double, double, char *);
Matrix GetNodalValuesFromGP(GaussPoint, Mesh, char [MAXW]);

Matrix GetNodalVelocity(Matrix, Matrix);

void GetGaussPointStrainIncrement(GaussPoint, Mesh, Matrix);

void UpdateGaussPointDensity(GaussPoint);

void UpdateGaussPointStressTensor(GaussPoint, Matrix);

void UpdateGridNodalMomentum(Mesh, Matrix, Matrix);

void UpdateVelocityAndPositionGP(GaussPoint, Mesh,
				 Matrix, Matrix, Matrix);
