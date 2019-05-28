
#ifndef GlobalVariables
#define GlobalVariables
#endif

#ifndef TypeDefinitions
#define TypeDefinitions
#endif

Mesh RectangularMesh(double, double,
		     double, double,
		     double, double, char *);
Matrix GetNodalValuesFromGP(GaussPoint, Mesh, char LisOfFields[]);
