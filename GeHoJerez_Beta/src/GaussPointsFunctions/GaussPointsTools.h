
#ifndef TypeDefinitions
#define TypeDefinitions
#endif

#ifndef GlobalVariables
#define GlobalVariables
#endif


GaussPoint Initialize_GP_Mesh(char *,Matrix, double, Matrix);

Matrix GetMassMatrix_L(Mesh,GaussPoint);

void GaussPointsToMesh(Mesh,GaussPoint,
		       Matrix,Matrix,Matrix);

void MeshToGaussPoints(Mesh,GaussPoint,
		       Matrix,Matrix,Matrix);

/* void Get_Lagrangian_CG_Tensor(GaussPoint *); */

/* void Get_Eulerian_CG_Tensor(GaussPoint *); */
