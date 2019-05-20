
#ifndef TypeDefinitions
#define TypeDefinitions
#endif

#ifndef GlobalVariables
#define GlobalVariables
#endif


GaussPoint Initialize_GP_Mesh(Matrix,Matrix,Element);

void LocateGP(GaussPoint,Element);

Matrix GetMassMatrix_L(Element,GaussPoint);

void GaussPointsToMesh(Element,GaussPoint,
		       Matrix,Matrix,Matrix);

void MeshToGaussPoints(Element,GaussPoint,
		       Matrix,Matrix,Matrix);

/* void Get_Lagrangian_CG_Tensor(GaussPoint *); */

/* void Get_Eulerian_CG_Tensor(GaussPoint *); */
