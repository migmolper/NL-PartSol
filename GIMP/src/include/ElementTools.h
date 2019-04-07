#include "TypeDefinitions.h"

double Get_Determinant(Tensor);

Tensor Get_Inverse(Tensor);

void Initialize_Element(Element *, char *, double **, int * );

void Initialize_GP(GaussPoint *, Vector *);

void Get_RefDeformation_Gradient(GaussPoint *, Element * );

void Get_Deformation_Gradient(Element *, GaussPoint *);

void Get_Lagrangian_CG_Tensor(GaussPoint *);

void Get_Eulerian_CG_Tensor(GaussPoint *); 

double ** Get_Stiffness_Matrix(double **, GaussPoint *, Element * );
