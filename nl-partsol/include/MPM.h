/*! \file MPM.h
    \brief File with the prototype of the material point functions/utilities
*/

#ifndef _MPM_H_
#define _MPM_H_

/*******************************************************/

/* Boundary conditions */
Curve BcDirichlet(char *);
Matrix Eval_Body_Forces(Load *, int, int, int);
Matrix Eval_Contact_Forces(Load *, int, int, int);
Matrix compute_Reactions(Mesh, Matrix);
void free_Load(Load);
void free_LoadTable(Load *);
void free_Boundaries(Boundaries);

/*******************************************************/

/* Forces-stress-strain functions */
Tensor compute_RateOfStrain(Matrix, Matrix);
Tensor update_Strain(Tensor, Tensor, double);

Tensor compute_increment_Deformation_Gradient(Matrix, Matrix);
void   update_Deformation_Gradient_n1(Tensor, Tensor, Tensor);
Tensor compute_RightCauchyGreen(Tensor);
Tensor compute_LagrangianStrain(Tensor);
  
double update_Density(double, double, Tensor);
Tensor compute_Stress(Tensor, Tensor, Material);
void   update_LocalState(Matrix, GaussPoint, Mesh, double);
double compute_InternalEnergy(Tensor, Tensor);
Matrix compute_InternalForces(Matrix, GaussPoint,Mesh);
Matrix compute_BodyForces(Matrix, GaussPoint, Mesh, int);
Matrix compute_ContacForces(Matrix, GaussPoint, Mesh, int);

/*******************************************************/

/* Equilibrium scheme */
Matrix compute_equilibrium_U(Matrix,GaussPoint,Mesh,double);

/*******************************************************/

/* Mesh utilities */
void GetInitialGaussPointPosition(Matrix, Mesh, int);
double GetMinElementSize(Mesh);
void GetNodalConnectivity(Mesh);
ChainPtr DiscardElements(ChainPtr, Matrix, Matrix, Mesh);
void LocalSearchGaussPoints(GaussPoint, Mesh);
void GPinCell(ChainPtr *, ChainPtr *, Matrix, int, double);

Element get_particle_Set(int, ChainPtr, int);
Matrix get_set_Coordinates(ChainPtr, Matrix, Matrix);
Matrix get_set_Field(Matrix, Element);

ChainPtr get_locality_of_node(int, Mesh);
int get_closest_node_to(Matrix, ChainPtr, Matrix);
bool InOut_Element(Matrix, ChainPtr, Matrix);
int search_particle_in(int, Matrix, ChainPtr, Mesh);
Matrix ElemCoordinates(ChainPtr, Matrix);
void asign_particle_to_nodes(int, ChainPtr, Mesh);

/*******************************************************/

#endif
