/*! \file MPM.h
    \brief File with the prototype of the material point functions/utilities
*/

#ifndef _MPM_H_
#define _MPM_H_

/*******************************************************/

/* Boundary conditions */
Curve BcDirichlet(char *);
void imposse_NodalMomentum(Mesh, Matrix, int);
void imposse_NodalVelocity(Mesh, Matrix, int);
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
double update_Density(double, double, Tensor);
Tensor compute_Stress(Tensor, Tensor, Material);
void update_LocalState(Matrix, GaussPoint, Mesh, double);
double compute_InternalEnergy(Tensor, Tensor);
Matrix compute_InternalForces(Matrix, GaussPoint,Mesh);
Matrix compute_BodyForces(Matrix, GaussPoint, Mesh, int);
Matrix compute_ContacForces(Matrix, GaussPoint, Mesh, int);

/*******************************************************/

/* Equilibrium scheme */
Matrix compute_equilibrium_U(Matrix,GaussPoint,Mesh,double);

/*******************************************************/

/* Forward-Euler */
Matrix compute_NodalMomentumMass(GaussPoint, Mesh);
Matrix compute_NodalVelocity(Mesh, Matrix);
void update_NodalMomentum(Mesh, Matrix, Matrix);
void update_Particles_FE(GaussPoint, Mesh, Matrix, Matrix, double);

/*******************************************************/

/* Explicit Predictor-Corrector */
Matrix compute_NodalMass(GaussPoint, Mesh);
Matrix compute_VelocityPredictor(GaussPoint, Mesh, Matrix,
				 Matrix, Time_Int_Params,double);
Matrix compute_VelocityCorrector(Mesh, Matrix, Matrix,
				 Matrix, Time_Int_Params,double);
void update_Particles_PCE(GaussPoint, Mesh, Matrix,
			  Matrix, Matrix, double);

/*******************************************************/

/* Generalized-alpha */
void GA_UpdateNodalKinetics(Mesh, Matrix, Matrix, Time_Int_Params);
Matrix GetNodalKinetics(GaussPoint, Mesh);
Matrix GetNodalVelocityDisplacement(GaussPoint, Mesh);
void update_Particles_GA(GaussPoint, Mesh, Matrix, Time_Int_Params);

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
