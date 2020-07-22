/*! \file MPM.h
    \brief File with the prototype of the material point functions/utilities
*/

#ifndef _MPM_H_
#define _MPM_H_

/*
  compute-Loads.c
 */
Matrix Eval_Body_Forces(Load *, int, int, int);
Matrix Eval_Contact_Forces(Load *, int, int, int);


Tensor compute_RateOfStrain(Matrix, Matrix);
Tensor update_Strain(Tensor, Tensor, double);
Tensor compute_increment_Deformation_Gradient(Matrix, Matrix);
void   update_Deformation_Gradient_n1(Tensor, Tensor, Tensor);
Tensor compute_RightCauchyGreen(Tensor);
Tensor compute_LagrangianStrain(Tensor);  
double update_Density(double, double, Tensor);
Tensor compute_Stress(Tensor, Tensor, Material);
double compute_InternalEnergy(Tensor, Tensor);

/*
  Particles utilities 
*/
void GetInitialGaussPointPosition(Matrix, Mesh, int);
void LocalSearchGaussPoints(GaussPoint, Mesh);
void GPinCell(ChainPtr *, ChainPtr *, Matrix, int, double);

Element get_particle_Set(int, ChainPtr, int);
int search_particle_in(int, Matrix, ChainPtr, Mesh);
void asign_particle_to_nodes(int, ChainPtr, Mesh);

/*******************************************************/

#endif
