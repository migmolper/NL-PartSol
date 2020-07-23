/*! \file Particles.h
    \brief File with the prototype of the material point functions/utilities
*/

#ifndef _PARTICLES_H_
#define _PARTICLES_H_

/*******************************************************/

/*!

*/
Matrix body_loads__Particles__(Load *, int, int, int);
/*******************************************************/

/*!

*/
Matrix contact_loads__Particles__(Load *, int, int, int);
/*******************************************************/

/*!

*/
Tensor rate_inifinitesimal_Strain__Particles__(Matrix, Matrix);
/*******************************************************/

/*!

*/
Tensor infinitesimal_Strain__Particles__(Tensor, Tensor, double);
/*******************************************************/

/*!

*/
Tensor increment_Deformation_Gradient__Particles__(Matrix, Matrix);
/*******************************************************/

/*!

*/
void   update_Deformation_Gradient_n1__Particles__(Tensor, Tensor, Tensor);
/*******************************************************/

/*!

*/
Tensor right_Cauchy_Green__Particles__(Tensor);
/*******************************************************/

/*!

*/
Tensor strain_Green_Lagrange__Particles__(Tensor);
/*******************************************************/

/*!

*/
double update_density__Particles__(double, double, Tensor);
/*******************************************************/

/*!

*/
Tensor explicit_integration_stress__Particles__(Tensor, Tensor, Material);
/*******************************************************/

/*!

*/
Tensor configurational_midpoint_integration_Stress__Particles__(Tensor,Tensor, Tensor, Material);
/*******************************************************/

/*!
  
*/
Tensor average_strain_integration_Stress__Particles__(Tensor, Tensor, Tensor, Material);
/*******************************************************/

/*!

*/
Tensor average_itegration_Stress__Particles__(Tensor, Tensor, Tensor, Material);
/*******************************************************/

/*!

*/
double internal_energy__Particles__(Tensor, Tensor);
/*******************************************************/

/*!

*/
void initial_position__Particles__(Matrix, Mesh, int);
/*******************************************************/

/*!

*/
void local_search__Particles__(GaussPoint, Mesh);
/*******************************************************/

/*!

*/
Element nodal_set__Particles__(int, ChainPtr, int);
/*******************************************************/

/*!

*/
int inout_element__Particles__(int, Matrix, ChainPtr, Mesh);
/*******************************************************/

/*!

*/
void asign_to_nodes__Particles__(int, ChainPtr, Mesh);
/*******************************************************/

#endif
