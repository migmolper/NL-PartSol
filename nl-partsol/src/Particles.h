/**
 * @file Particles.h
 * @author Miguel Molinos (@migmolper)
 * @brief File with the prototype of the material point functions/utilities
 * @version 0.1
 * @date 2022-05-31
 * 
 * @copyright Copyright (c) 2022
 * 
 */

#ifndef _PARTICLES_H_
#define _PARTICLES_H_

// clang-format off
#include "Nodes/H8.h"
#include "Nodes/Q4.h"
#include "Nodes/T3.h"
#include "Nodes/T4.h"
// clang-format on

/*******************************************************/

/*!
  \fn Tensor rate_of_deformation__Particles__(Tensor dFdt, Tensor Fm1)  

  \param dFdt Temporal derivative of the deformation gradient
  \param Fm1 Inverse of the deformation gradient  
*/
Tensor rate_of_deformation__Particles__(Tensor, Tensor);
/*******************************************************/


/*!

*/
void initial_position__Particles__(Matrix, Mesh, int);
/*******************************************************/

/*!

*/
Element nodal_set__Particles__(int, ChainPtr, int);
/*******************************************************/

/*!

*/
int search_particle_in_surrounding_elements__Particles__(int, Matrix, ChainPtr, Mesh);
/*******************************************************/

/*!

*/
void asign_to_nodes__Particles__(int, int, int, ChainPtr, Mesh);
/*******************************************************/


#endif
