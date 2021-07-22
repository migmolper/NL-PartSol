/*! \file Q4.h
  \brief Linear quadrilateral shape function

  (3)     (2)  
  o-------o   
  |       |   
  |       |   
  o-------o   
  (0)     (1)  

*/

#ifndef _Q4_H_
#define _Q4_H_

void   initialize__Q4__(Particle, Mesh);
Matrix N__Q4__(Matrix);
Matrix dN__Q4__(Matrix,Matrix);
Matrix dN_Ref__Q4__(Matrix);
void   X_to_Xi__Q4__(Matrix,Matrix,Matrix);
bool   in_out__Q4__(Matrix,Matrix);
void   element_to_particles__Q4__(Matrix, Mesh, int);
double min_DeltaX__Q4__(ChainPtr, Matrix);
double volume__Q4__(Matrix);
void   local_search__Q4__(Particle, Mesh);
double compute_Jacobian_patch__Q4__(int,Particle,ChainPtr *,double *,double *);
#endif
