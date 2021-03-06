/*! \file T3.h
  \brief Linear triangle shape function

  (2) 
  o      
  |\     
  | \    
  o--o   
  (0) (1) 
*/

#ifndef _T3_H_
#define _T3_H_

void   initialize__T3__(Particle, Mesh);
Matrix N__T3__(Matrix);
Matrix dN__T3__(Matrix,Matrix);
Matrix dN_Ref__T3__(Matrix);
void   X_to_Xi__T3__(Matrix,Matrix,Matrix);
bool   in_out__T3__(Matrix,Matrix);
void   element_to_particles__T3__(Matrix, Mesh, int);
double min_DeltaX__T3__(ChainPtr, Matrix);
double volume__T3__(Matrix);
void   local_search__T3__(Particle, Mesh);
double compute_Jacobian_patch__T3__(int,Particle,ChainPtr *,ChainPtr *);
#endif
