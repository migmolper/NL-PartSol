/*! \file T3.h
  \brief Linear triangle shape function

  (2) 
  o      
  |\     
  | \    
  o--o   
  (0) (1) 
*/

#ifndef _T6_H_
#define _T6_H_

void   initialize__T6__(Particle, Mesh);
Matrix N__T6__(Matrix);
Matrix dN__T6__(Matrix,Matrix);
Matrix dN_Ref__T6__(Matrix);
bool   in_out__T6__(Matrix,Matrix);
void   element_to_particles__T6__(Matrix, Mesh, int);
double min_DeltaX__T6__(ChainPtr, Matrix);
double volume__T6__(Matrix);
void   local_search__T6__(Particle, Mesh);
#endif
