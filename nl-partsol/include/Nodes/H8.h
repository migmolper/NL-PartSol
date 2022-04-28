/*! \file H8.h
  \brief Linear quadrilateral shape function

*/

/***********************************************/
/******* 3D cuadrilateral linear element *******/
/***********************************************/

/*     (7)      (6)   */
/*       o-------o    */
/*      /.      /|    */
/* (4) / . (5) / |    */
/*    o--o----o..o    */
/*    | / (3) | / (2) */
/*    |/      |/      */
/*    o-------o       */
/*   (0)     (1)      */


#ifndef _H8_H_
#define _H8_H_


// Shape functions auxilar tools
#include "Nodes/Nodes-Tools.h"


/****************************************************************************/

void   initialize__H8__(Particle, Mesh);
Matrix N__H8__(Matrix);
Matrix dN_Ref__H8__(Matrix);
Matrix dN__H8__(Matrix,Matrix);
bool   in_out__H8__(Matrix, Matrix);
void   element_to_particles__H8__(Matrix, Mesh, int);
double min_DeltaX__H8__(ChainPtr, Matrix);
double volume__H8__(Matrix);
void   local_search__H8__(Particle, Mesh);
#endif