/*! \file Q4.h
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


#ifndef _Q8_H_
#define _Q8_H_

void   initialize__Q8__(GaussPoint, Mesh);
Matrix N__Q8__(Matrix);
Matrix dN_Ref__Q8__(Matrix);
Matrix dN__Q8__(Matrix,Matrix);
void   X_to_Xi__Q8__(Matrix,Matrix,Matrix);
void   element_to_particles__Q8__(Matrix, Mesh, int);

#endif