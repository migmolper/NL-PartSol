/***********************************************/
/********* 3D triangle linear element **********/
/***********************************************/

/*        o (3)    */
/*       /|\       */
/*      / | \      */
/* (0) /  |  \     */
/*    o-------o    */
/*     \  |  / (2) */
/*      \ | /      */ 
/*       \|/       */
/*    (1) o        */

#ifndef _T4_H_
#define _T4_H_

void   initialize__T4__(GaussPoint, Mesh);
Matrix N__T4__(Matrix);
Matrix dN__T4__(Matrix,Matrix);
Matrix dN_Ref__T4__(Matrix);
void   X_to_Xi__T4__(Matrix,Matrix,Matrix);
void   element_to_particles__T4__(Matrix, Mesh, int);

#endif