#include <stdio.h>
#include <stdlib.h>

#ifndef TypeDefinitions
#define TypeDefinitions
#endif

/*******************************************************/

typedef struct {

  /* Elastic modulus */
  double E;
  /* Poisson ratio */
  double mu;
  
} Material;

/*******************************************************/

/* Material Linear elastic */
Matrix LinearElastic(Matrix,Matrix,double,double);
double W_LinearElastic(Matrix,Matrix);

/*******************************************************/

/* Library */

typedef struct { /* Constitutive models */

  /* Linear elastic material */
  Matrix (* LE)(Matrix,Matrix,double,double);
  
} ConstLib;

ConstLib Contitutive(void);
