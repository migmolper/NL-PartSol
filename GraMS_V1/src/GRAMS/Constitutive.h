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

Matrix LinearElastic2D(Matrix,double,double);

/*******************************************************/

/* Library */

typedef struct { /* Constitutive models */

  /* 2D linear elastic material */
  Matrix (* LE2D)(Matrix,double,double);
  
} ConstLib;

ConstLib Contitutive(void);
