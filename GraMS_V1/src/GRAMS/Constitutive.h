#include <stdio.h>
#include <stdlib.h>

#ifndef TypeDefinitions
#define TypeDefinitions
#endif

/*******************************************************/

typedef struct {

  /* Density field */
  Matrix rho;
  /* Mass field */
  Matrix mass;
  /* Elastic modulus */
  Matrix E;
  /* Poisson ratio */
  Matrix mu;
  /* Normalizing constant (Fracture) */
  Matrix Ceps;
  /* Failure energy (Fracture)  */
  Matrix Gf;
  /* Damage parameter (Fracture) */
  Matrix ji;
  
} Material;

/*******************************************************/

/* Material Linear-Elastic */
Matrix LinearElastic(Matrix,Matrix,
		     double,double);
double W_LinearElastic(Matrix,Matrix,double);

/* Fracture */
Matrix EigenerosionAlgorithm(Matrix, Matrix,
			     Matrix, Matrix,
			     Matrix, ChainPtr *);
Matrix ComputeDamage(Matrix, Material, ChainPtr *);

/*******************************************************/

/* Library */

typedef struct { /* Constitutive models */

  /* Linear elastic material */
  Matrix (* LE)(Matrix,Matrix,double,double);
  
} ConstLib;

ConstLib Contitutive(void);
