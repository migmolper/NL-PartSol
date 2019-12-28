#include <stdio.h>
#include <stdlib.h>

#ifndef TypeDefinitions
#define TypeDefinitions
#endif

/*******************************************************/

typedef struct {
  
  /* Constitutive model */
  char Type [100];
  /* Elastic modulus */
  double E;
  /* Poisson ratio */
  double mu;
  /* Activate fracture modulus */
  bool Fracture;
  /* Normalizing constant (Fracture) */
  double Ceps;
  /* Failure energy (Fracture)  */
  double Gf;
  
} Material;

/*******************************************************/

/* Material Linear-Elastic */
Matrix LinearElastic(Matrix,Matrix,
		     double,double);
double W_LinearElastic(Matrix,Matrix,double);

/* Fracture */
Matrix EigenerosionAlgorithm(Matrix, Matrix, Matrix,
			     int *, Material *,
			     ChainPtr *, double);
Matrix ComputeDamage(Matrix, Matrix, Matrix,
		     int *, Material *,
		     ChainPtr *, double);

/*******************************************************/

/* Library */

typedef struct { /* Constitutive models */

  /* Linear elastic material */
  Matrix (* LE)(Matrix,Matrix,double,double);
  
} ConstLib;

ConstLib Contitutive(void);
