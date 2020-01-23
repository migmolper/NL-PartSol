#include <stdio.h>
#include <stdlib.h>

#ifndef TypeDefinitions
#define TypeDefinitions
#endif

/*******************************************************/

typedef struct {

  /* Linear elastic material */
  Matrix (* LE)(Matrix,Matrix,double,double);
  
} ConstLib;

typedef struct {

  /* Constitutive model */
  ConstLib D; /* Constitutive response */
  char Type [100];
  /* Initial density */
  double rho;
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
Matrix EigenerosionAlgorithm(Matrix, Matrix,
			     Matrix, Matrix,
			     int *, Material *,
			     ChainPtr *, double);
Matrix ComputeDamage(Matrix, Matrix,
		     Matrix, Matrix,
		     int *, Material *,
		     ChainPtr *, double);

/*******************************************************/

/* Library */

ConstLib Contitutive(char *);
