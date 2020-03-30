#include <stdio.h>
#include <stdlib.h>

#ifndef TypeDefinitions
#define TypeDefinitions
#endif

/*******************************************************/

typedef struct {

  /* Name and id of the material */
  int Id;
  char Type [100];
  /* Initial density */
  double rho;
  /* Elastic modulus */
  double E;
  /* Poisson ratio */
  double mu;
  /* Activate fracture modulus */
  bool Eigenerosion;
  bool Eigensoftening;
  /* Normalizing constant (Eigenerosion) */
  double Ceps;
  /* Failure energy (Eigenerosion)  */
  double Gf;
  /* Tensile strengt of the material (Eigensoftening) */
  double ft;
  /* Bandwidth of the cohesive fracture (Eigensoftening) */
  double heps;
  /* Critical opening displacement */
  double Wc;
  
} Material;

/*******************************************************/

/* Material Linear-Elastic */
Matrix LinearElastic(Matrix, Matrix, Material);
double W_LinearElastic(Matrix, Matrix, double);

/* Fracture */
void EigenerosionAlgorithm(Matrix, Matrix,
			   Matrix, Matrix,
			   int *, Material *,
			   ChainPtr *, double);
void EigensofteningAlgorithm(Matrix, Matrix,
			     Matrix, Matrix,
			     int *, Material *,
			     ChainPtr *, double);

/*******************************************************/
