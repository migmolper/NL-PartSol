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
  /* Material celerity */
  double Cel;
  /* Initial density */
  double rho;
  /* Elastic modulus */
  double E;
  /* Poisson ratio */
  double mu;
  /* Thickness of the Material */
  double thickness;
  /* Activate fracture modulus */
  bool Eigenerosion;
  bool Eigensoftening;
  /*Normalizing constant (Eigenerosion/Eigensoftening) */
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

/* Solid rigid */
Tensor SolidRigid(Tensor);

/* Material Linear-Elastic */
Tensor LinearElastic(Tensor, Tensor, Material);


/* Fracture */
void EigenerosionAlgorithm(Matrix, Matrix,
			   Matrix, Matrix,
			   int *, Material *,
			   ChainPtr *, double);
void EigensofteningAlgorithm(Matrix, Matrix,
			     Matrix, Matrix,
			     Matrix, int *,
			     Material *, ChainPtr *);

/*******************************************************/
