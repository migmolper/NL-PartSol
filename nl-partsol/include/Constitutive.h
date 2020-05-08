#ifndef _CONSTITUTIVE_H_
#define _CONSTITUTIVE_H_

/*******************************************************/

/* Solid rigid */
Tensor SolidRigid(Tensor);

/* Material Linear-Elastic */
Tensor LinearElastic(Tensor, Tensor, Material);


/* Fracture */
void EigenerosionAlgorithm(int, Matrix, Matrix,
			   Matrix, Matrix, Matrix,
			   Material, ChainPtr *, double);
void EigensofteningAlgorithm(int, Matrix, Matrix,
			     Matrix, Matrix, Matrix,
			     Material, ChainPtr *);

/*******************************************************/


#endif
