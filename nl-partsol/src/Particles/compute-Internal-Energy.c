#include <string.h>
#include "nl-partsol.h"

double internal_energy__Particles__(Tensor Strain, Tensor Stress) {

  /* Internal energy for the Gauss-Point */
  double W = 0;
  /*Check in the input its is ok */
  if ((Strain.Order == 2) && (Stress.Order == 2)) {
    /* Calcule the internal work */
    W = 0.5 * inner_product__TensorLib__(Strain, Stress);
  } else {
    fprintf(stderr, "%s : %s !!! \n", "Error in internal_energy__Particles__()",
            "The input should be 2nd tensor and a 2nd tensor");
    exit(EXIT_FAILURE);
  }
  return W;
}

/*********************************************************************/

void finite_strains_internal_energy__Particles__(
  double * W,
  const double * P,
  const double * dFdt,
  double dt) {

  unsigned Ndim = NumberDimensions;
  double P__x__dFdt = 0.0;
  
  for (unsigned i = 0; i < Ndim*Ndim; i++)
  {
    P__x__dFdt += P[i] * dFdt[i];
  }
  
  *W += dt*P__x__dFdt;
}

/*********************************************************************/
