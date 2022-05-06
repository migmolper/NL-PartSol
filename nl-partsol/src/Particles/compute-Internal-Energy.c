#include <string.h>
#include "nl-partsol.h"

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
