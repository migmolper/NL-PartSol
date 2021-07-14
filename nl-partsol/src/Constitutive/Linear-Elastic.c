#include "nl-partsol.h"

/*************************************************************/

State_Parameters compute_kirchhoff_isotropic_linear_elasticity(
  State_Parameters Intput_SP,
  Material MatProp_p)
{

  int Ndim = NumberDimensions;

  /* Get information from the state parameter */
  Tensor Strain = memory_to_tensor__TensorLib__(Intput_SP.Strain,2);
  Tensor Stress = memory_to_tensor__TensorLib__(Intput_SP.Stress,2);

  double nu = MatProp_p.nu; 
  double E = MatProp_p.E;
  double K = E/(3*(1-2*nu));
  double G = E/(2*(1+nu));
  double traceStrain = I1__TensorLib__(Strain);

  for(int i = 0 ; i<Ndim ; i++)
  {
    for(int j = 0 ; j<Ndim ; j++)
    {
      Stress.N[i][j] = 2*G*Strain.N[i][j] + (K - 2*G/3.0)*(i==j)*traceStrain;
    }
  }

  /* Plain stress condition */
  if(Ndim == 2)
  {
    Intput_SP.Stress[4] = (K - 2*G/3.0)*traceStrain;
  }

  /* 
    Define output
  */
  State_Parameters Output_SP;

  return Output_SP;  
}

/**************************************************************/

