#include "nl-partsol.h"

Tensor LinearElastic(Tensor Strain, Tensor Stress, Material Mat)
{

  int Ndim = NumberDimensions;
  
  /*! 
    Check in the input its is ok 
  */
  if ((Strain.Order == 2) && (Stress.Order == 2))
    {
      /*!
	Define material and other properties 
      */
      double nu = Mat.nu; 
      double E = Mat.E;
      double lambda = nu*E/((1-nu*2)*(1+nu));
      double G = E/(2*(1+nu));
      double traceStrain = I1__TensorLib__(Strain);
      Tensor I = Identity__TensorLib__();

      for(int i = 0 ; i<Ndim ; i++)
	{
	  /*!
	    Add first component 
	  */
	  for(int j = 0 ; j<Ndim ; j++)
	    {
	      Stress.N[i][j] = 2*G*Strain.N[i][j];
	    }
	  /*!
	    Add diagonal component 
	  */
	  Stress.N[i][i] += lambda*traceStrain*I.N[i][i];
	}

      free__TensorLib__(I);
    }
  else
    {
      fprintf(stderr,"%s : %s !!! \n",
	      "Error in LinearElastic()",
	      "The input should be 2nd tensor and a 2nd tensor");
      exit(EXIT_FAILURE);
    }
  return Stress;  
}

/**************************************************************/

Tensor volumetric_stress__LinearElastic__(double E_elastic_vol, Material MatProp)
{

  int Ndim = NumberDimensions;
  double nu = MatProp.nu; 
  double E = MatProp.E;
  double K = E/(3*(1-2*nu));
  double aux = K*E_elastic_vol;
  Tensor p_trial = alloc__TensorLib__(2);

  /*
    Compute deviatoric stress tensor
  */
  for(int i = 0 ; i<Ndim ; i++)
  {
    p_trial.N[i][i] = aux;
  }

  return p_trial;
}

/**************************************************************/

Tensor deviatoric_stress__LinearElastic__(Tensor E_elastic_dev, Material MatProp)
{

  int Ndim = NumberDimensions;
  double nu = MatProp.nu; 
  double E = MatProp.E;
  double G = E/(2*(1+nu));
  Tensor s_trial = alloc__TensorLib__(2);
  /*
    Compute deviatoric stress tensor
  */

  for(int i = 0 ; i<Ndim ; i++)
    {
      for(int j = 0 ; j<Ndim ; j++)
      {
        s_trial.N[i][j] = 2*G*E_elastic_dev.N[i][j];
      }
    }

  return s_trial;
}

/**************************************************************/
