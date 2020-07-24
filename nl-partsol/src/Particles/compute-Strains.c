#include "nl-partsol.h"

/*************************************************************/

Tensor rate_inifinitesimal_Strain__Particles__(Matrix Velocity, Matrix Gradient)
{
  int Ndim = NumberDimensions;
  Tensor Rate_Strain = alloc__TensorLib__(2);
  Tensor Velocity_I;
  Tensor Gradient_I;
  Tensor VoG_I;

  int NodesElem = Gradient.N_rows;

  /* Compute rate of strain */
  for(int I = 0 ; I<NodesElem ; I++){
    /* Assign from matrix to tensor */
    Velocity_I = memory_to_tensor__TensorLib__(Velocity.nM[I], 1);
    Gradient_I = memory_to_tensor__TensorLib__(Gradient.nM[I], 1);
   
    /* Compute the dyadic product of the nodal velocity and the
       gradient of the shape functions */
    VoG_I = dyadic_Product__TensorLib__(Velocity_I, Gradient_I);
    
    /* Ad the nodal contribution to the train tensor */
    for(int i = 0 ; i<Ndim ; i++){
      for(int j = 0 ; j<Ndim ; j++){
	Rate_Strain.N[i][j] +=
	  0.5*(VoG_I.N[i][j] + VoG_I.N[j][i]);
      }
    }
    /* Free memory */
    free__TensorLib__(VoG_I);
  }
  
  return Rate_Strain;
}

/*******************************************************/

Tensor infinitesimal_Strain__Particles__(Tensor Strain, Tensor Rate_Strain, double TimeStep)
{
  int Ndim = NumberDimensions;
  /* Check in the input its is ok */
  if ((Strain.Order == 2) && (Rate_Strain.Order == 2)){
    /* Update strain tensor with the rate of strain tensor */
    for(int i = 0 ; i<Ndim ; i++){
      for(int j = 0 ; j<Ndim ; j++){
	Strain.N[i][j] += TimeStep*Rate_Strain.N[i][j];
      }
    }
  }
  else{
    fprintf(stderr,"%s : %s %s !!! \n",
	    "Error in infinitesimal_Strain__Particles__()",
	    "The input should be",
	    "two tensors of 2nd order and a scalar");
    exit(EXIT_FAILURE);    
  }

  return Strain;
}

/*******************************************************/

Tensor increment_Deformation_Gradient__Particles__(Matrix DeltaU, Matrix gradient_p)
{

  /* Variable definition */
  int Ndim = NumberDimensions;
  int Nnodes_p = DeltaU.N_rows;  
  Tensor f_n1;
  Tensor DeltaU_I;
  Tensor gradient_I;
  Tensor gradient_DeltaU_I;
  
  /*
    Compute increment of the deformation gradient 
    f_n1 = I + Delta_u 0 gradient_N
  */
  
  /* Add identity tensor */
  f_n1 = Identity__TensorLib__();
  
  for(int I = 0 ; I<Nnodes_p ; I++)
    {

      /* Assign from matrix to tensor */
      DeltaU_I = memory_to_tensor__TensorLib__(DeltaU.nM[I], 1);
      gradient_I = memory_to_tensor__TensorLib__(gradient_p.nM[I], 1);
      
      /* Compute the dyadic product of the nodal velocity and the
	 gradient of the shape functions */
      gradient_DeltaU_I = dyadic_Product__TensorLib__(DeltaU_I, gradient_I);
      
      /* Ad the nodal contribution to the train tensor */
      for(int i = 0 ; i<Ndim ; i++)
	{
	  for(int j = 0 ; j<Ndim ; j++)
	    {
	      f_n1.N[i][j] += gradient_DeltaU_I.N[i][j];
	    }
	}
      
      /* Free memory */
      free__TensorLib__(gradient_DeltaU_I);
    }
  
  return f_n1;
}

/*******************************************************/

void update_Deformation_Gradient_n1__Particles__(Tensor F_n1, Tensor F_n, Tensor f_n1)
{
  int Ndim = NumberDimensions;
  double aux;
    
  for(int i = 0 ; i < Ndim  ; i++)
    {
      for(int j = 0 ; j < Ndim  ; j++)
	{
	  /*
	    Set to zero the deformation gradient at t = n + 1 
	   */
	  F_n1.N[i][j] = 0;

	  /*
	    Compute row-column multiplication
	   */
	  aux = 0;
	  for(int k = 0 ; k < Ndim  ; k++)
	    {
	      aux += f_n1.N[i][k]*F_n.N[k][j];
	    }

	  /*
	    New value
	   */
	  F_n1.N[i][j] = aux;
	}
    }
}


/*******************************************************/

Tensor right_Cauchy_Green__Particles__(Tensor F)
{
  /* Define output */
  Tensor C = alloc__TensorLib__(2);
  /* Define the number of dimensions */
  int Ndim = NumberDimensions;

  /* Compute C = F^T F */
  for(int i = 0 ; i < Ndim  ; i++)
    {
      for(int j = 0 ; j < Ndim  ; j++)
	{
	  for(int k = 0 ; k < Ndim  ; k++)
	    {
	      C.N[i][j] += F.N[k][i]*F.N[k][j];
	    }
	}
    }

  return C;
}


/*******************************************************/

Tensor strain_Green_Lagrange__Particles__(Tensor C)
{
  /* Define output */
  Tensor E = alloc__TensorLib__(2);
  /* Define eye tensor */
  Tensor I = Identity__TensorLib__();
  /* Define the number of dimensions */
  int Ndim = NumberDimensions;

  /* Compute E = 1/2 * [ C - I]  */
  for(int i = 0 ; i < Ndim  ; i++)
    {
      for(int j = 0 ; j < Ndim  ; j++)
	{
	  E.N[i][j] = 0.5 * (C.N[i][j] - I.N[i][j]);
	}
    }

  free__TensorLib__(I);

  return E;
}

/*******************************************************/