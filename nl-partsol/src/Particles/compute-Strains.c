#include "nl-partsol.h"

/*************************************************************/

Tensor rate_inifinitesimal_Strain__Particles__(
  Matrix Velocity, 
  Matrix Gradient)
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
    for(int i = 0 ; i<Ndim ; i++)
    {
      for(int j = 0 ; j<Ndim ; j++)
      {
        Rate_Strain.N[i][j] += 0.5*(VoG_I.N[i][j] + VoG_I.N[j][i]);
      }
    }
    /* Free memory */
    free__TensorLib__(VoG_I);
  }
  
  return Rate_Strain;
}

/*******************************************************/

Tensor infinitesimal_Strain__Particles__(
  Tensor Strain, 
  Tensor Rate_Strain, 
  double TimeStep)
{
  int Ndim = NumberDimensions;
  /* Check in the input its is ok */
  if ((Strain.Order == 2) && (Rate_Strain.Order == 2))
  {
    /* Update strain tensor with the rate of strain tensor */
    for(int i = 0 ; i<Ndim ; i++)
    {
      for(int j = 0 ; j<Ndim ; j++)
      {
        Strain.N[i][j] += TimeStep*Rate_Strain.N[i][j];
      }
    }
  }
  else
  {
    fprintf(stderr,"%s : %s %s !!! \n",
	    "Error in infinitesimal_Strain__Particles__()",
	    "The input should be",
	    "two tensors of 2nd order and a scalar");
    exit(EXIT_FAILURE);    
  }

  return Strain;
}

/*******************************************************/

void update_increment_Deformation_Gradient__Particles__(
  Tensor DF_p, 
  Matrix DeltaU, 
  Matrix gradient_p)
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
  
  /* Initialise with the identity tensor */
  for(int i = 0 ; i<Ndim ; i++)
  {
    for(int j = 0 ; j<Ndim ; j++)
    {
      DF_p.N[i][j] = 1*(i==j);
    }
  }
  
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
        DF_p.N[i][j] += gradient_DeltaU_I.N[i][j];
      }
    }
      
    /* Free memory */
    free__TensorLib__(gradient_DeltaU_I);
  }
  
}

/*******************************************************/

void update_rate_increment_Deformation_Gradient__Particles__(
  Tensor dt_DF_p,
  Matrix DeltaV,
  Matrix gradient_p)
{

  /* Variable definition */
  int Ndim = NumberDimensions;
  int Nnodes_p = DeltaV.N_rows;  
  Tensor f_n1;
  Tensor DeltaV_I;
  Tensor gradient_I;
  Tensor gradient_DeltaV_I;
  
  /*
    Compute increment of the deformation gradient 
    dt_f_n1 = I + (Delta_V o gradient_N)
  */
  
  /* Initialise with the identity tensor */
  for(int i = 0 ; i<Ndim ; i++)
  {
    for(int j = 0 ; j<Ndim ; j++)
    {
      dt_DF_p.N[i][j] = 0.0;
    }
  }
  
  for(int I = 0 ; I<Nnodes_p ; I++)
  {

    /* Assign from matrix to tensor */
    DeltaV_I = memory_to_tensor__TensorLib__(DeltaV.nM[I], 1);
    gradient_I = memory_to_tensor__TensorLib__(gradient_p.nM[I], 1);
      
    /* Compute the dyadic product of the nodal velocity and the
        gradient of the shape functions */
    gradient_DeltaV_I = dyadic_Product__TensorLib__(DeltaV_I, gradient_I);
      
    /* Ad the nodal contribution to the train tensor */
    for(int i = 0 ; i<Ndim ; i++)
    {
      for(int j = 0 ; j<Ndim ; j++)
      {
        dt_DF_p.N[i][j] += gradient_DeltaV_I.N[i][j];
      }
    }
      
    /* Free memory */
    free__TensorLib__(gradient_DeltaV_I);
  }
} 

/*******************************************************/

void update_Deformation_Gradient_n1__Particles__(
  Tensor F_n1,
  Tensor F_n,
  Tensor f_n1)
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

void get_locking_free_Deformation_Gradient_n1__Particles__(
  int p,
  Particle MPM_Mesh)
{

  int Ndim = NumberDimensions;
  int MatIndx_p = MPM_Mesh.MatIdx[p];

  double J_n1_patch = MPM_Mesh.Phi.Jbar.nV[p];
  double J_p;
  double J_averaged;
  double averaged_F_vol;

  double alpha = MPM_Mesh.Mat[MatIndx_p].alpha_Fbar;

  Tensor F_n1_p = memory_to_tensor__TensorLib__(MPM_Mesh.Phi.F_n1.nM[p],2);
  Tensor Fbar   = memory_to_tensor__TensorLib__(MPM_Mesh.Phi.Fbar.nM[p],2);

  // Compute the averaged jacobian of the deformation gradient
  J_p = MPM_Mesh.Phi.J.nV[p];
  J_averaged = J_n1_patch/J_p;

  // Compute the averaged volume of the deformation gradient
  averaged_F_vol = pow(J_averaged,(double)1/Ndim);

  // Update the deformation gradient to avoid locking (F-bar)
  for(int i = 0 ; i<Ndim ; i++)
  {
    for(int j = 0 ; j<Ndim ; j++)
    {
      Fbar.N[i][j] = alpha*F_n1_p.N[i][j] + (1 - alpha)*averaged_F_vol*F_n1_p.N[i][j];
    }
  }

}

/*******************************************************/

void update_rate_Deformation_Gradient_n1__Particles__(
  Tensor dt_F_n1, 
  Tensor dt_f_n1, 
  Tensor F_n, 
  Tensor f_n1, 
  Tensor dt_F_n)
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
      dt_F_n1.N[i][j] = 0;

      /*
        Compute row-column multiplication
      */
      aux = 0;
      for(int k = 0 ; k < Ndim  ; k++)
      {
        aux += dt_f_n1.N[i][k]*F_n.N[k][j] + f_n1.N[i][k]*dt_F_n.N[k][j];
      }

      /*
        New value
      */
      dt_F_n1.N[i][j] = aux;
    }
  }
}

/*******************************************************/

double compute_Jacobian_Rate__Particles__(
  double J_p,
  Tensor F_p,
  Tensor dt_F_p)
{
  /* Variable definition */
  double dt_J_p = 0;
  Tensor inverse_F_p = Inverse__TensorLib__(F_p);
  Tensor transpose_inverse_F_p = transpose__TensorLib__(inverse_F_p);

  dt_J_p = J_p*inner_product__TensorLib__(transpose_inverse_F_p, dt_F_p);

  free__TensorLib__(inverse_F_p);
  free__TensorLib__(transpose_inverse_F_p);

  return dt_J_p;
}

/*******************************************************/

Tensor right_Cauchy_Green__Particles__(
  Tensor F)
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

Tensor strain_Green_Lagrange__Particles__(
  Tensor C)
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


/**************************************************************/

Matrix compute_B_matrix__Particles__(
  Tensor F, 
  Tensor GRAD_N)
/*
B=[F(1,1)*gradNptog(1,nn1) F(2,1)*gradNptog(1,nn1);...
   F(1,2)*gradNptog(2,nn1) F(2,2)*gradNptog(2,nn1);...
   (F(1,1)*gradNptog(2,nn1)+F(1,2)*gradNptog(1,nn1)) (F(2,1)*gradNptog(2,nn1)+F(2,2)*gradNptog(1,nn1))];
*/
{
  Matrix B = allocZ__MatrixLib__(3,2);

  B.nM[0][0] = F.N[0][0]*GRAD_N.n[0]; 
  B.nM[0][1] = F.N[1][0]*GRAD_N.n[0];
  B.nM[1][0] = F.N[0][1]*GRAD_N.n[1]; 
  B.nM[1][1] = F.N[1][1]*GRAD_N.n[1];
  B.nM[2][0] = F.N[0][0]*GRAD_N.n[1] + F.N[0][1]*GRAD_N.n[0]; 
  B.nM[2][1] = F.N[1][0]*GRAD_N.n[1] + F.N[1][1]*GRAD_N.n[0];

  return B;
}

/**************************************************************/

Matrix compute_BT_matrix__Particles__(
  Tensor F, 
  Tensor GRAD_N)
{
  Matrix BT = allocZ__MatrixLib__(2,3);

  BT.nM[0][0] = F.N[0][0]*GRAD_N.n[0];
  BT.nM[1][0] = F.N[1][0]*GRAD_N.n[0];
  BT.nM[0][1] = F.N[0][1]*GRAD_N.n[1];
  BT.nM[1][1] = F.N[1][1]*GRAD_N.n[1];
  BT.nM[0][2] = F.N[0][0]*GRAD_N.n[1] + F.N[0][1]*GRAD_N.n[0];
  BT.nM[1][2] = F.N[1][0]*GRAD_N.n[1] + F.N[1][1]*GRAD_N.n[0];

  return BT;
}

/**************************************************************/



/**************************************************************/
