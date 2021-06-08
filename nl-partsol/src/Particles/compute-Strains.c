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
  Tensor F_bar,
  Particle MPM_Mesh,
  Mesh FEM_Mesh)
{

  int q;
  int Ndim = NumberDimensions;
  int I0_p = MPM_Mesh.I0[p];
  int IdxElement;
  double Simplex_Radius = FEM_Mesh.DeltaX;
  double J_n_p_patch;
  double J_n1_q_patch;
  double Vol_0_q;
  double Vol_n1_q;
  double V0_patch;
  double Vn1_patch;
  double J_p;
  double J_n1_patch;
  double J_averaged;
  double averaged_F_vol;

  double alpha = 0.5;

  Matrix X_p;
  Matrix X_q;
  Matrix Coordinates_Patch_p;
  ChainPtr Elements_Near_I0;
  ChainPtr Node_Patch_p;
  ChainPtr Particle_Patch_p;

  Tensor F_n1_p;

  V0_patch = 0.0;
  Vn1_patch = 0.0;

  /* Get the coordinates of the particle */
  X_p = memory_to_matrix__MatrixLib__(Ndim,1,MPM_Mesh.Phi.x_GC.nM[p]);

  /* List of elements near the particle */
  Elements_Near_I0 = FEM_Mesh.NodeNeighbour[I0_p];

  /* Get the element inside of the */
  IdxElement = search_particle_in_surrounding_elements__Particles__(p,X_p,Elements_Near_I0,FEM_Mesh);

  if(IdxElement != -999)
  {
    // Get the surrounding nodes  
    Node_Patch_p = FEM_Mesh.Connectivity[IdxElement];

    Coordinates_Patch_p = get_nodes_coordinates__MeshTools__(Node_Patch_p, FEM_Mesh.Coordinates);

    // Loop in the sourrounding particles to get the deformed and reference volumes
    while(Node_Patch_p != NULL)
    {
      Particle_Patch_p = NULL;

      // Get the list of particles close to this node
      Particle_Patch_p = FEM_Mesh.I_particles[Node_Patch_p->I];

      while(Particle_Patch_p != NULL)
      {
        q = Particle_Patch_p->I;
        X_q = memory_to_matrix__MatrixLib__(Ndim,1,MPM_Mesh.Phi.x_GC.nM[q]);

        if(FEM_Mesh.In_Out_Element(X_q,Coordinates_Patch_p))
        {
          Vol_0_q = MPM_Mesh.Phi.Vol_0.nV[q];

          J_n1_q_patch = MPM_Mesh.Phi.J.nV[q];
          Vol_n1_q = Vol_0_q*J_n1_q_patch;

          V0_patch += Vol_0_q;
          Vn1_patch += Vol_n1_q;
        }

        Particle_Patch_p = Particle_Patch_p->next;
      }

      Node_Patch_p = Node_Patch_p->next; 

    }


    free__MatrixLib__(Coordinates_Patch_p);

  }

  // Compute the averaged jacobian of the deformation gradient
  J_p = MPM_Mesh.Phi.J.nV[p];
  J_n1_patch = Vn1_patch/V0_patch;
  J_averaged = J_n1_patch/J_p;

  // Compute the averaged volume of the deformation gradient
  averaged_F_vol = pow(J_averaged,(double)1/Ndim);

  F_n1_p = memory_to_tensor__TensorLib__(MPM_Mesh.Phi.F_n1.nM[p],2);

  // Update the deformation gradient to avoid locking (F-bar)
  for(int i = 0 ; i<Ndim ; i++)
  {
    for(int j = 0 ; j<Ndim ; j++)
    {
      F_bar.N[i][j] = alpha*F_n1_p.N[i][j] + (1 - alpha)*averaged_F_vol*F_n1_p.N[i][j];
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

Tensor logarithmic_strains__Particles__(
  Tensor C)
/*
    Use the approach of Ortiz and Camacho to compute the elastic infinitesimal strain tensor.
*/
{

  int Ndim = NumberDimensions;
  EigenTensor Eigen_C;
  Tensor logC_spectral = alloc__TensorLib__(2);
  Tensor logC;


  /*
    Compute the spectral descomposition of the tensor logC
  */
  Eigen_C = Eigen_analysis__TensorLib__(C);
  
  for(int i = 0 ; i < Ndim  ; i++)
  {
    logC_spectral.N[i][i] = 0.5*log(Eigen_C.Value.n[i]);
  }

  /*
    Rotate the spectral descomposition of the tensor logC
  */
  logC = rotate__TensorLib__(logC_spectral, Eigen_C.Vector);


  /*
    Free memory
  */
  free__TensorLib__(Eigen_C.Value);
  free__TensorLib__(Eigen_C.Vector);
  free__TensorLib__(logC_spectral);

  return logC; 
}

/*******************************************************/

Tensor increment_Deformation_Gradient_exponential_strains__Particles__(
  Tensor D_E)
{

  int Ndim = NumberDimensions;
  EigenTensor Eigen_D_E;
  Tensor exp_D_E_spectral = alloc__TensorLib__(2);
  Tensor exp_D_E;


  /*
    Compute the spectral descomposition of the tensor exp(D_E)
  */
  Eigen_D_E = Eigen_analysis__TensorLib__(D_E);

  
  for(int i = 0 ; i < Ndim  ; i++)
  {
    exp_D_E_spectral.N[i][i] = exp(Eigen_D_E.Value.n[i]);
  }

  /*
    Rotate the spectral descomposition of the tensor exp(D_E)
  */
  exp_D_E = rotate__TensorLib__(exp_D_E_spectral, Eigen_D_E.Vector);

  /*
    Free memory
  */
  free__TensorLib__(Eigen_D_E.Value);
  free__TensorLib__(Eigen_D_E.Vector);
  free__TensorLib__(exp_D_E_spectral);

  return exp_D_E; 
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

/*******************************************************/

void update_plastic_deformation_gradient__Particles__(
  Tensor Increment_E_plastic, 
  Tensor F_m1_plastic)
{
  int Ndim = NumberDimensions;
  EigenTensor Eigen_Increment_E_plastic;
  Tensor D_F_elastic_spectral = alloc__TensorLib__(2);
  Tensor D_F_elastic;
  Tensor F_m1_plastic_new;

  /*
    Compute the spectral descomposition of the tensor exp(D_E)
  */
  Eigen_Increment_E_plastic = Eigen_analysis__TensorLib__(Increment_E_plastic);

  
  for(int i = 0 ; i < Ndim  ; i++)
  {
    D_F_elastic_spectral.N[i][i] = exp(- Eigen_Increment_E_plastic.Value.n[i]);
  }

  /*
    Rotate the spectral descomposition of the tensor exp(D_E)
  */
  D_F_elastic = rotate__TensorLib__(D_F_elastic_spectral, Eigen_Increment_E_plastic.Vector);

  /*
    Compute the new value of the plastic deformation gradient 
  */
  F_m1_plastic_new = matrix_product__TensorLib__(F_m1_plastic,D_F_elastic);

  /*
    Update value of the plastic deformation gradient
  */
  for(int i = 0 ; i < Ndim  ; i++)
  {
    for(int j = 0 ; j < Ndim  ; j++)
    {
      F_m1_plastic.N[i][j] = F_m1_plastic_new.N[i][j];
    }
  }

  /*
    Free memory
  */
  free__TensorLib__(Eigen_Increment_E_plastic.Value);
  free__TensorLib__(Eigen_Increment_E_plastic.Vector);
  free__TensorLib__(D_F_elastic_spectral);
  free__TensorLib__(D_F_elastic);
  free__TensorLib__(F_m1_plastic_new);

}

/**************************************************************/

void update_elastic_deformation_gradient__Particles__(
  Tensor F_elastic,
  Tensor F_trial_elastic,
  Tensor Increment_E_plastic)
{
  int Ndim = NumberDimensions;
  EigenTensor Eigen_Increment_E_plastic;
  Tensor D_F_elastic_spectral = alloc__TensorLib__(2);
  Tensor D_F_elastic;
  Tensor F_elastic_new;

  /*
    Compute the spectral descomposition of the tensor exp(D_E)
  */
  Eigen_Increment_E_plastic = Eigen_analysis__TensorLib__(Increment_E_plastic);

  
  for(int i = 0 ; i < Ndim  ; i++)
  {
    D_F_elastic_spectral.N[i][i] = exp(-Eigen_Increment_E_plastic.Value.n[i]);
  }

  /*
    Rotate the spectral descomposition of the tensor exp(D_E)
  */
  D_F_elastic = rotate__TensorLib__(D_F_elastic_spectral, Eigen_Increment_E_plastic.Vector);


  /*
    Compute the new value of the plastic deformation gradient 
  */
  F_elastic_new = matrix_product__TensorLib__(F_trial_elastic,D_F_elastic);

  /*
    Update value of the plastic deformation gradient
  */
  for(int i = 0 ; i < Ndim  ; i++)
  {
    for(int j = 0 ; j < Ndim  ; j++)
    {
      F_elastic.N[i][j] = F_elastic_new.N[i][j];
    }
  }

  /*
    Free memory 
  */
  free__TensorLib__(Eigen_Increment_E_plastic.Value);
  free__TensorLib__(Eigen_Increment_E_plastic.Vector);
  free__TensorLib__(D_F_elastic_spectral);
  free__TensorLib__(D_F_elastic);
  free__TensorLib__(F_elastic_new);
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
