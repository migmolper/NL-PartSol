#include "nl-partsol.h"

#ifdef __linux__
#include <lapacke.h>

#elif __APPLE__
#include <Accelerate/Accelerate.h>

#endif

/*
  Auxiliar functions 
*/
static void update_StressField_SV(Matrix,Matrix,Matrix,GaussPoint,Mesh,double);
static void compute_NodalMagnitudes_SV(Matrix,Matrix,Matrix,GaussPoint,Mesh);
static void update_VelocityField_SV(Matrix,Matrix,Matrix,GaussPoint,Mesh,double);
static void update_Particles_SV(Matrix,Matrix,GaussPoint,Mesh,double);


/**************************************************************/

void SV_Two_Steps_Taylor_Galerkin(Mesh FEM_Mesh, GaussPoint MPM_Mesh, int InitialStep)
{

  /* Some auxiliar variables for the outputs */
  Matrix List_Fields;
  /* Time step */
  int TimeStep;

  int Ndim = NumberDimensions;
  int Nnodes = FEM_Mesh.NumNodesMesh;

  Matrix V_I, S_I, Div_S_I, M_I;

  for(TimeStep = InitialStep ; TimeStep<NumTimeStep ; TimeStep++ )
    {
      print_Status("*************************************************",TimeStep);
      DeltaTimeStep = DeltaT_CFL(MPM_Mesh, FEM_Mesh.DeltaX);
      print_step(TimeStep,DeltaTimeStep);
    
      print_Status("*************************************************",TimeStep);
      print_Status("First step : Get the nodal fields ... WORKING",TimeStep);

      /* Allocate the output list of fields */
      V_I = MatAllocZ(Nnodes,Ndim);
      strcpy(V_I.Info,"VELOCITY");
      S_I = MatAllocZ(Nnodes,Ndim*Ndim);
      strcpy(S_I.Info,"STRESS");
      M_I = MatAllocZ(Nnodes,1);
      strcpy(M_I.Info,"MASS");
      
      compute_NodalMagnitudes_SV(S_I,V_I,M_I,MPM_Mesh,FEM_Mesh);
      /* imposse_BoudaryConditions_SV(S_I,V_I,FEM_Mesh,TimeStep); */
      print_Status("DONE !!!",TimeStep);

      print_Status("*************************************************",TimeStep);
      print_Status("Second step : Compute equilibrium ... WORKING",TimeStep);
      update_StressField_SV(S_I,V_I,M_I,MPM_Mesh,FEM_Mesh,DeltaTimeStep);
      update_VelocityField_SV(S_I,V_I,M_I,MPM_Mesh,FEM_Mesh,DeltaTimeStep);
      print_Status("DONE !!!",TimeStep);
 
      print_Status("*************************************************",TimeStep);
      print_Status(" Third step : Update lagrangian ... WORKING",TimeStep);
      update_Particles_SV(S_I,V_I,MPM_Mesh,FEM_Mesh,DeltaTimeStep);
      LocalSearchGaussPoints(MPM_Mesh,FEM_Mesh);
      print_Status("DONE !!!",TimeStep);

      if(TimeStep % ResultsTimeStep == 0)
      	{
      	  /* Print Nodal values after appling the BCCs */
      	  WriteVtk_FEM("Mesh",FEM_Mesh,V_I,
      		       (int)TimeStep/ResultsTimeStep);
      	  /* Print GPs results */
      	  WriteVtk_MPM("MPM_VALUES",MPM_Mesh,List_Fields,
      		       (int)TimeStep/ResultsTimeStep);
      	}
    
      print_Status("*************************************************",TimeStep);
      print_Status(" Five step : Reset nodal values ... WORKING",TimeStep);
      FreeMat(S_I);
      FreeMat(V_I);
      FreeMat(M_I);
      print_Status("DONE !!!",TimeStep);

    }

}


/*************************************************************/

static void update_StressField_SV(Matrix S_I,Matrix V_I,Matrix M_I,
         GaussPoint MPM_Mesh,Mesh FEM_Mesh,double Dt)
{ 
  /* Define some dimensions */
  int Ndim = NumberDimensions;
  int Nnodes = FEM_Mesh.NumNodesMesh;
  int Np = MPM_Mesh.NumGP;
  int Nn_p;
 
  /* Auxiliar values for magnitudes of each particle */
  int Ip, Idx_Mat_p;
  double rho_p, m_p, V_p, ShapeFunction_pI;
  Tensor Rate_Strain_p, Rate_Stress_p; 
  Element Nodes_p;
  Matrix Nodal_Velocity_p, ShapeFunction_p, Gradient_p;
  Material Material_p;

  /* Geometrical mass matrix */
  double G_I;

  /* Allocate the RHS for the calculus */
  Matrix RHS_1 = MatAllocZ(Nnodes,Ndim*Ndim);
  Matrix RHS_2 = MatAllocZ(Nnodes,Ndim*Ndim);

  /* Loop in the particles */
  for(int p = 0 ; p<Np ; p++)
    {
      m_p = MPM_Mesh.Phi.mass.nV[p];
      Idx_Mat_p = MPM_Mesh.MatIdx[p];
      Material_p = MPM_Mesh.Mat[Idx_Mat_p];

      /* Define a list of nodes for each particle */
      Nn_p = MPM_Mesh.NumberNodes[p];
      Nodes_p = get_Element(p,MPM_Mesh.ListNodes[p],Nn_p);
    
      /* Compute shape functions */
      ShapeFunction_p = compute_ShapeFunction(Nodes_p, MPM_Mesh,FEM_Mesh);
      /* Compute gradient of the shape function in each node */
      Gradient_p = compute_ShapeFunction_Gradient(Nodes_p, MPM_Mesh,FEM_Mesh);
   
      /* Compute rate of stress tensor */
      Nodal_Velocity_p = get_Element_Field(V_I, Nodes_p);
      Rate_Strain_p = compute_RateOfStrain(Nodal_Velocity_p,Gradient_p);
      Rate_Stress_p = alloc_Tensor(2);
      Rate_Stress_p = compute_Stress(Rate_Strain_p,Rate_Stress_p,Material_p);
        
      /* Compute the RHS with the rate of stress tensor */
      for(int I = 0 ; I<Nn_p ; I++)
  {

    /* Nodal value of the shape function for the particle */
    ShapeFunction_pI = ShapeFunction_p.nV[I];

    /* Get the node of the mesh for the contribution */
    Ip = Nodes_p.Connectivity[I];
      
    /* Upgrade the value of the stress tensor */
    for(int i = 0 ; i<Ndim ; i++)
      {
        for(int j = 0 ; j<Ndim ; j++)
    {
      RHS_1.nM[Ip][i*j] +=
        Dt*ShapeFunction_pI*Rate_Stress_p.N[i][j]*m_p;
      /* RHS_2.nM[Ip][i*j] = 0.5*Dt*Dt* */
    }
      }
  }

      /* Free the matrix with the nodal velocity of the particle */
      free(Nodes_p.Connectivity);
      FreeMat(Nodal_Velocity_p);
      FreeMat(Gradient_p);
      FreeMat(ShapeFunction_p);
      free_Tensor(Rate_Stress_p);
      free_Tensor(Rate_Strain_p);    
    }

  /* Update the nodal value of the stress tensor */
  for(int I = 0 ; I<Nnodes ; I++)
    {
      if(M_I.nV[I] > TOL_zero)
  {
    for(int i = 0 ; i<Ndim*Ndim ; i++)
      {      
        S_I.nM[I][i] += (RHS_1.nM[I][i] + RHS_2.nM[I][i])/M_I.nV[I];
      }
  }    
    }

  /* PrintMatrix(RHS_1,Nnodes,Ndim*Ndim); */

  /* Free the matrix with the RHS */
  FreeMat(RHS_1);
  FreeMat(RHS_2);
}

/*******************************************************/

static void compute_NodalMagnitudes_SV(Matrix S_I,Matrix V_I, Matrix M_I,
        GaussPoint MPM_Mesh,Mesh FEM_Mesh)
/*
  This function performs a information trasference from
  the lagrangian particles to the nodes of the eulerian mesh.
  Output : {Stress_I | V_I | M_I }
*/
{

  /* Define some dimensions */
  int Ndim = NumberDimensions;
  int Nnodes = FEM_Mesh.NumNodesMesh;
  int Np = MPM_Mesh.NumGP;
  int Ip;

  /* Auxiliar variables */
  Matrix ShapeFunction_p;  /* Value of the shape-function */
  double ShapeFunction_pI; /* Evaluation of the GP in the node */
  double M_Ip;
  double m_p; /* Mass of the particle */
  Element Nodes_p; /* Element for each Gauss-Point */

  /* Iterate over the particles to get the nodal values */
  for(int p = 0 ; p<Np ; p++)
    {
      /*! Define element of the particle */
      Nodes_p = get_Element(p, MPM_Mesh.ListNodes[p], MPM_Mesh.NumberNodes[p]);

      /*! Evaluate the shape function in the coordinates of the particle */
      ShapeFunction_p = compute_ShapeFunction(Nodes_p, MPM_Mesh, FEM_Mesh);
      /*! Get the mass of the particle */
      m_p = MPM_Mesh.Phi.mass.nV[p];

      /*! Get the nodal mass and mommentum */
      for(int I = 0 ; I<Nodes_p.NumberNodes ; I++)
  {
      
    /*! Get the node for the GP */
    Ip = Nodes_p.Connectivity[I];
      
    /*! Evaluate the GP function in the node */
    ShapeFunction_pI = ShapeFunction_p.nV[I];
      
    /*! If this node has a null Value of the SHF continue */
    if(ShapeFunction_pI == 0)
      {
        continue;
      }

    /*! Particle contribution to node I */
    M_Ip = m_p*ShapeFunction_pI;
      
    /*! Nodal stress */
    for(int i = 0 ; i<Ndim*Ndim ; i++)
      {
    S_I.nM[Ip][i] += M_Ip*MPM_Mesh.Phi.Stress.nM[p][i];
      }

    /*! Nodal velocity */
    for(int i = 0 ; i<Ndim ; i++)
      {
        V_I.nM[Ip][i] += M_Ip*MPM_Mesh.Phi.vel.nM[p][i];
      }
      
    /*! Nodal mass matrix */
    M_I.nV[Ip] += M_Ip;
      
  }

      /*! Free the value of the shape functions */
      FreeMat(ShapeFunction_p);
      free(Nodes_p.Connectivity);
    }

  /*! Get nodal values of the velocity and stress */
  for(int I = 0 ; I<Nnodes ; I++)
    {
      if(M_I.nV[I] > TOL_zero)
  {
    for(int i = 0 ; i<Ndim*Ndim ; i++)
      {      
        S_I.nM[I][i] = S_I.nM[I][i]/M_I.nV[I];
      }
    for(int i = 0 ; i<Ndim ; i++)
      {      
        V_I.nM[I][i] = V_I.nM[I][i]/M_I.nV[I];
      }
  }    
    }
 
}

/*******************************************************/

static void update_VelocityField_SV(Matrix S_I,Matrix V_I,Matrix M_I,
           GaussPoint MPM_Mesh,Mesh FEM_Mesh,double Dt)
{
  int Nnodes = FEM_Mesh.NumNodesMesh;
  int Np = MPM_Mesh.NumGP;
  int Ndim = NumberDimensions;
  int Nn_p,Ip,Idx_Mat_p;
  Element Nodes_p;
  Matrix Nodal_Velocity_p; /* Velocity of the element nodes */
  Matrix Gradient_p; /* Shape functions gradients */
  Material Material_p; /* Properties of the particle */
  Tensor Stress_p; 
  Tensor Rate_Stress_p;
  Tensor Rate_Strain_p; 
  Tensor Gradient_pI;
  Tensor InternalForcesDensity_Ip;
  Tensor Rate_InternalForcesDensity_Ip;
  double ShapeFunction_pI;
  double rho_p, m_p, V_p;
  Matrix RHS_1 = MatAllocZ(Nnodes,Ndim);
  Matrix RHS_2 = MatAllocZ(Nnodes,Ndim);

  /* Loop in the particles */
  for(int p = 0 ; p<Np ; p++)
    {
      /* Get the value of the density */
      rho_p = MPM_Mesh.Phi.rho.nV[p];
      m_p = MPM_Mesh.Phi.mass.nV[p];
      V_p = m_p/rho_p;
      Stress_p = memory_to_Tensor(MPM_Mesh.Phi.Stress.nM[p], 2);
      Idx_Mat_p = MPM_Mesh.MatIdx[p];
      Material_p = MPM_Mesh.Mat[Idx_Mat_p];
      
      /* Define a list of nodes for each particle */
      Nn_p = MPM_Mesh.NumberNodes[p];
      Nodes_p = get_Element(p,MPM_Mesh.ListNodes[p],Nn_p);
       
      /* Compute shape functions gradient */
      Gradient_p = compute_ShapeFunction_Gradient(Nodes_p, MPM_Mesh,FEM_Mesh);

      /* Compute rate of stress tensor */
      Nodal_Velocity_p = get_Element_Field(V_I, Nodes_p);
      Rate_Strain_p = compute_RateOfStrain(Nodal_Velocity_p,Gradient_p);
      Rate_Stress_p = alloc_Tensor(2);
      Rate_Stress_p = compute_Stress(Rate_Strain_p,Rate_Stress_p,Material_p);
        
      /* Compute the RHS with the rate of stress tensor */
      for(int I = 0 ; I<Nn_p ; I++)
  {
    /* Pass by reference the nodal gradient to the tensor */
    Gradient_pI = memory_to_Tensor(Gradient_p.nM[I], 1);
    
    /* Compute internal force desnsity and its rate */
    InternalForcesDensity_Ip =
      get_firstOrderContraction_Of(Stress_p, Gradient_pI);
    Rate_InternalForcesDensity_Ip =
      get_firstOrderContraction_Of(Rate_Stress_p, Gradient_pI);

    /* Get the node of the mesh for the contribution */
    Ip = Nodes_p.Connectivity[I];
      
    /* Upgrade the value of the stress tensor */
    for(int i = 0 ; i<Ndim ; i++)
      {
        RHS_1.nM[Ip][i] -= Dt*InternalForcesDensity_Ip.n[i]*V_p;
        RHS_2.nM[Ip][i] -= 0.5*Dt*Dt*Rate_InternalForcesDensity_Ip.n[i]*V_p;
      }
    
    /* Free the internal forces density */
    free_Tensor(InternalForcesDensity_Ip);
    free_Tensor(Rate_InternalForcesDensity_Ip);
  }
      
      /* Free the pointer with the nodal connectivity of the particle */
      free(Nodes_p.Connectivity);
      FreeMat(Gradient_p);
      FreeMat(Nodal_Velocity_p);
      free_Tensor(Rate_Strain_p);   
      free_Tensor(Rate_Stress_p);
    
    }
  
  /* Update the grid nodal momentum */
  for(int I = 0 ; I<Nnodes ; I++)
    {
      if(M_I.nV[I] > TOL_zero)
        {   
    for(int i = 0 ; i<Ndim ; i++)
      {
        V_I.nM[I][i] +=
    (RHS_1.nM[I][i] + RHS_2.nM[I][i])/M_I.nV[I];
      }
  }
    }
  /* Free the matrix with the RHS */
  FreeMat(RHS_1);
  FreeMat(RHS_2);

}

/*******************************************************/

static void update_Particles_SV(Matrix S_I, Matrix V_I,
       GaussPoint MPM_Mesh, Mesh FEM_Mesh,
       double Dt)
{

  int Ndim = NumberDimensions;
  Element Nodes_p; /* Element for each Gauss-Point */
  Matrix N_p; /* Value of the shape-function in the GP */
  double N_pI; /* Nodal value for the GP */
  int Np = MPM_Mesh.NumGP;
  int Nnodes;
  int Ip; /* Index of each tributary node for the GP */

  /* 1º iterate over the Gauss-Points */
  for(int p = 0 ; p<Np ; p++)
    {
      /* 2º Define element of the GP */
      Nnodes = MPM_Mesh.NumberNodes[p];
      Nodes_p = get_Element(p, MPM_Mesh.ListNodes[p], Nnodes);

      /* 3º Evaluate shape function in the GP i */
      N_p = compute_ShapeFunction(Nodes_p, MPM_Mesh, FEM_Mesh);

      /* Set to zero velocity for interpolation */
      for(int i = 0 ; i<Ndim ; i++)
  {
    MPM_Mesh.Phi.vel.nM[p][i] = 0;
  }
      /* Set to zero velocity for interpolation */
      for(int i = 0 ; i<Ndim*Ndim ; i++)
  {
    MPM_Mesh.Phi.Stress.nM[p][i] = 0;
  }
    
      /* 4º Iterate over the nodes of the element */
      for(int I = 0; I<Nnodes; I++)
  {
    /* Node of the GP */
    Ip = Nodes_p.Connectivity[I];
    /* Evaluate the GP function in the node */
    N_pI = N_p.nV[I];
    /* If this node has a null Value of the SHF continue */
    if(fabs(N_pI) <= TOL_zero)
      {
        continue;
      }
    /* Update the particles velocities and position */
    for(int i = 0 ; i<Ndim ; i++)
      {
        MPM_Mesh.Phi.vel.nM[p][i] += N_pI*V_I.nM[Ip][i];
        MPM_Mesh.Phi.x_GC.nM[p][i] += Dt*N_pI*V_I.nM[Ip][i];    
      }
    /* Update partiucle stress */
    for(int i = 0 ; i<Ndim*Ndim; i++)
      {
        MPM_Mesh.Phi.Stress.nM[p][i] += N_pI*S_I.nM[Ip][i];
      }
      
  }
    
      /* 5º Free memory */
      free(Nodes_p.Connectivity);
      FreeMat(N_p);
    }  
}

/*******************************************************/
