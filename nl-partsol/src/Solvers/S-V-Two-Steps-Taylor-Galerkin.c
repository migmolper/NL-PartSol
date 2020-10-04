#include "nl-partsol.h"


/*
  Auxiliar functions 
*/
static void update_StressField_SV(Matrix,Matrix,Matrix,GaussPoint,Mesh,double);



void SV_Two_Steps_Taylor_Galerkin(Mesh FEM_Mesh, GaussPoint MPM_Mesh)
{

  /* Some auxiliar variables for the outputs */
  Matrix List_Fields;
  /* Time step */
  int TimeStep;

  int Ndim = NumberDimensions;
  int Nnodes = FEM_Mesh.NumNodesMesh;

  Matrix V_I, S_I, Div_S_I, M_I;

  for(TimeStep = 0 ; TimeStep<NumTimeStep ; TimeStep++ )
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

void update_StressField_SV(Matrix S_I,Matrix V_I,Matrix M_I,
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
