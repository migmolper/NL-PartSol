#include "nl-partsol.h"

void U_FE(Mesh FEM_Mesh, GaussPoint MPM_Mesh)
{

  /*!
    Some auxiliar variables for the outputs 
  */
  Matrix List_Fields;

  /*!
    Integer variables 
  */
  int Ndim = NumberDimensions;
  int Nnodes = FEM_Mesh.NumNodesMesh;

  /*!
    Auxiliar variable for the mass and momentum 
  */
  Matrix Phi_I;
  strcpy(Phi_I.Info,"MOMENTUM;MASS");
  Matrix V_I;
  Matrix F_I;
  Matrix R_I;
  

  for(int TimeStep = 0 ; TimeStep<NumTimeStep ; TimeStep++ )
    {

      print_Status("*************************************************",TimeStep);
      DeltaTimeStep = DeltaT_CFL(MPM_Mesh, FEM_Mesh.DeltaX);
      print_step(TimeStep,DeltaTimeStep);

      print_Status("*************************************************",TimeStep);
      print_Status(" First step : Get the nodal fields ... WORKING",TimeStep);
      Phi_I = compute_NodalMomentumMass(MPM_Mesh,FEM_Mesh);
      imposse_NodalMomentum(FEM_Mesh,Phi_I,TimeStep);    
      V_I = compute_NodalVelocity(FEM_Mesh, Phi_I);
      print_Status("DONE !!!",TimeStep);

      if(TimeStep % ResultsTimeStep == 0)
	{
	  /*!
	    Print Nodal values after appling the BCCs
	  */
	  WriteVtk_FEM("Mesh",FEM_Mesh,Phi_I,
		       (int)TimeStep/ResultsTimeStep);
	  /*!
	    Print particle results 
	  */
	  WriteVtk_MPM("MPM_VALUES",MPM_Mesh,List_Fields,
		       (int)TimeStep/ResultsTimeStep);
	}

      print_Status("*************************************************",TimeStep);
      print_Status("Second step : Compute equilibrium ... WORKING",TimeStep);
      F_I = compute_equilibrium_U(V_I,MPM_Mesh,FEM_Mesh,TimeStep); 
      R_I = compute_Reactions(FEM_Mesh, F_I);
      print_Status("DONE !!!",TimeStep);

      print_Status("*************************************************",TimeStep);
      print_Status(" Third step : Update nodal momentum ... WORKING",TimeStep);
      update_NodalMomentum(FEM_Mesh,Phi_I,F_I);
      print_Status("DONE !!!",TimeStep);
    
      print_Status("*************************************************",TimeStep);
      print_Status(" Four step : Update lagrangian ... WORKING",TimeStep);
      update_Particles_FE(MPM_Mesh, FEM_Mesh, Phi_I, F_I, DeltaTimeStep);
      LocalSearchGaussPoints(MPM_Mesh,FEM_Mesh);
      print_Status("DONE !!!",TimeStep);
    
      print_Status("*************************************************",TimeStep);
      print_Status(" Five step : Reset nodal values ... WORKING",TimeStep);
      FreeMat(Phi_I);
      FreeMat(V_I);
      FreeMat(F_I);
      FreeMat(R_I);
      print_Status("DONE !!!",TimeStep);

    } 

}
