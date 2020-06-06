#include "nl-partsol.h"

void U_PCE(Mesh FEM_Mesh, GaussPoint MPM_Mesh, int InitialStep)
{

  int Ndim = NumberDimensions;
  int Nnodes = FEM_Mesh.NumNodesMesh;

  /*!
    Control parameters of the generalized-alpha algorithm 
    all the parameters are controled by a simple parameter :
    SpectralRadius 
  */
  Time_Int_Params Params;
  Params.GA_gamma = 0.5;

  /*!
    Auxiliar variable for the mass and momentum 
  */
  Matrix M_I;
  Matrix V_I;
  Matrix F_I;
  Matrix R_I;

  
  for(int TimeStep = InitialStep ; TimeStep<NumTimeStep ; TimeStep++ )
    {
      print_Status("*************************************************",TimeStep);
      DeltaTimeStep = DeltaT_CFL(MPM_Mesh, FEM_Mesh.DeltaX);
      print_step(TimeStep,DeltaTimeStep);

      print_Status("*************************************************",TimeStep);
      print_Status("First step : Predictor stage ... WORKING",TimeStep);
      M_I = compute_NodalMass(MPM_Mesh, FEM_Mesh);      
      V_I = compute_VelocityPredictor(MPM_Mesh,FEM_Mesh,V_I, M_I,
				      Params, DeltaTimeStep);
      imposse_NodalVelocity(FEM_Mesh,V_I,TimeStep);
      print_Status("DONE !!!",TimeStep);
    
      print_Status("*************************************************",TimeStep);
      print_Status("Second step : Compute equilibrium ... WORKING",TimeStep);
      F_I = compute_equilibrium_U(V_I,MPM_Mesh,FEM_Mesh,TimeStep); 
      R_I = compute_Reactions(FEM_Mesh, F_I);
      print_Status("DONE !!!",TimeStep);
    
      print_Status("*************************************************",TimeStep);
      print_Status(" Third step : Corrector stage ... WORKING",TimeStep);
      V_I = compute_VelocityCorrector(FEM_Mesh,V_I,F_I, M_I,Params,DeltaTimeStep);
      print_Status("DONE !!!",TimeStep);

      print_Status("*************************************************",TimeStep);
      print_Status("Four step : Update particles lagrangian ... WORKING",TimeStep);
      update_Particles_PCE(MPM_Mesh, FEM_Mesh, M_I, V_I, F_I,DeltaTimeStep);
      LocalSearchGaussPoints(MPM_Mesh,FEM_Mesh);
      print_Status("DONE !!!",TimeStep);

      if(TimeStep % ResultsTimeStep == 0)
	{
	  /*!
	    Print Nodal values after appling the BCCs
	  */
	  WriteVtk_FEM("Mesh",FEM_Mesh,R_I,TimeStep);
	  /*!
	    Print particle results 
	  */
	  WriteVtk_MPM("MPM_VALUES",MPM_Mesh,"ALL",TimeStep,ResultsTimeStep);
	}

      print_Status("*************************************************",TimeStep);
      print_Status("Five step : Reset nodal values ... WORKING",TimeStep);
      FreeMat(M_I);
      FreeMat(V_I);
      FreeMat(F_I);
      FreeMat(R_I);
      print_Status("DONE !!!",TimeStep);

    } 
}
