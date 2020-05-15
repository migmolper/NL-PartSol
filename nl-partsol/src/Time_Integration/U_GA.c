#include "nl-partsol.h"

void U_GA(Mesh FEM_Mesh, GaussPoint MPM_Mesh)
{

  /*!
    Some auxiliar variables for the outputs
  */
  Matrix List_Fields;

  /*!
    Time step 
  */
  int TimeStep;

  /*!
    Control parameters of the generalized-alpha algorithm 
    all the parameters are controled by a simple parameter :
    SpectralRadius 
  */
  Time_Int_Params Params;
  Params.GA_alpha =
    (2*SpectralRadius-1)/(1+SpectralRadius);
  Params.GA_beta = (5-3*SpectralRadius)/
    (pow((1+SpectralRadius),2)*(2-SpectralRadius));
  Params.GA_gamma =
    3/2 - Params.GA_alpha;

  /*!
    Auxiliar variable for the nodal kinetics
    Nodal_Kinetics = {mass, a0, a1, v}
  */
  int Nnodes = FEM_Mesh.NumNodesMesh;
  int Ndim = NumberDimensions;
  Matrix V_I;
  Matrix Nodal_Kinetics;

  /*!
    Nodal forces for the balance 
  */
  Matrix F_I = get_RowFrom(Ndim,Nnodes,NULL);
  Matrix R_I = get_RowFrom(Ndim,Nnodes,NULL);

  
  puts("*************************************************");
  puts(" First step : Get the nodal kinetics");
  puts(" \t WORKING ...");
  Nodal_Kinetics = GetNodalKinetics(MPM_Mesh,FEM_Mesh);
  /* V_I =  MatAssign(Nnodes,Ndim,NAN,NULL, */
  /* 		   (double**)malloc(Nnodes*sizeof(double *))); */
  for(int i = 0 ; i<Ndim ; i++)
    {
      V_I.nM[i] = Nodal_Kinetics.nM[1+2*Ndim+i];
    }
  puts(" \t DONE !!! \n");


  for(TimeStep = 0 ; TimeStep<NumTimeStep ; TimeStep++ )
    {

      print_Status("*************************************************",TimeStep);
      DeltaTimeStep = DeltaT_CFL(MPM_Mesh, FEM_Mesh.DeltaX);
      print_step(TimeStep,DeltaTimeStep);
      
      print_Status("*************************************************",TimeStep);
      print_Status("First step : Compute nodal kinetics ... WORKING",TimeStep);
      /* BCC_Nod_VALUE(FEM_Mesh, V_I, TimeStep); */
      print_Status("DONE !!!",TimeStep);

      if(TimeStep % ResultsTimeStep == 0)
	{
	  /*! 
	    Print Nodal values after appling the BCCs 
	  */
	  WriteVtk_FEM("Mesh",FEM_Mesh,Nodal_Kinetics,
		       (int)TimeStep/ResultsTimeStep);
	  /*!
	    Print GPs results
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
      print_Status(" Third step : Update kinetics ... WORKING",TimeStep);
      GA_UpdateNodalKinetics(FEM_Mesh, Nodal_Kinetics, F_I, Params);
      print_Status("DONE !!!",TimeStep);

      print_Status("*************************************************",TimeStep);
      print_Status("Four step : Update particles lagrangian ... WORKING",TimeStep);
      update_Particles_GA(MPM_Mesh, FEM_Mesh, Nodal_Kinetics, Params);
      LocalSearchGaussPoints(MPM_Mesh, FEM_Mesh);
      print_Status("DONE !!!",TimeStep);
      
      print_Status("*************************************************",TimeStep);
      print_Status("Five step : Reset nodal values ... WORKING",TimeStep);
      FreeMat(F_I);
      print_Status("DONE !!!",TimeStep);

    }
  
}
