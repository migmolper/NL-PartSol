#include "nl-partsol.h"

void SV_TG(Mesh FEM_Mesh, GaussPoint MPM_Mesh)
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
