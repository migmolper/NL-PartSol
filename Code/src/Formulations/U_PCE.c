#include "grams.h"

void U_PCE(Mesh FEM_Mesh, GaussPoint MPM_Mesh)
/*!
 * Explicit predictor-corrector gamma = 0.5 and beta = 0.25
 */
{

  /* Some auxiliar variables for the outputs */
  Matrix List_Fields;

  int Ndim = NumberDimensions;
  int Nnodes = FEM_Mesh.NumNodesMesh;

  /* Control parameters of the generalized-alpha algorithm 
     all the parameters are controled by a simple parameter :
     SpectralRadius */
  Time_Int_Params Params;
  Params.GA_gamma = 0.5;

  /*********************************************************************/
  /***** INITIALIZE AUXILIAR STRUCTURES TO STORE NODAL INFORMATION *****/
  /*********************************************************************/

  /* Auxiliar variable for the mass and momentum */
  Matrix M_I;
  Matrix V_I;
  Matrix F_I = MatAllocZ(Nnodes,Ndim);
  Matrix R_I;
  
  /*********************************************************************/
  /*********************************************************************/

  for(int TimeStep = 0 ; TimeStep<NumTimeStep ; TimeStep++ ){

    puts("*************************************************");
    DeltaTimeStep = DeltaT_CFL(MPM_Mesh, FEM_Mesh.DeltaX);
    printf("***************** STEP : %i , DeltaT : %f \n",
	   TimeStep,DeltaTimeStep);
    
    puts("*****************************************");
    puts(" First step : Predictor stage ... WORKING");
    M_I = compute_NodalMass(MPM_Mesh, FEM_Mesh);
    V_I = compute_VelocityPredictor(MPM_Mesh,FEM_Mesh,V_I, M_I,
				    Params, DeltaTimeStep);
    imposse_NodalVelocity(FEM_Mesh,V_I,TimeStep);
    puts(" DONE !!!");
    

    puts("*************************************************");
    puts(" Second step : Update local state ... WORKING");
    update_LocalState(V_I, MPM_Mesh, FEM_Mesh, DeltaTimeStep);
    puts(" DONE !!!");
       
    puts("*************************************************");
    puts(" Third step : Compute forces ... WORKING");
    F_I = MatAllocZ(Nnodes,Ndim);    
    puts(" \t Compute internal forces");
    F_I = compute_InternalForces(F_I, MPM_Mesh, FEM_Mesh);    
    puts(" \t DONE !!!");
    puts(" \t Compute body forces");
    F_I = compute_BodyForces(F_I, MPM_Mesh, FEM_Mesh, TimeStep);
    puts(" \t DONE !!!");
    puts(" \t Compute contact forces");
    F_I = compute_ContacForces(F_I, MPM_Mesh, FEM_Mesh, TimeStep);
    puts(" \t DONE !!!");
    puts(" \t Compute reactions");
    R_I = compute_Reactions(FEM_Mesh, F_I);
    puts(" \t DONE !!!");
    
    puts("****************************************");
    puts(" Four step : Corrector stage ... WORKING");
    V_I = compute_VelocityCorrector(FEM_Mesh,V_I,F_I, M_I,Params,DeltaTimeStep);
    puts(" DONE !!!"); 

    puts("***************************************************");
    puts(" Five step : Update particles lagrangian ... WORKING");
    update_Particles_PCE(MPM_Mesh, FEM_Mesh, M_I, V_I, F_I,DeltaTimeStep);
    LocalSearchGaussPoints(MPM_Mesh,FEM_Mesh);
    puts(" DONE !!!");

    if(TimeStep % ResultsTimeStep == 0){
      /* Print Nodal values after appling the BCCs */
      WriteVtk_FEM("Mesh",FEM_Mesh,R_I,
      		   (int)TimeStep/ResultsTimeStep);
      /* Print GPs results */
      WriteVtk_MPM("MPM_VALUES",MPM_Mesh,List_Fields,
      		   (int)TimeStep/ResultsTimeStep);
    }

    puts("********************************************************");
    puts(" Six step : Reset nodal values of the mesh ... WORKING");
    FreeMat(M_I);
    FreeMat(V_I);
    FreeMat(F_I);
    FreeMat(R_I);
    puts(" DONE !!!");

  } /* End of temporal integration */
}
