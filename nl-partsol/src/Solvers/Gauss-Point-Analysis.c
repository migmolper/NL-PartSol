#include "nl-partsol.h"


/*
  Call global variables
*/
int NumTimeStep;
Event * Out_Gauss_Point_evolution_csv;
int Number_Out_Gauss_Point_evolution_csv;

/*
  Auxiliar functions and variables
*/

static char Error_message[MAXW];

static void standard_error();

/**************************************************************/

void NonLinear_Gauss_Point_Analysis(Particle PointAnalysis)
{

  /*
    Integer variables 
  */
  int Ndim = NumberDimensions;

  /* Particles variables */
  Tensor S_k; /* 2ª Piola-Kirchhoff stress tensor at k step */
  Tensor F_k; /* Deformation gradient tensor at k step */
  Tensor C_k; /* Right Cauchy-Green tensor at k step */
  Tensor F_m1_plastic_k; /* Plastic deformation gradient at k step */
  double J_k; /* Jacobian */

  State_Parameters Input_SP;
  State_Parameters Output_SP;

  for(int k = 1 ; k<NumTimeStep ; k++)
  {   
  	if(strcmp(PointAnalysis.Mat[0].Type,"Saint-Venant-Kirchhoff") == 0)
  	{
      Input_SP.Stress = PointAnalysis.Phi.Stress.nM[k];
      Output_SP = compute_1PK_Stress_Tensor_Saint_Venant_Kirchhoff(Input_SP, PointAnalysis.Mat[0]);
    }
    else if(strcmp(PointAnalysis.Mat[0].Type,"Neo-Hookean-Wriggers") == 0)
    {
      Input_SP.Stress = PointAnalysis.Phi.Stress.nM[k];
      Input_SP.J = PointAnalysis.Phi.J.nV[k];
      Output_SP = compute_1PK_Stress_Tensor_Neo_Hookean_Wriggers(Input_SP, PointAnalysis.Mat[0]);
    }
    else if(strcmp(PointAnalysis.Mat[0].Type,"Von-Mises") == 0)
    {
      Input_SP.Stress = PointAnalysis.Phi.Stress.nM[k];
      Input_SP.F_m1_plastic_p = PointAnalysis.Phi.F_m1_plastic.nM[k-1];
      Input_SP.EPS = PointAnalysis.Phi.EPS.nV[k-1];
      Input_SP.Back_stress = PointAnalysis.Phi.Back_stress.nM[k-1];
      Input_SP.F_n1_p = PointAnalysis.Phi.F_n.nM[k];
  
      if(strcmp(PointAnalysis.Mat[0].Plastic_Solver,"Backward-Euler") == 0)
      {
        Output_SP = finite_strain_plasticity(Input_SP,PointAnalysis.Mat[0],Von_Mises_backward_euler);
      }
      else if(strcmp(PointAnalysis.Mat[0].Plastic_Solver,"Forward-Euler") == 0)
      {
        Output_SP = finite_strain_plasticity(Input_SP,PointAnalysis.Mat[0],Von_Mises_forward_euler);
      }
      else
      {
        fprintf(stderr,"%s : %s %s %s \n","Error in stress_integration__Particles__()",
      "The solver",PointAnalysis.Mat[0].Type,"has not been yet implemented");
        exit(EXIT_FAILURE);
      }

      for(int i = 0 ; i<Ndim*Ndim ; i++)
      {
        PointAnalysis.Phi.F_m1_plastic.nM[k][i] = PointAnalysis.Phi.F_m1_plastic.nM[k-1][i];
        PointAnalysis.Phi.Back_stress.nM[k][i] = PointAnalysis.Phi.Back_stress.nM[k-1][i];
      }

      PointAnalysis.Phi.EPS.nV[k] = Output_SP.EPS;

    }
    else if((strcmp(PointAnalysis.Mat[0].Type,"Drucker-Prager-Plane-Strain") == 0) || 
              (strcmp(PointAnalysis.Mat[0].Type,"Drucker-Prager-Outer-Cone") == 0))
    {
      Input_SP.F_m1_plastic_p = PointAnalysis.Phi.F_m1_plastic.nM[k];
      Input_SP.Cohesion = PointAnalysis.Phi.cohesion.nV[k];
      Input_SP.EPS = PointAnalysis.Phi.EPS.nV[k];
      Input_SP.F_n1_p = PointAnalysis.Phi.F_n1.nM[k];
  
      if(strcmp(PointAnalysis.Mat[0].Plastic_Solver,"Backward-Euler") == 0)
      {
        Output_SP = finite_strain_plasticity(Input_SP,PointAnalysis.Mat[0],Drucker_Prager_backward_euler);
      }
      else
      {
        fprintf(stderr,"%s : %s %s %s \n","Error in stress_integration__Particles__()",
      "The solver",PointAnalysis.Mat[0].Type,"has not been yet implemented");
        exit(EXIT_FAILURE);
      }

      PointAnalysis.Phi.cohesion.nV[k] = Output_SP.Cohesion;
      PointAnalysis.Phi.EPS.nV[k] = Output_SP.EPS;
    }
	else
	{

	  sprintf(Error_message,"%s %s %s","The material",PointAnalysis.Mat[0].Type,"has not been yet implemnented");
    standard_error(); 

  }

	/* Output stress trajectory */
	for(int i = 0 ; i<Number_Out_Gauss_Point_evolution_csv; i++)
  	{

    	Gauss_Point_evolution__InOutFun__(PointAnalysis,Out_Gauss_Point_evolution_csv[i],"Particle_evolution", k,i);
       
	  }

  }
      
}

/***************************************************************************/

static void standard_error()
{
  fprintf(stderr,"%s : %s !!! \n",
     "Error in NonLinear_Gauss_Point_Analysis",Error_message);
    exit(EXIT_FAILURE);
}

/**************************************************************/

