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

void NonLinear_Gauss_Point_Analysis(GaussPoint PointAnalysis)
{

  /*
    Integer variables 
  */
  int Ndim = NumberDimensions;

  /* Particles variables */
  Tensor S_k; /* 2ª Piola-Kirchhoff stress tensor at k step */
  Tensor F_k; /* Deformation gradient tensor at k step */
  Tensor C_k; /* Right Cauchy-Green tensor at k step */
  Tensor F_plastic_k; /* Plastic deformation gradient at k step */
  double J_k; /* Jacobian */
  Plastic_status Input_Plastic_Parameters;
  Plastic_status Output_Plastic_Parameters;

  for(int k = 0 ; k<NumTimeStep ; k++)
  {

	  S_k = memory_to_tensor__TensorLib__(PointAnalysis.Phi.Stress.nM[k],2);     
  	F_k = memory_to_tensor__TensorLib__(PointAnalysis.Phi.F_n.nM[k],2);
  	C_k = right_Cauchy_Green__Particles__(F_k);

  	if(strcmp(PointAnalysis.Mat[0].Type,"Saint-Venant-Kirchhoff") == 0)
  	{

  		S_k = compute_SPK__SaintVenantKirchhoff__(S_k, C_k, PointAnalysis.Mat[0]);

    }
    else if(strcmp(PointAnalysis.Mat[0].Type,"Neo-Hookean-Wriggers") == 0)
    {

    	J_k = I3__TensorLib__(F_k);
      S_k = grad_energy_Neo_Hookean_Wriggers(S_k, C_k, J_k, PointAnalysis.Mat[0]);

    }
    else if(strcmp(PointAnalysis.Mat[0].Type,"Von-Mises") == 0)
    {

    	J_k = I3__TensorLib__(F_k);
      F_plastic_k = memory_to_tensor__TensorLib__(PointAnalysis.Phi.F_plastic.nM[k],2);
      Input_Plastic_Parameters.Cohesion = PointAnalysis.Phi.cohesion.nV[k];
      Input_Plastic_Parameters.EPS = PointAnalysis.Phi.EPS.nV[k];

      Output_Plastic_Parameters = finite_strains_plasticity_Von_Mises(S_k, C_k, F_plastic_k, F_k, 
                                                                          Input_Plastic_Parameters, PointAnalysis.Mat[0], J_k);

      /* Update variables (cohesion and EPS) */
      PointAnalysis.Phi.cohesion.nV[k] = Output_Plastic_Parameters.Yield_stress;
      PointAnalysis.Phi.EPS.nV[k] = Output_Plastic_Parameters.EPS;

    }
    else if((strcmp(PointAnalysis.Mat[0].Type,"Drucker-Prager-Plane-Strain") == 0) || 
              (strcmp(PointAnalysis.Mat[0].Type,"Drucker-Prager-Outer-Cone") == 0))
    {

    	J_k = I3__TensorLib__(F_k);
    	F_plastic_k = memory_to_tensor__TensorLib__(PointAnalysis.Phi.F_plastic.nM[k],2);
    	Input_Plastic_Parameters.Cohesion = PointAnalysis.Phi.cohesion.nV[k];
    	Input_Plastic_Parameters.EPS = PointAnalysis.Phi.EPS.nV[k];

      Output_Plastic_Parameters = finite_strains_plasticity_Drucker_Prager_Sanavia(S_k, C_k, F_plastic_k,
                                                                                    Input_Plastic_Parameters, PointAnalysis.Mat[0]);

      /* Update variables (cohesion and EPS) */
      PointAnalysis.Phi.cohesion.nV[k] = Output_Plastic_Parameters.Cohesion;
      PointAnalysis.Phi.EPS.nV[k] = Output_Plastic_Parameters.EPS;

    }
	else
	{

	  sprintf(Error_message,"%s %s %s","The material",PointAnalysis.Mat[0].Type,"has not been yet implemnented");
    standard_error(); 

  }

  free__TensorLib__(C_k);

	/* Output stress trajectory */
	for(int i = 0 ; i<Number_Out_Gauss_Point_evolution_csv; i++)
  	{

    	Gauss_Point_evolution__InOutFun__(PointAnalysis,Out_Gauss_Point_evolution_csv[i],"GaussPoint_evolution", k,i);
       
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

