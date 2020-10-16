#include "nl-partsol.h"

/*
  Auxiliar functions 
*/
static Tensor compute_small_strain_tensor(Tensor, Tensor);
static void   compute_volumetric_deviatoric_stress_tensor(double *, Tensor, Tensor, Material);
static double compute_yield_surface(Tensor, Tensor, Material);
static double compute_hardening_modulus(double, double, double, double);

/**************************************************************/

Tensor viscoplastic_Drucker_Prager_Sanavia(Tensor grad_e, Tensor C, Tensor F_plastic, Tensor F, double J, Material MatProp)
	/*	
		Radial returning algorithm
	*/
{

	
	Tensor E_elastic;
	double p_trial;
	Tensor s_trial;
	double Phi;
	double H;
	double p_lim;
	double d_Phi;
	Tensor sigma_k1;
	Tensor D_E_plastic;

	double cohesion_reference = MatProp.cohesion_reference;
	double cohesion_reference = MatProp.cohesion_reference;
	double hardening_exp = MatProp.hardening_exp;

	/*	
		Calculation of the small strain tensor
	*/
	E_elastic = compute_small_strain_tensor(Tensor C, Tensor F_plastic);

	/*
		Elastic predictor : Volumetric and deviatoric stress measurements
	*/
	compute_volumetric_deviatoric_stress_tensor(&p_trial, s_trial, E_elastic, MatProp);

	/*
		Yield condition : Starting from incremental plastic strain equal to zero
	*/
	Phi = compute_yield_surface(p_trial, s_trial, MatProp);
	if(Phi >= 0)
	{
		H = compute_hardening_modulus();

		p_lim = compute_limit_between_classic_apex_algorithm();

		/*
			Classical plastic iterator
		*/
		if(p_trial <= p_lim)
		{
			while(iter < iter_max_DP)
			{

				d_Phi = compute_derivative_yield_surface_classical();

				delta_Gamma = compute_increment_plastic_strain();

				EPS_k1 = compute_equivalent_plastic_strain_classical();

				c_k1 = compute_cohesion_modulus();

				H = compute_hardening_modulus(EPS_k1, E_plastic_reference, cohesion_reference, hardening_exp);

				Phi = compute_yield_surface_classical();

				/*
					Check convergence
				*/
				if(Phi < TOL_Drucker_Prager)
				{
					iter++;
				}
				else
				{
					break;
				}

			}

			/*
				update
			*/
			D_E_plastic = compute_increment_plastic_strain_classical();
			sigma_k1 = compute_finite_stress_tensor_classical();

		}
		/*
			Apex plastic iterator
		*/
		else
		{

			while(iter < iter_max_DP)
			{

				d_Phi = compute_derivative_yield_surface_apex();

				delta_Gamma = compute_increment_plastic_strain();

				Phi = compute_yield_surface_apex();

				/*
					Check convergence
				*/
				if(Phi < TOL_Drucker_Prager)
				{
					iter++;
				}
				else
				{
					break;
				}

			}

			/*
				update
			*/

			EPS_k1 = compute_equivalent_plastic_strain_apex();

			D_E_plastic = compute_increment_plastic_strain_apex();

			sigma_k1 = compute_finite_stress_tensor_apex();

		}

	}

	/*
		Step 
	*/
	update_plastic_deformation_gradient(D_E_plastic,F_plastic);

	/*
		Get the stress tensor in the deformed configuration
	*/
	push_backward_tensor__Particles__(grad_e, sigma_k1, F);
	
	for(int i = 0 ; i<Ndim ; i++)
	{
		for(int j = 0 ; j<Ndim ; j++)
		{
			grad_e.N[i][j] *= J;
		}
	}

	/*
		Free memory
	*/
	free__TensorLib__(E_elastic);
	free__TensorLib__(s_trial);
	free__TensorLib__(sigma_k1);
	free__TensorLib__(Increment_E_plastic);

	/*
		Return stress tensor
	*/
	return grad_e;
}

/**************************************************************/

static Tensor compute_small_strain_tensor(Tensor C, Tensor F_plastic)
{
	Tensor C_elastic;
	Tensor E_elastic;

	/*
		Compute the trial elastic right Cauchy-Green tensor using the intermediate configuration.
	*/
	push_backward_tensor__Particles__(C_elastic, C, F_plastic);

	/*
		Use the approach of Ortiz and Camacho to compute the elastic infinitesimal strain tensor.
	*/
	Eps_elastic = logarithmic_elastic_strains(C_elastic);

	/*	
		Free memory
	*/
	free__TensorLib__(C_elastic);


	return E_elastic;
}

/**************************************************************/

static void compute_volumetric_deviatoric_stress_tensor(double * p_trial, Tensor s_tria, Tensor E_elastic, Material MatProp)
{
	volumetric_desviatoric_decomposition__TensorLib__(Eps_elastic,Eps_elastic_vol,Eps_elastic_dev);
}

/**************************************************************/

static double compute_hardening_modulus(double Equiv_plastic_strain, double E_plastic_reference, double cohesion_reference, double hardening_exp)
{
	double c0 = cohesion_reference;
	double N_eps = hardening_exp;
	double EPS = Equiv_plastic_strain;
	double E0 = E_plastic_reference;

	double H = (c0/(N_eps*E0))*pow((1 + EPS/E0),(1/N_eps - 1));

	return H;
}

/**************************************************************/


