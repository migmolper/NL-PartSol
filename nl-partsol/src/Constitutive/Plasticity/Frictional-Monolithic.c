
#ifdef __linux__
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <stdbool.h>
#include <lapacke.h>
#elif __APPLE__
#include <stdio.h>
#include <string.h>
#include <stdbool.h>
#include <Accelerate/Accelerate.h>
#endif

#include "Macros.h"
#include "Types.h"
#include "Globals.h"

/*
  Auxiliar variables
*/
typedef struct
{
  double alpha;
  double m; 
  double pa; 
  double c0; 
  double a1; 
  double a2; 
  double a3;
  double CC[9];

} Model_Parameters;

/*
  Define local global variables
*/
int Particle_Idx;
double Error0;
bool Is_Matsuoka_Nakai;
bool Is_Lade_Duncan;
bool Is_Modified_Lade_Duncan;

/*
  Auxiliar functions 
*/
static Model_Parameters fill_model_paramters(Material Material);
static double eval_K1(double,double,double,double,double);
static double eval_K2(double,double,double,double,double);
static void eval_d_K2_d_stress(double *,double,double,double,double);
static double eval_b1(double,double,double,double,double);
static double eval_b2(double,double,double,double,double);
static void eval_d_b2_d_stress(double *,double *,double,double,double,double,double);
static void eval_kappa(double *,double,double,double,double,double,double);
static void eval_d_kappa1_d_stress(double *,double,double,double,double,double);
static double eval_d_kappa1_d_lambda(double,double,double,double,double);
static double eval_f(double,double);
static void eval_d_f_d_stress(double *,double *,double,double);
static double eval_g(double,double);
static void eval_d_g_d_stress(double *,double *,double,double);
static void eval_dd_g_dd_stress(double *,double *,double,double);
static double eval_Yield_Function(double,double,double,double,double,double,double);
static void eval_d_Yield_Function_d_stress(double *,double *,double,double,double,double,double,double,double);
static double eval_d_Yield_Function_d_kappa1(double,double,double,double,double,double);
static double eval_Plastic_Potential(double,double,double,double,double,double,double);
static void eval_d_Plastic_Potential_d_stress(double *,double *,double,double,double,double,double,double,double);
static void eval_dd_Plastic_Potential_dd_stress(double *,double *,double,double,double,double,double,double,double);
static void eval_dd_Plastic_Potential_d_stress_d_kappa2(double *, double *,double,double,double,double,double,double);
static void eval_strain(double *,double *,double *);
static double assemble_residual(double *,double *,double *,double *,double *,double,double,Model_Parameters);
static double compute_condition_number(double *);
static bool check_convergence(double,int,int);
static void assemble_tangent_matrix(double *,double *,double *,double *,double,double,Model_Parameters);
static void solver(double *,double *,double *);
static void update_variables(double *,double *,double *,double *,double *,double,double,double);
static State_Parameters fill_Outputs(double *,double *,double *,double,double,double);

/**************************************************************/

State_Parameters Frictional_Monolithic(
  State_Parameters Inputs_SP,
  Material MatProp)
/*	
	Monolithic algorithm for a smooth Mohr-Coulomb plastic criterium
*/
{

  /*
    Define output state parameters
  */
  State_Parameters Outputs_VarCons;

  /*
    Get material variables
  */
  Model_Parameters Params = fill_model_paramters(MatProp);

  /*
    Get the index of the particle
  */
  Particle_Idx = Inputs_SP.Particle_Idx;

  /*
    Assign values from the inputs state parameters
  */
  double * Increment_E_plastic = Inputs_SP.Increment_E_plastic;
  double * Stress_k  = Inputs_SP.Stress;
  double Plastic_Flow_k[3] = {0.0, 0.0, 0.0};
  double kappa_k[2] = {Inputs_SP.Kappa, Params.alpha*Inputs_SP.Kappa};
  double Lambda_n = Inputs_SP.Equiv_Plast_Str;
  double Lambda_k = Lambda_n;
  double delta_lambda_k = 0.0;

  /*
    Duplicated variables for the line search algorithm
  */
  double Stress_k2[3];
  double kappa_k2[2];
  double Lambda_k2 = Lambda_k;
  double delta_lambda_k2 = delta_lambda_k;

  /*
    Initialize Newton-Raphson solver
  */
  double Strain_e_tri[3];
  double Residual[5];
  double D_Residual[5];
  double Tangent_Matrix[25]; 
  double Norm_Residual_1;
  double Norm_Residual_2;
  double delta = 1;
  int MaxIter = Max_Iterations_Radial_Returning;
  int Iter = 0;
  bool Convergence = false;

  /*
    Compute invariants
  */
  double I1 = Stress_k[0] + Stress_k[1] + Stress_k[2];
  double I2 = Stress_k[0]*Stress_k[1] + Stress_k[1]*Stress_k[2] + Stress_k[0]*Stress_k[2];
  double I3 = Stress_k[0]*Stress_k[1]*Stress_k[2];

  /*
    Check yield condition
  */
  if(eval_Yield_Function(Params.c0,kappa_k[0],Params.pa,I1,I2,I3,Params.m) > 0.0)
  {
    /*
      Compute elastic trial with the inverse elastic relation
    */
    eval_strain(Strain_e_tri,Stress_k,Params.CC);
    
    /*
      Newton-Rapson with line search
    */
    Norm_Residual_1 = assemble_residual(Residual,Stress_k,Strain_e_tri,Plastic_Flow_k,kappa_k,delta_lambda_k,Lambda_k,Params); 

    Convergence = check_convergence(Norm_Residual_1,Iter,MaxIter);
        
    while(Convergence == false)
    {

      if(Convergence == false)
      {

        delta = 1;

        Iter++;
        
        assemble_tangent_matrix(Tangent_Matrix,Plastic_Flow_k,Stress_k,kappa_k,delta_lambda_k,Lambda_k,Params);
        
        solver(Tangent_Matrix,Residual,D_Residual);

        Stress_k2[0] = Stress_k[0]; 
        Stress_k2[1] = Stress_k[1];
        Stress_k2[2] = Stress_k[2];
        kappa_k2[0] = kappa_k[0];
        kappa_k2[1] = kappa_k[1];
        Lambda_k2 = Lambda_k;
        delta_lambda_k2 = delta_lambda_k;

        update_variables(D_Residual,Stress_k2,kappa_k2,&delta_lambda_k2,&Lambda_k2,Lambda_n,Params.alpha,delta);

        Norm_Residual_2 = assemble_residual(Residual,Stress_k2,Strain_e_tri,Plastic_Flow_k,kappa_k2,delta_lambda_k2,Lambda_k2,Params);
          
        if((Norm_Residual_2 - Norm_Residual_1) > TOL_Radial_Returning)
        {
          while((Norm_Residual_2 - Norm_Residual_1) > TOL_Radial_Returning)
          {

            delta = pow(delta,2)*0.5*Norm_Residual_1/(Norm_Residual_2 - delta*Norm_Residual_1 + Norm_Residual_1);

            if(delta < TOL_Radial_Returning)
            {
              Stress_k2[0] = Stress_k[0]; 
              Stress_k2[1] = Stress_k[1];
              Stress_k2[2] = Stress_k[2];
              kappa_k2[0] = kappa_k[0];
              kappa_k2[1] = kappa_k[1];
              Lambda_k2 = Lambda_k;
              delta_lambda_k2 = delta_lambda_k;
              break;
            }
            else
            {
              Stress_k2[0] = Stress_k[0]; 
              Stress_k2[1] = Stress_k[1];
              Stress_k2[2] = Stress_k[2];
              kappa_k2[0] = kappa_k[0];
              kappa_k2[1] = kappa_k[1];
              Lambda_k2 = Lambda_k;
              delta_lambda_k2 = delta_lambda_k;

              update_variables(D_Residual,Stress_k2,kappa_k2,&delta_lambda_k2,&Lambda_k2,Lambda_n,Params.alpha,delta);

              Norm_Residual_2 = assemble_residual(Residual,Stress_k2,Strain_e_tri,Plastic_Flow_k,kappa_k2,delta_lambda_k2,Lambda_k2,Params);
            }

          }

        }

        Stress_k[0] = Stress_k2[0]; 
        Stress_k[1] = Stress_k2[1];
        Stress_k[2] = Stress_k2[2];
        kappa_k[0] = kappa_k2[0];
        kappa_k[1] = kappa_k2[1];
        Lambda_k = Lambda_k2;
        delta_lambda_k = delta_lambda_k2;
        Norm_Residual_1 = Norm_Residual_2;

        Convergence = check_convergence(Norm_Residual_1,Iter,MaxIter);
        
      }

    }

    /*
      Update equivalent plastic strain and increment of plastic deformation
    */
    Outputs_VarCons = fill_Outputs(Increment_E_plastic,Stress_k,Plastic_Flow_k,Lambda_k,delta_lambda_k,kappa_k[0]);

  }
  else
  {
    Outputs_VarCons = fill_Outputs(Increment_E_plastic,Stress_k,Plastic_Flow_k,Lambda_k,delta_lambda_k,kappa_k[0]);
  }

  
  return Outputs_VarCons;
}

/**************************************************************/

static Model_Parameters fill_model_paramters(Material MatProp)
{
  Model_Parameters Params;

  if(strcmp(MatProp.Yield_Function_Frictional,"Matsuoka-Nakai") == 0)
  {
    Params.m = 0.0;
    Params.c0 = 9;
    Is_Matsuoka_Nakai = true;
    Is_Lade_Duncan = false;
    Is_Modified_Lade_Duncan = false;
  }
  else if(strcmp(MatProp.Yield_Function_Frictional,"Lade-Duncan") == 0)
  {
    Params.m = 0.0;
    Params.c0 = 27;
    Is_Matsuoka_Nakai = false;
    Is_Lade_Duncan = true;
    Is_Modified_Lade_Duncan = false;
  }
  else if(strcmp(MatProp.Yield_Function_Frictional,"Modified-Lade-Duncan") == 0)
  {
    Params.m = MatProp.m_Frictional;
    Params.c0 = 27;
    Is_Matsuoka_Nakai = false;
    Is_Lade_Duncan = false;
    Is_Modified_Lade_Duncan = true;
  }

  Params.pa = MatProp.atmospheric_pressure;
  Params.alpha =  MatProp.alpha_Hardening_Borja;
  Params.a1 = MatProp.a_Hardening_Borja[0];
  Params.a2 = MatProp.a_Hardening_Borja[1];
  Params.a3 = MatProp.a_Hardening_Borja[2];

  double E = MatProp.E;
  double nu = MatProp.nu;

  Params.CC[0] =   1.0/E;
  Params.CC[1] = - nu/E;
  Params.CC[2] = - nu/E;

  Params.CC[3] = - nu/E;
  Params.CC[4] =   1.0/E;
  Params.CC[5] = - nu/E;

  Params.CC[6] = - nu/E;
  Params.CC[7] = - nu/E;
  Params.CC[8] =   1.0/E;

  return Params;
}

/**************************************************************/

static double eval_K1(
  double kappa1,
  double I1,
  double c0,
  double m,
  double pa)
{
  return c0 + kappa1*pow(pa/I1,m);
}

/**************************************************************/

static double eval_K2(
  double kappa2,
  double I1,
  double c0,
  double m,
  double pa)
{
  return c0 + kappa2*pow(pa/I1,m);
}

/**************************************************************/

static void eval_d_K2_d_stress(
  double * d_K2_d_stress,
  double kappa2,
  double I1,  
  double m,
  double pa)
{
  d_K2_d_stress[0] = - (m*kappa2/I1)*pow(pa/I1,m);
  d_K2_d_stress[1] = - (m*kappa2/I1)*pow(pa/I1,m);
  d_K2_d_stress[2] = - (m*kappa2/I1)*pow(pa/I1,m);
}

/**************************************************************/

static double eval_b1(
  double kappa1,
  double I1,
  double I3,  
  double m,
  double pa)
{
  return m*kappa1*(pow(pa/I1,m))*(cbrt(I3)/I1);
}

/**************************************************************/

static double eval_b2(
  double kappa2,
  double I1,
  double I3, 
  double m,
  double pa)
{
  return m*kappa2*(pow(pa/I1,m))*(cbrt(I3)/I1);
}

/**************************************************************/

static void eval_d_b2_d_stress(
  double * d_b2_d_stress,
  double * Stress,
  double kappa2,
  double I1,
  double I3, 
  double m,
  double pa)
  {

    double b2 = eval_b2(kappa2,I1,I3,m,pa);

    d_b2_d_stress[0] = (b2/I1)*(I1/(3*Stress[0]) - m - 1.0);
    d_b2_d_stress[1] = (b2/I1)*(I1/(3*Stress[1]) - m - 1.0);
    d_b2_d_stress[2] = (b2/I1)*(I1/(3*Stress[2]) - m - 1.0);
  } 
    
/**************************************************************/

static void eval_kappa(
  double * kappa,
  double Lambda,
  double I1,
  double a1,
  double a2,
  double a3,
  double alpha)
{
  kappa[0] = a1*Lambda*exp(a2*I1)*exp(-a3*Lambda);
  kappa[1] = alpha*kappa[0];
}

/**************************************************************/

static void eval_d_kappa1_d_stress(
  double * d_kappa1_d_stress,
  double Lambda,
  double I1,
  double a1,
  double a2,
  double a3)
{
  d_kappa1_d_stress[0] = a1*a2*Lambda*exp(a2*I1)*exp(-a3*Lambda);
  d_kappa1_d_stress[1] = a1*a2*Lambda*exp(a2*I1)*exp(-a3*Lambda);
  d_kappa1_d_stress[2] = a1*a2*Lambda*exp(a2*I1)*exp(-a3*Lambda);
}

/**************************************************************/

static double eval_d_kappa1_d_lambda(
  double Lambda,
  double I1,
  double a1,
  double a2,
  double a3)
{
  return (1 - a3*Lambda)*a1*exp(a2*I1)*exp(-a3*Lambda);
}
    
/**************************************************************/

static double eval_f(
  double I1,
  double I2)
  {

    if(Is_Matsuoka_Nakai)
    {
      return cbrt(I1*I2);
    }
    else if(Is_Lade_Duncan)
    {
      return I1;
    }
    else if(Is_Modified_Lade_Duncan)
    {
      return I1;
    }
    else
    {
      fprintf(stderr,"%s : %s !!! \n",
	    "Error in Frictional_Monolithic()",
	    "Undefined kind of function for f");
      exit(EXIT_FAILURE);    
    }  
    
  } 

/**************************************************************/

static void eval_d_f_d_stress(
  double * d_f_Matsuoka_Nakai_d_stress,
  double * Stress,
  double I1,
  double I2)
  {
    if(Is_Matsuoka_Nakai)
    {
      d_f_Matsuoka_Nakai_d_stress[0] = (I1*(I1 - Stress[0]) + I2)/(3*pow(cbrt(I1*I2),2));
      d_f_Matsuoka_Nakai_d_stress[1] = (I1*(I1 - Stress[1]) + I2)/(3*pow(cbrt(I1*I2),2));
      d_f_Matsuoka_Nakai_d_stress[2] = (I1*(I1 - Stress[2]) + I2)/(3*pow(cbrt(I1*I2),2));
    }
    else if(Is_Lade_Duncan)
    {
      d_f_Matsuoka_Nakai_d_stress[0] = 1.0;
      d_f_Matsuoka_Nakai_d_stress[1] = 1.0;
      d_f_Matsuoka_Nakai_d_stress[2] = 1.0;
    }
    else if(Is_Modified_Lade_Duncan)
    {
      d_f_Matsuoka_Nakai_d_stress[0] = 1.0;
      d_f_Matsuoka_Nakai_d_stress[1] = 1.0;
      d_f_Matsuoka_Nakai_d_stress[2] = 1.0;
    }
    else
    {
      fprintf(stderr,"%s : %s !!! \n",
	    "Error in Frictional_Monolithic()",
	    "Undefined kind of function for df");
      exit(EXIT_FAILURE);    
    } 
  } 

/**************************************************************/

static double eval_g(
  double I1,
  double I2)
  {
    if(Is_Matsuoka_Nakai)
    {
      return cbrt(I1*I2);
    }
    else if(Is_Lade_Duncan)
    {
      return I1;
    }
    else if(Is_Modified_Lade_Duncan)
    {
      return I1;
    }
    else
    {
      fprintf(stderr,"%s : %s !!! \n",
	    "Error in Frictional_Monolithic()",
	    "Undefined kind of function for g");
      exit(EXIT_FAILURE);    
    }  
  } 

/**************************************************************/

static void eval_d_g_d_stress(
  double * d_g_d_stress,
  double * Stress,
  double I1,
  double I2)
  {
    if(Is_Matsuoka_Nakai)
    {
      d_g_d_stress[0] = (I1*(I1 - Stress[0]) + I2)/(3*pow(cbrt(I1*I2),2));
      d_g_d_stress[1] = (I1*(I1 - Stress[1]) + I2)/(3*pow(cbrt(I1*I2),2));
      d_g_d_stress[2] = (I1*(I1 - Stress[2]) + I2)/(3*pow(cbrt(I1*I2),2));
    }
    else if(Is_Lade_Duncan)
    {
      d_g_d_stress[0] = 1.0;
      d_g_d_stress[1] = 1.0;
      d_g_d_stress[2] = 1.0;
    }
    else if(Is_Modified_Lade_Duncan)
    {
      d_g_d_stress[0] = 1.0;
      d_g_d_stress[1] = 1.0;
      d_g_d_stress[2] = 1.0;
    }
    else
    {
      fprintf(stderr,"%s : %s !!! \n",
	    "Error in Frictional_Monolithic()",
	    "Undefined kind of function for dg");
      exit(EXIT_FAILURE);    
    } 

  } 


/**************************************************************/

static void eval_dd_g_dd_stress(
  double * dd_g_dd_stress,
  double * Stress,
  double I1,
  double I2)
  {
    if(Is_Matsuoka_Nakai)
    {
      double dg_d_stress[3]; eval_d_g_d_stress(dg_d_stress,Stress,I1,I2);

      for (int A = 0; A<3; A++)
      {
        for (int B = 0; B<3; B++)
        {  
          dd_g_dd_stress[A*3 + B] = (3.0*I1 - Stress[A] - Stress[B] - I1*(A==B))/(3.0*pow(cbrt(I1*I2),2)) 
	    - (2.0/cbrt(I1*I2))*dg_d_stress[A]*dg_d_stress[B];
        }
      }
    }
    else if(Is_Lade_Duncan)
    {
      for (int A = 0; A<3; A++)
      {
        for (int B = 0; B<3; B++)
        {
          dd_g_dd_stress[A*3 + B] = 0.0;
        }
      }
    }
    else if(Is_Modified_Lade_Duncan)
    {
      for (int A = 0; A<3; A++)
      {
        for (int B = 0; B<3; B++)
        {
          dd_g_dd_stress[A*3 + B] = 0.0;
        }
      }
    }

  }


/**************************************************************/

static double eval_Yield_Function(
  double c0,
  double kappa1,
  double pa,
  double I1,
  double I2,
  double I3,
  double m)
{

  double K1 = eval_K1(kappa1,I1,c0,m,pa);

  double f = eval_f(I1,I2);

  return cbrt(K1*I3) - f;
}

/**************************************************************/

static void eval_d_Yield_Function_d_stress(
  double * d_F_d_stress,
  double * Stress,
  double I1,
  double I2,
  double I3,
  double c0,
  double kappa1,
  double pa,
  double m)
  {
    double Grad_f[3]; eval_d_f_d_stress(Grad_f,Stress,I1,I2);
    double K1 = eval_K1(kappa1,I1,c0,m,pa);
    double b1 = eval_b1(kappa1,I1,I3,m,pa);

    d_F_d_stress[0] = cbrt(K1*I3)/(3*Stress[0]) - b1/(3*pow(cbrt(K1),2)) - Grad_f[0];
    d_F_d_stress[1] = cbrt(K1*I3)/(3*Stress[1]) - b1/(3*pow(cbrt(K1),2)) - Grad_f[1];
    d_F_d_stress[2] = cbrt(K1*I3)/(3*Stress[2]) - b1/(3*pow(cbrt(K1),2)) - Grad_f[2];
  } 

/**************************************************************/

static double eval_d_Yield_Function_d_kappa1(
  double I1,
  double I3,
  double c0,
  double m,
  double pa,
  double kappa1)
  {
    double K1 = eval_K1(kappa1,I1,c0,m,pa);
    return cbrt(I3)/(3*pow(cbrt(K1),2))*pow(pa/I1,m);
  }

/**************************************************************/

static double eval_Plastic_Potential(
  double c0,
  double kappa2,
  double pa,
  double I1,
  double I2,
  double I3,
  double m)
  {
    double K2 = eval_K2(kappa2,I1,c0,m,pa);
    double g = eval_g(I1,I2);
    return cbrt(K2*I3) - g;
  } 

/**************************************************************/

static void eval_d_Plastic_Potential_d_stress(
  double * d_G_d_stress,
  double * Stress,
  double I1,
  double I2,
  double I3,
  double c0,
  double kappa2,
  double pa,
  double m)
  {
    double Grad_g[3]; eval_d_g_d_stress(Grad_g,Stress,I1,I2);
    double K2 = eval_K2(kappa2,I1,c0,m,pa);
    double b2 = eval_b2(kappa2,I1,I3,m,pa);

    d_G_d_stress[0] = cbrt(K2*I3)/(3*Stress[0]) - b2/(3*pow(cbrt(K2),2)) - Grad_g[0];
    d_G_d_stress[1] = cbrt(K2*I3)/(3*Stress[1]) - b2/(3*pow(cbrt(K2),2)) - Grad_g[1];
    d_G_d_stress[2] = cbrt(K2*I3)/(3*Stress[2]) - b2/(3*pow(cbrt(K2),2)) - Grad_g[2];
  } 

/**************************************************************/

static void eval_dd_Plastic_Potential_dd_stress(
  double * Hess_G,
  double * Stress,
  double kappa2,
  double I1,
  double I2,
  double I3,
  double m,
  double pa,
  double c0)
  {

    double K2 = eval_K2(kappa2,I1,c0,m,pa);
    double b2 = eval_b2(kappa2,I1,I3,m,pa);
    double Grad_K2[3]; eval_d_K2_d_stress(Grad_K2,kappa2,I1,m,pa);
    double Grad_b2[3]; eval_d_b2_d_stress(Grad_b2,Stress,kappa2,I1,I3,m,pa);
    double Hess_g[9]; eval_dd_g_dd_stress(Hess_g,Stress,I1,I2);

    for (int A = 0; A<3; A++)
    {
      for (int B = 0; B<3; B++)
      {
        Hess_G[A*3 + B] = (1.0/3.0)*cbrt(K2*I3)*(1.0/(3*Stress[A]*Stress[B]) - 1.0*(A==B)/pow(Stress[A],2))
        + (cbrt(I3)/Stress[A] + 2.0*b2/K2)*Grad_K2[B]/(9.0*pow(cbrt(K2),2))
        - Grad_b2[B]/(3.0*pow(cbrt(K2),2))
        - Hess_g[A*3 + B];
      }
    }

  }

/**************************************************************/

static void eval_dd_Plastic_Potential_d_stress_d_kappa2(
  double * dd_G_d_stress_d_kappa2, 
  double * Stress,
  double I1,
  double I3,
  double m,
  double pa,
  double c0,
  double kappa2)
{

  double K2 = eval_K2(kappa2,I1,c0,m,pa);
  double b2 = eval_b2(kappa2,I1,I3,m,pa);

  dd_G_d_stress_d_kappa2[0] = pow((pa/I1),m)*(cbrt(I3)/(3.0*Stress[0]) + 2.0*b2/(3.0*K2) - m*cbrt(I3)/I1)/(3.0*pow(cbrt(K2),2));
  dd_G_d_stress_d_kappa2[1] = pow((pa/I1),m)*(cbrt(I3)/(3.0*Stress[1]) + 2.0*b2/(3.0*K2) - m*cbrt(I3)/I1)/(3.0*pow(cbrt(K2),2));
  dd_G_d_stress_d_kappa2[2] = pow((pa/I1),m)*(cbrt(I3)/(3.0*Stress[2]) + 2.0*b2/(3.0*K2) - m*cbrt(I3)/I1)/(3.0*pow(cbrt(K2),2));

}

/**************************************************************/

static void eval_strain(
  double * Strain,
  double * Stress,
  double * CC)
{
  Strain[0] = CC[0]*Stress[0] + CC[1]*Stress[1] + CC[2]*Stress[2];
  Strain[1] = CC[3]*Stress[0] + CC[4]*Stress[1] + CC[5]*Stress[2];
  Strain[2] = CC[6]*Stress[0] + CC[7]*Stress[1] + CC[8]*Stress[2];
}

/**************************************************************/

static double assemble_residual(
  double * Residual,
  double * Stress_k,
  double * Strain_e_tri,
  double * Grad_G,
  double * kappa,
  double delta_lambda,
  double Lambda_k,
  Model_Parameters Params)
  {
    double alpha = Params.alpha;
    double m =  Params.m;
    double pa = Params.pa;
    double c0 = Params.c0;
    double a1 = Params.a1;
    double a2 = Params.a2;
    double a3 = Params.a3;
    double I1 = Stress_k[0] + Stress_k[1] + Stress_k[2];
    double I2 = Stress_k[0]*Stress_k[1] + Stress_k[1]*Stress_k[2] + Stress_k[0]*Stress_k[2];
    double I3 = Stress_k[0]*Stress_k[1]*Stress_k[2];
    double Strain_e_k[3];
    double kappa_hat[3];
    double F;
    double Error;
    
    if(fabs(I1/3.0) < TOL_Radial_Returning)
    {
      fprintf(stderr,"%s: %s %i %s \n",
      "Error in Frictional_Monolithic()",
      "The stress state of particle",Particle_Idx,"has reach the apex");
      exit(EXIT_FAILURE);
    }

    eval_kappa(kappa_hat,Lambda_k,I1,a1,a2,a3,alpha);
    eval_d_Plastic_Potential_d_stress(Grad_G,Stress_k,I1,I2,I3,c0,kappa[1],pa,m);
    F = eval_Yield_Function(c0,kappa[0],pa,I1,I2,I3,m);

    Residual[0] = Strain_e_k[0] - Strain_e_tri[0] + delta_lambda*Grad_G[0];
    Residual[1] = Strain_e_k[1] - Strain_e_tri[1] + delta_lambda*Grad_G[1];
    Residual[2] = Strain_e_k[2] - Strain_e_tri[2] + delta_lambda*Grad_G[2];
    Residual[3] = kappa[0] - kappa_hat[0];
    Residual[4] = F;
    
    /*
      Compute absolute error from the residual
    */
    for(int A = 0 ; A<5 ; A++)
    {
      Error += DSQR(Residual[A]);
    }

    Error = pow(Error,0.5);

    return Error;

  }

/**************************************************************/

static bool check_convergence(
  double Error,
  int Iter,
  int MaxIter)
{
  bool convergence;
  double Error_relative = 0.0;

  /*
    Compute relative error
  */
  if(Iter == 0)
  {
    Error0 = Error;
    Error_relative = Error/Error0;      

    if(Error0 < TOL_Radial_Returning)
    {
      return true;
    }
  }

 /*
    Compute residual error
   */
  Error_relative = Error/Error0;


  if((Error > TOL_Radial_Returning*100) 
  && (Error_relative > TOL_Radial_Returning) 
  && (Iter < MaxIter))
  {
    return false;
  }
  else
  {

    if(Iter >= MaxIter) 
    {
      fprintf(stderr,"%s: %s %i \n",
      "Error in Frictional_Monolithic()",
      "Maximm number of iteration reached in particle",Particle_Idx);
      fprintf(stderr,"Iter: %i | Error0: %e | Error: %e | Erro_rel: %f \n",Iter,Error0,Error,Error_relative);
    }

    return true;
  }

}

/**************************************************************/

static void assemble_tangent_matrix(
  double * Tangent_Matrix,
  double * d_G_d_stress,
  double * Stress_k, 
  double * kappa_k, 
  double delta_lambda,
  double Lambda_k, 
  Model_Parameters Params)
  {

    double alpha = Params.alpha;
    double m =  Params.m;
    double pa = Params.pa;
    double c0 = Params.c0;
    double a1 = Params.a1;
    double a2 = Params.a2;
    double a3 = Params.a3;
    double I1 = Stress_k[0] + Stress_k[1] + Stress_k[2];
    double I2 = Stress_k[0]*Stress_k[1] + Stress_k[1]*Stress_k[2] + Stress_k[0]*Stress_k[2];
    double I3 = Stress_k[0]*Stress_k[1]*Stress_k[2];
    double dd_G_dd_stress[9]; eval_dd_Plastic_Potential_dd_stress(dd_G_dd_stress,Stress_k,kappa_k[1],I1,I2,I3,m,pa,c0);
    double dd_G_d_stress_d_kappa2[3]; eval_dd_Plastic_Potential_d_stress_d_kappa2(dd_G_d_stress_d_kappa2,Stress_k,I1,I3,m,pa,c0,kappa_k[1]);
    double d_kappa1_d_stress[3]; eval_d_kappa1_d_stress(d_kappa1_d_stress,Lambda_k,I1,a1,a2,a3);
    double d_kappa1_d_lambda = eval_d_kappa1_d_lambda(Lambda_k,I1,a1,a2,a3);
    double d_F_d_stress[3]; eval_d_Yield_Function_d_stress(d_F_d_stress,Stress_k,I1,I2,I3,c0,kappa_k[0],pa,m);
    double d_F_d_kappa1 = eval_d_Yield_Function_d_kappa1(I1,I3,c0,m,pa,kappa_k[0]);

    /* First row */
    Tangent_Matrix[0] = Params.CC[0] + delta_lambda*dd_G_dd_stress[0];
    Tangent_Matrix[1] = Params.CC[1] + delta_lambda*dd_G_dd_stress[1];
    Tangent_Matrix[2] = Params.CC[2] + delta_lambda*dd_G_dd_stress[2];
    Tangent_Matrix[3] = alpha*delta_lambda*dd_G_d_stress_d_kappa2[0];
    Tangent_Matrix[4] = d_G_d_stress[0];

    /* Second row */
    Tangent_Matrix[5] = Params.CC[3] + delta_lambda*dd_G_dd_stress[3];
    Tangent_Matrix[6] = Params.CC[4] + delta_lambda*dd_G_dd_stress[4];
    Tangent_Matrix[7] = Params.CC[5] + delta_lambda*dd_G_dd_stress[5];
    Tangent_Matrix[8] = alpha*delta_lambda*dd_G_d_stress_d_kappa2[1];
    Tangent_Matrix[9] = d_G_d_stress[1];

    /* Third row */    
    Tangent_Matrix[10] = Params.CC[6] + delta_lambda*dd_G_dd_stress[6];
    Tangent_Matrix[11] = Params.CC[7] + delta_lambda*dd_G_dd_stress[7];
    Tangent_Matrix[12] = Params.CC[8] + delta_lambda*dd_G_dd_stress[8];
    Tangent_Matrix[13] = alpha*delta_lambda*dd_G_d_stress_d_kappa2[2];
    Tangent_Matrix[14] = d_G_d_stress[2];

    /* Four row */    
    Tangent_Matrix[15] = - d_kappa1_d_stress[0];
    Tangent_Matrix[16] = - d_kappa1_d_stress[1];
    Tangent_Matrix[17] = - d_kappa1_d_stress[2];
    Tangent_Matrix[18] = 1.0;
    Tangent_Matrix[19] = - d_kappa1_d_lambda;

    /* Five row */    
    Tangent_Matrix[20] = d_F_d_stress[0];
    Tangent_Matrix[21] = d_F_d_stress[1];
    Tangent_Matrix[22] = d_F_d_stress[2];
    Tangent_Matrix[23] = d_F_d_kappa1;
    Tangent_Matrix[24] = 0.0;

  } 

/**************************************************************/

static double compute_condition_number(double * Tangent_Matrix)
/*
  C = rcond(Tangent_Matrix) returns an estimate for the reciprocal condition of Tangent_Matrix in 1-norm. 
*/
{

  double RCOND; 
  double ANORM;
  int INFO;
  int N_rows = 5;
  int N_cols = 5;
  int LDA = IMAX(N_rows,N_cols);
  double * AUX_MEMORY = (double *)calloc(N_rows*N_cols,sizeof(double));
  double * WORK_ANORM = (double *)calloc(IMAX(1,N_rows),sizeof(double));
  double * WORK_RCOND = (double *)calloc(4*N_rows,sizeof(double));
  int * IPIV = (int *)calloc(IMIN(N_rows,N_cols),sizeof(int));
  int * IWORK_RCOND = (int *)calloc(N_rows,sizeof(int));

  // Copy matrix because dgetrf_ is a destructive operation
  memcpy(AUX_MEMORY, Tangent_Matrix, N_rows*N_cols*sizeof(double));

  // Compute 1-norm
  ANORM = dlange_("1", &N_rows, &N_cols, AUX_MEMORY, &LDA, WORK_ANORM);

  // The factors L and U from the factorization A = P*L*U
  dgetrf_(&N_rows,&N_cols,AUX_MEMORY,&LDA,IPIV,&INFO);

  // Check output of dgetrf
  if(INFO != 0)
  {
    if(INFO < 0)
    {
      printf("%s : \n","Error in Frictional_Monolithic()");
      printf("the %i-th argument had an illegal value",abs(INFO));
    }
    else if(INFO > 0)
    {
      printf("%s :\n","Error in Frictional_Monolithic()");
      printf(" M(%i,%i) %s \n %s \n %s \n %s \n",INFO,INFO,"is exactly zero. The factorization",
        "has been completed, but the factor M is exactly",
        "singular, and division by zero will occur if it is used",
        "to solve a system of equations.");
    
    }    
    exit(EXIT_FAILURE);
  }

  // Compute the Reciprocal condition number
  dgecon_("1", &N_rows, AUX_MEMORY, &LDA, &ANORM, &RCOND, WORK_RCOND, IWORK_RCOND, &INFO);

 if(INFO<0)
 {
  printf("Error in Frictional_Monolithic() : the %i-th argument of dgecon_ had an illegal value\n",abs(INFO));
  exit(EXIT_FAILURE);
 }

 // Free auxiliar memory
 free(IPIV);
 free(WORK_ANORM);
 free(WORK_RCOND);
 free(IWORK_RCOND);
 free(AUX_MEMORY);

 return RCOND;
}

/**************************************************************/
static void solver(
  double * Tangent_Matrix,
  double * Residual,
  double * D_Residual)
{
  int Order = 5;
  int LDA   = 5;
  int LDB   = 5;
  char  TRANS = 'T'; /* (Transpose) */
  int   INFO = 3;
  int * IPIV = (int *)malloc(Order*sizeof(int));
  int NRHS = 1;
  
  double rcond = compute_condition_number(Tangent_Matrix);
  
  if(rcond < 1E-10)
  {
    fprintf(stderr,"%s: %s %i, %s: %e\n",
    "Error in Frictional_Monolithic",
    "Tangent_Matrix is near to singular matrix for particle",
    Particle_Idx,"rcond",rcond);
  }

  /*
    Generate auxiliar copy of the mass matrix to avoid destructive operations
  */
  memcpy(D_Residual, Residual, 5*sizeof(double));

  /*
    Compute the LU factorization 
  */
  dgetrf_(&Order,&Order,Tangent_Matrix,&LDA,IPIV,&INFO);

  /*
    Check error messages in the LAPACK LU descomposition  
  */
  if(INFO)
  {
    fprintf(stderr,"%s : %s %s %s \n",
      "Error in Frictional_Monolithic",
      "The function","dgetrf_","returned an error message !!!" );
    exit(EXIT_FAILURE);
  }

  /*
    Solve
  */
  dgetrs_(&TRANS,&Order,&NRHS,Tangent_Matrix,&LDA,IPIV,D_Residual,&LDB,&INFO);
  free(IPIV);

  /*
    Check error messages in the LAPACK solver  
  */
  if(INFO)
  {
    fprintf(stderr,"%s : %s %s %s \n","Error in Frictional_Monolithic",
      "The function","dgetrs_","returned an error message !!!" );
    exit(EXIT_FAILURE);
  }
}


/**************************************************************/

static void update_variables(
  double * D_Residual,
  double * Stress_k,
  double * kappa_k,
  double * delta_lambda,
  double * lambda_k,
  double lambda_n,
  double alpha,
  double delta)
{

  Stress_k[0] -= delta*D_Residual[0];
  Stress_k[1] -= delta*D_Residual[1];
  Stress_k[2] -= delta*D_Residual[2];
  kappa_k[0] -= delta*D_Residual[3];
  kappa_k[1] = alpha*kappa_k[0];
  *delta_lambda -= delta*D_Residual[4];
  *lambda_k = *delta_lambda + lambda_n; 
}


/**************************************************************/

static State_Parameters fill_Outputs(
  double * Increment_E_plastic,
  double * Stress_k,
  double * Plastic_Flow, 
  double Lambda_k,
  double delta_lambda,
  double kappa_k1)
{
  State_Parameters Outputs_VarCons;

  Outputs_VarCons.Equiv_Plast_Str = Lambda_k;
  Outputs_VarCons.Kappa = kappa_k1;
  Outputs_VarCons.Stress = Stress_k;
  Outputs_VarCons.Increment_E_plastic = Increment_E_plastic;
  Outputs_VarCons.Increment_E_plastic[0] = delta_lambda*Plastic_Flow[0];
  Outputs_VarCons.Increment_E_plastic[1] = delta_lambda*Plastic_Flow[1];
  Outputs_VarCons.Increment_E_plastic[2] = delta_lambda*Plastic_Flow[2];

  return Outputs_VarCons;
}

/**************************************************************/

