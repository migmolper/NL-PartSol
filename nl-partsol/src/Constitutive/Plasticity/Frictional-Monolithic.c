
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
  double cohesion;
  double friction_angle;
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
static int eval_d_kappa1_d_stress(double *,double,double,double,double,double);
static double eval_d_kappa1_d_lambda(double,double,double,double,double);
static int eval_f(double *,double,double);
static int eval_d_f_d_stress(double *,double *,double,double);
static int eval_g(double *,double,double);
static int eval_d_g_d_stress(double *,double *,double,double);
static int eval_dd_g_dd_stress(double *,double *,double,double);
static int eval_Yield_Function(double *,double,double,double,double,double,double,double,double,double);
static int eval_d_Yield_Function_d_stress(double *,double *,double,double,double,double,double,double,double);
static double eval_d_Yield_Function_d_kappa1(double,double,double,double,double,double);
static int eval_Plastic_Potential(double*,double,double,double,double,double,double,double);
static int eval_d_Plastic_Potential_d_stress(double *,double *,double,double,double,double,double,double,double);
static int eval_dd_Plastic_Potential_dd_stress(double *,double *,double,double,double,double,double,double,double);
static void eval_dd_Plastic_Potential_d_stress_d_kappa2(double *, double *,double,double,double,double,double,double);
static void eval_strain(double *,double *,double *);
static int assemble_residual(double *,double *,double *,double *,double *,double *,double,double,Model_Parameters);
static int compute_condition_number(double *,double *);
static bool check_convergence(double,int,int);
static int assemble_tangent_matrix(double *,double *,double *,double *,double,double,Model_Parameters);
static int solver(double *,double *,double *);
static void update_variables(double *,double *,double *,double *,double *,double,double);
static void update_state_variables(State_Parameters *,double *,double *,double *,double,double,double);

/**************************************************************/

int Frictional_Monolithic__Constitutive__(
  State_Parameters * ptr_SP_p,
  Material MatProp)
/*	
	Monolithic algorithm for a smooth Mohr-Coulomb plastic criterium
*/
{
  int status = 0;

  /*
    Get material variables
  */
  Model_Parameters Params = fill_model_paramters(MatProp);

  /*
    Get the index of the particle
  */
  Particle_Idx = ptr_SP_p->Particle_Idx;

  /*
    Assign values from the inputs state parameters
  */
  double * Increment_E_plastic = ptr_SP_p->Increment_E_plastic;
  double * Stress_k  = ptr_SP_p->Stress;
  double Plastic_Flow_k[3] = {0.0, 0.0, 0.0};
  double kappa_k[2] = {ptr_SP_p->Kappa, Params.alpha*ptr_SP_p->Kappa};
  double Lambda_n = ptr_SP_p->Equiv_Plast_Str;
  double Lambda_k = Lambda_n;
  double delta_lambda_k = 0.0;

  /*
    Initialize Newton-Raphson solver
  */
  double F = 0;
  double Strain_e_tri[3];
  double Residual[5];
  double D_Residual[5];
  double Tangent_Matrix[25]; 
  double Norm_Residual;
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
  status = eval_Yield_Function(&F,Params.c0,kappa_k[0],Params.pa,I1,I2,I3,Params.m,Params.cohesion,Params.friction_angle);
  if(status)
  {
    fprintf(stderr,"%s %s \n%s %s\n",
    "Error in the function",__func__,
    "File",__FILE__);
    return EXIT_FAILURE;
  }

  if(F > 1E-8)
  {

    /*
      Compute elastic trial with the inverse elastic relation
    */
    eval_strain(Strain_e_tri,Stress_k,Params.CC);
    
    status = assemble_residual(&Norm_Residual,Residual,Stress_k,Strain_e_tri,Plastic_Flow_k,kappa_k,delta_lambda_k,Lambda_k,Params); 

    Convergence = check_convergence(Norm_Residual,Iter,MaxIter);

    /*
      Newton-Rapson
    */
    while(Convergence == false)
    {
      Iter++;
             
      status = assemble_tangent_matrix(Tangent_Matrix,Plastic_Flow_k,Stress_k,kappa_k,delta_lambda_k,Lambda_k,Params);
      if(status)
      {
        fprintf(stderr,"%s %s \n%s %s\n",
        "Error in the function",__func__,
        "File",__FILE__);
        return EXIT_FAILURE;
      }

      status = solver(Tangent_Matrix,Residual,D_Residual);
      if(status)
      {
        fprintf(stderr,"%s %s \n%s %s\n",
        "Error in the function",__func__,
        "File",__FILE__);
        return EXIT_FAILURE;
      }

      update_variables(D_Residual,Stress_k,kappa_k,&delta_lambda_k,&Lambda_k,Lambda_n,Params.alpha);

      status = assemble_residual(&Norm_Residual,Residual,Stress_k,Strain_e_tri,Plastic_Flow_k,kappa_k,delta_lambda_k,Lambda_k,Params); 

      Convergence = check_convergence(Norm_Residual,Iter,MaxIter);
    }

    /*
      Update equivalent plastic strain and increment of plastic deformation
    */
    update_state_variables(ptr_SP_p,Increment_E_plastic,Stress_k,Plastic_Flow_k,Lambda_k,delta_lambda_k,kappa_k[0]);

  }
  else
  {
    update_state_variables(ptr_SP_p,Increment_E_plastic,Stress_k,Plastic_Flow_k,Lambda_k,delta_lambda_k,kappa_k[0]);
  }

  
  return EXIT_SUCCESS;
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
  Params.cohesion = MatProp.yield_stress_0;
  Params.friction_angle = MatProp.phi_Frictional;

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
  if(m == 0)
  {
    return c0 + kappa1;
  }
  else
  {
    return c0 + kappa1*pow(pa/I1,m);
  }
  
}

/**************************************************************/

static double eval_K2(
  double kappa2,
  double I1,
  double c0,
  double m,
  double pa)
{
  if(m == 0)
  {
    return c0 + kappa2;
  }
  else
  {
    return c0 + kappa2*pow(pa/I1,m);
  }
}

/**************************************************************/

static void eval_d_K2_d_stress(
  double * d_K2_d_stress,
  double kappa2,
  double I1,  
  double m,
  double pa)
{
  if(m == 0)
  {
    d_K2_d_stress[0] = 0.0;
    d_K2_d_stress[1] = 0.0;
    d_K2_d_stress[2] = 0.0;
  }
  else
  {
    d_K2_d_stress[0] = - (m*kappa2/I1)*pow(pa/I1,m);
    d_K2_d_stress[1] = - (m*kappa2/I1)*pow(pa/I1,m);
    d_K2_d_stress[2] = - (m*kappa2/I1)*pow(pa/I1,m);
  }
}

/**************************************************************/

static double eval_b1(
  double kappa1,
  double I1,
  double I3,  
  double m,
  double pa)
{
  if(m == 0)
  {
    return 0.0;
  }
  else
  {
    return m*kappa1*(pow(pa/I1,m))*(cbrt(I3)/I1);
  }
}

/**************************************************************/

static double eval_b2(
  double kappa2,
  double I1,
  double I3, 
  double m,
  double pa)
{
  if(m == 0)
  {
    return 0.0;
  }
  else
  {
    return m*kappa2*(pow(pa/I1,m))*(cbrt(I3)/I1);
  }
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

    if(m == 0)
    {
      d_b2_d_stress[0] = 0.0;
      d_b2_d_stress[1] = 0.0;
      d_b2_d_stress[2] = 0.0;
    }
    else
    {
      double b2 = eval_b2(kappa2,I1,I3,m,pa);

      d_b2_d_stress[0] = (b2/I1)*(I1/(3*Stress[0]) - m - 1.0);
      d_b2_d_stress[1] = (b2/I1)*(I1/(3*Stress[1]) - m - 1.0);
      d_b2_d_stress[2] = (b2/I1)*(I1/(3*Stress[2]) - m - 1.0);
    }
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

static int eval_d_kappa1_d_stress(
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

  return EXIT_SUCCESS;
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

static int eval_f(
  double * f,
  double I1,
  double I2)
  {

    if(Is_Matsuoka_Nakai)
    {
      (*f) = cbrt(I1*I2);
    }
    else if(Is_Lade_Duncan)
    {
      (*f) = I1;
    }
    else if(Is_Modified_Lade_Duncan)
    {
      (*f) = I1;
    }
    else
    {
      fprintf(stderr,"%s %s: %s \n%s %s\n",
      "Error in the function",__func__,
      "Undefined kind of function for f",
      "File",__FILE__);
      return EXIT_FAILURE; 
    }  
    
    return EXIT_SUCCESS;
  } 

/**************************************************************/

static int eval_d_f_d_stress(
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
      fprintf(stderr,"%s %s: %s \n%s %s\n",
      "Error in the function",__func__,
      "Undefined kind of function for df",
      "File",__FILE__);
      return EXIT_FAILURE;
    } 

    return EXIT_SUCCESS;

  } 

/**************************************************************/

static int eval_g(
  double * g,
  double I1,
  double I2)
  {
    if(Is_Matsuoka_Nakai)
    {
      (*g) = cbrt(I1*I2);
    }
    else if(Is_Lade_Duncan)
    {
      (*g) = I1;
    }
    else if(Is_Modified_Lade_Duncan)
    {
      (*g) = I1;
    }
    else
    {
      fprintf(stderr,"%s %s: %s \n%s %s\n",
      "Error in the function",__func__,
      "Undefined kind of function for g",
      "File",__FILE__);
      return EXIT_FAILURE;  
    }  

    return EXIT_SUCCESS;
  } 

/**************************************************************/

static int eval_d_g_d_stress(
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
      fprintf(stderr,"%s %s: %s \n%s %s\n",
      "Error in the function",__func__,
      "Undefined kind of function for dg",
      "File",__FILE__);
      return EXIT_FAILURE;
    } 

    return EXIT_SUCCESS;
  } 


/**************************************************************/

static int eval_dd_g_dd_stress(
  double * dd_g_dd_stress,
  double * Stress,
  double I1,
  double I2)
  {
    int status = 0;

    if(Is_Matsuoka_Nakai)
    {

      double dg_d_stress[3]; 
      
      status = eval_d_g_d_stress(dg_d_stress,Stress,I1,I2);
      if(status)
      {
        fprintf(stderr,"%s %s \n%s %s\n",
        "Error in the function",__func__,
        "File",__FILE__);
        return EXIT_FAILURE;
      }

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

  
  return EXIT_SUCCESS;
}


/**************************************************************/

static int eval_Yield_Function(
  double * Yield_Function,
  double c0,
  double kappa1,
  double pa,
  double I1,
  double I2,
  double I3,
  double m,
  double cohesion,
  double friction_angle)
{
  int status = 0;
  double f = 0; 
  double K1 = 0;

  K1 = eval_K1(kappa1,I1,c0,m,pa);

  status = eval_f(&f,I1,I2);
  if(status)
  {
    fprintf(stderr,"%s %s \n%s %s\n",
    "Error in the function",__func__,
    "File",__FILE__);
    return EXIT_FAILURE;
  }

  (*Yield_Function) = cbrt(K1*I3) - f - 2*cohesion*cos(friction_angle);

  return EXIT_SUCCESS;
}

/**************************************************************/

static int eval_d_Yield_Function_d_stress(
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
    int status = 0;
    double Grad_f[3]; 
    
    double K1 = eval_K1(kappa1,I1,c0,m,pa);
    double b1 = eval_b1(kappa1,I1,I3,m,pa);

    status = eval_d_f_d_stress(Grad_f,Stress,I1,I2);
    if(status)
    {
      fprintf(stderr,"%s %s \n%s %s\n",
      "Error in the function",__func__,
      "File",__FILE__);
      return EXIT_FAILURE;
    }

    d_F_d_stress[0] = cbrt(K1*I3)/(3*Stress[0]) - b1/(3*pow(cbrt(K1),2)) - Grad_f[0];
    d_F_d_stress[1] = cbrt(K1*I3)/(3*Stress[1]) - b1/(3*pow(cbrt(K1),2)) - Grad_f[1];
    d_F_d_stress[2] = cbrt(K1*I3)/(3*Stress[2]) - b1/(3*pow(cbrt(K1),2)) - Grad_f[2];

    return EXIT_SUCCESS;
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

static int eval_Plastic_Potential(
  double * Plastic_Potential,
  double c0,
  double kappa2,
  double pa,
  double I1,
  double I2,
  double I3,
  double m)
  {
    int status = 0;
    double K2 = eval_K2(kappa2,I1,c0,m,pa);
    double g = 0;
    
    status = eval_g(&g,I1,I2);
    if(status)
    {
      fprintf(stderr,"%s %s \n%s %s\n",
      "Error in the function",__func__,
      "File",__FILE__);
      return EXIT_FAILURE;
    }

    (*Plastic_Potential) = cbrt(K2*I3) - g;
    
    return EXIT_SUCCESS;
  } 

/**************************************************************/

static int eval_d_Plastic_Potential_d_stress(
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
    int status = 0;
    double Grad_g[3]; 
    
    status = eval_d_g_d_stress(Grad_g,Stress,I1,I2);
    if(status)
    {
      fprintf(stderr,"%s %s \n%s %s\n",
      "Error in the function",__func__,
      "File",__FILE__);
      return EXIT_FAILURE;
    }    
    
    double K2 = eval_K2(kappa2,I1,c0,m,pa);
    double b2 = eval_b2(kappa2,I1,I3,m,pa);

    d_G_d_stress[0] = cbrt(K2*I3)/(3*Stress[0]) - b2/(3*pow(cbrt(K2),2)) - Grad_g[0];
    d_G_d_stress[1] = cbrt(K2*I3)/(3*Stress[1]) - b2/(3*pow(cbrt(K2),2)) - Grad_g[1];
    d_G_d_stress[2] = cbrt(K2*I3)/(3*Stress[2]) - b2/(3*pow(cbrt(K2),2)) - Grad_g[2];

    return EXIT_SUCCESS;
  } 

/**************************************************************/

static int eval_dd_Plastic_Potential_dd_stress(
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
    int status = 0;
    double K2 = eval_K2(kappa2,I1,c0,m,pa);
    double b2 = eval_b2(kappa2,I1,I3,m,pa);
    double Grad_K2[3]; eval_d_K2_d_stress(Grad_K2,kappa2,I1,m,pa);
    double Grad_b2[3]; eval_d_b2_d_stress(Grad_b2,Stress,kappa2,I1,I3,m,pa);
    double Hess_g[9]; 
    
    status = eval_dd_g_dd_stress(Hess_g,Stress,I1,I2);
    if(status)
    {
      fprintf(stderr,"%s %s \n%s %s\n",
      "Error in the function",__func__,
      "File",__FILE__);
      return EXIT_FAILURE;
    }

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

    return EXIT_SUCCESS;
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

static int assemble_residual(
  double * Norm_Residual,
  double * Residual,
  double * Stress_k,
  double * Strain_e_tri,
  double * Grad_G,
  double * kappa,
  double delta_lambda,
  double Lambda_k,
  Model_Parameters Params)
  {
    int status = 0;
    double alpha = Params.alpha;
    double m =  Params.m;
    double pa = Params.pa;
    double c0 = Params.c0;
    double a1 = Params.a1;
    double a2 = Params.a2;
    double a3 = Params.a3;
    double cohesion = Params.cohesion;
    double friction_angle = Params.friction_angle;
    double I1 = Stress_k[0] + Stress_k[1] + Stress_k[2];
    double I2 = Stress_k[0]*Stress_k[1] + Stress_k[1]*Stress_k[2] + Stress_k[0]*Stress_k[2];
    double I3 = Stress_k[0]*Stress_k[1]*Stress_k[2];
    double Strain_e_k[3];
    double kappa_hat[3];
    double F = 0;

    if((fabs(I1/3.0) < TOL_Radial_Returning) && (m != 0))
    {
      fprintf(stderr,"%s: %s %i %s \n",
      "Error in Frictional_Monolithic()",
      "The stress state of particle",Particle_Idx,"has reach the apex");
      exit(EXIT_FAILURE);
    }

    eval_strain(Strain_e_k,Stress_k,Params.CC);
    eval_kappa(kappa_hat,Lambda_k,I1,a1,a2,a3,alpha);
    
    status = eval_d_Plastic_Potential_d_stress(Grad_G,Stress_k,I1,I2,I3,c0,kappa[1],pa,m);
    if(status)
    {
      fprintf(stderr,"%s %s \n%s %s\n",
      "Error in the function",__func__,
      "File",__FILE__);
      return EXIT_FAILURE;
    }

    status = eval_Yield_Function(&F,c0,kappa[0],pa,I1,I2,I3,m,cohesion,friction_angle);
    if(status)
    {
      fprintf(stderr,"%s %s \n%s %s\n",
      "Error in the function",__func__,
      "File",__FILE__);
      return EXIT_FAILURE;
    }

    Residual[0] = Strain_e_k[0] - Strain_e_tri[0] + delta_lambda*Grad_G[0];
    Residual[1] = Strain_e_k[1] - Strain_e_tri[1] + delta_lambda*Grad_G[1];
    Residual[2] = Strain_e_k[2] - Strain_e_tri[2] + delta_lambda*Grad_G[2];
    Residual[3] = kappa[0] - kappa_hat[0];
    Residual[4] = F;

    /*
      Compute absolute error from the residual
    */
   (*Norm_Residual) = 0.0;
   
    for(int A = 0 ; A<5 ; A++)
    {
      (*Norm_Residual) += DSQR(Residual[A]);
    }

    (*Norm_Residual) = pow((*Norm_Residual),0.5);

    return EXIT_SUCCESS;
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

    if(Error0 < TOL_Radial_Returning*100)
    {
      return true;
    }
  }
  else
  {
    Error_relative = Error/Error0;
  }


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
      fprintf(stderr,"%s %s: %s %i \n",
      "Error in",__func__,
      "Maximm number of iteration reached in particle",Particle_Idx);
      fprintf(stderr,"Iter: %i | Error0: %e | Error: %e | Erro_rel: %e \n",Iter,Error0,Error,Error_relative);
    }

    return true;
  }

}

/**************************************************************/

static int assemble_tangent_matrix(
  double * Tangent_Matrix,
  double * d_G_d_stress,
  double * Stress_k, 
  double * kappa_k, 
  double delta_lambda,
  double Lambda_k, 
  Model_Parameters Params)
  {
    int status = 0;
    double rcond = 0;
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
    double dd_G_dd_stress[9]; 
    double dd_G_d_stress_d_kappa2[3] = {0.0, 0.0, 0.0}; 
    double d_kappa1_d_stress[3] = {0.0, 0.0, 0.0};
    double d_F_d_stress[3] = {0.0, 0.0, 0.0};

    eval_dd_Plastic_Potential_d_stress_d_kappa2(dd_G_d_stress_d_kappa2,Stress_k,I1,I3,m,pa,c0,kappa_k[1]);
    
    double d_kappa1_d_lambda = eval_d_kappa1_d_lambda(Lambda_k,I1,a1,a2,a3);
    
    double d_F_d_kappa1 = eval_d_Yield_Function_d_kappa1(I1,I3,c0,m,pa,kappa_k[0]);

    status = eval_dd_Plastic_Potential_dd_stress(dd_G_dd_stress,Stress_k,kappa_k[1],I1,I2,I3,m,pa,c0); 
    if(status)
    {
      fprintf(stderr,"%s %s \n%s %s\n",
      "Error in the function",__func__,
      "File",__FILE__);
      return EXIT_FAILURE;
    }

    status = eval_d_kappa1_d_stress(d_kappa1_d_stress,Lambda_k,I1,a1,a2,a3);
    if(status)
    {
      fprintf(stderr,"%s %s \n%s %s\n",
      "Error in the function",__func__,
      "File",__FILE__);
      return EXIT_FAILURE;
    }

    status = eval_d_Yield_Function_d_stress(d_F_d_stress,Stress_k,I1,I2,I3,c0,kappa_k[0],pa,m);
    if(status)
    {
      fprintf(stderr,"%s %s \n%s %s\n",
      "Error in the function",__func__,
      "File",__FILE__);
      return EXIT_FAILURE;
    }

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
    
    
    status = compute_condition_number(&rcond,Tangent_Matrix);
    if(status)
    {
      fprintf(stderr,"%s %s \n%s %s\n",
      "Error in the function",__func__,
      "File",__FILE__);
      return EXIT_FAILURE;
    }
    
    if(rcond < 1E-12)
    {
      fprintf(stderr,"%s %s: %s %i \n%s %s\n",
      "Error in the function",__func__,
      "Tangent_Matrix is near to singular matrix for particle",
      Particle_Idx,"File",__FILE__);
      
      fprintf(stderr,"Tangent Matrix: \n");
      for(int i = 0 ; i<5 ; i++)
      {
        for(int j = 0 ; j<5 ; j++)
        {
          printf("%e ",Tangent_Matrix[i*5+j]);
        }
        printf("\n");
      }

      printf("delta_lambda: %e, Lambda_k: %e\n",delta_lambda,Lambda_k);
      
      return EXIT_FAILURE;
    }


    return EXIT_SUCCESS;
  } 

/**************************************************************/

static int compute_condition_number(double * RCOND, double * Tangent_Matrix)
/*
  C = rcond(Tangent_Matrix) returns an estimate for the reciprocal condition of Tangent_Matrix in 1-norm. 
*/
{

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
  if(INFO<0)
  {
    free(IPIV); free(WORK_ANORM);
    free(WORK_RCOND); free(IWORK_RCOND); free(AUX_MEMORY);
    fprintf(stderr,"%s %s: %s %i %s \n%s %s\n",
    "Error in the function",__func__,
    "The",-INFO,"argument of dgetrf_ had an illegal value",
    "File",__FILE__);
  }
  else if(INFO > 0)
  {
    free(IPIV); free(WORK_ANORM);
    free(WORK_RCOND); free(IWORK_RCOND); free(AUX_MEMORY);
    fprintf(stderr,"%s %s: %s M[%i][%i] %s \n%s %s\n",
    "Error in the function",__func__,
    "The function dgetrf_ return",INFO,INFO,"is singular,",
    "File",__FILE__);
    return EXIT_FAILURE;
  }

  // Compute the Reciprocal condition number
  dgecon_("1", &N_rows, AUX_MEMORY, &LDA, &ANORM, RCOND, WORK_RCOND, IWORK_RCOND, &INFO);
  if(INFO < 0)
  {
    free(IPIV); free(WORK_ANORM);
    free(WORK_RCOND); free(IWORK_RCOND); free(AUX_MEMORY);
    fprintf(stderr,"%s %s: %s %i %s \n%s %s\n",
    "Error in the function",__func__,
    "The",-INFO,"argument of dgecon_ had an illegal value",
    "File",__FILE__);
  }
  
  // Free auxiliar memory
  free(IPIV); free(WORK_ANORM);
  free(WORK_RCOND); free(IWORK_RCOND); free(AUX_MEMORY);

 return EXIT_SUCCESS;
}

/**************************************************************/
static int solver(
  double * Tangent_Matrix,
  double * Residual,
  double * D_Residual)
{
  int status = 0;
  int Order = 5;
  int LDA   = 5;
  int LDB   = 5;
  char  TRANS = 'N';
  int   INFO = 3;
  int * IPIV = (int *)malloc(Order*sizeof(int));
  int NRHS = 1;

  /*
    Generate auxiliar copy of the mass matrix to avoid destructive operations
  */
  memcpy(D_Residual, Residual, 5*sizeof(double));

  /*
    Compute the LU factorization 
  */
  dgetrf_(&Order,&Order,Tangent_Matrix,&LDA,IPIV,&INFO);
  if(INFO<0)
  {
    free(IPIV);
    fprintf(stderr,"%s %s: %s %i %s \n%s %s\n",
    "Error in the function",__func__,
    "The",-INFO,"argument of dgetrf_ had an illegal value",
    "File",__FILE__);
    return EXIT_FAILURE;
  }
  else if(INFO>0)
  {
    free(IPIV);
    fprintf(stderr,"%s %s: %s M[%i][%i] %s \n%s %s\n",
    "Error in the function",__func__,
    "The function dgetrf_ return",INFO,INFO,"is singular,",
    "File",__FILE__);
    return EXIT_FAILURE;
  }

  /*
    Solve
  */
  dgetrs_(&TRANS,&Order,&NRHS,Tangent_Matrix,&LDA,IPIV,D_Residual,&LDB,&INFO);
  if(INFO<0)
  {
    free(IPIV);
    fprintf(stderr,"%s %s: %s %i %s \n%s %s\n",
    "Error in the function",__func__,
    "The",-INFO,"argument of dgetrs_ had an illegal value",
    "File",__FILE__);
    return EXIT_FAILURE;
  }

  free(IPIV);

  return EXIT_SUCCESS;
}


/**************************************************************/

static void update_variables(
  double * D_Residual,
  double * Stress_k,
  double * kappa_k,
  double * delta_lambda,
  double * lambda_k,
  double lambda_n,
  double alpha)
{

  Stress_k[0] -= D_Residual[0];
  Stress_k[1] -= D_Residual[1];
  Stress_k[2] -= D_Residual[2];
  kappa_k[0] -= D_Residual[3];
  kappa_k[1] = alpha*kappa_k[0];
  *delta_lambda -= D_Residual[4];
  *lambda_k = *delta_lambda + lambda_n; 
}


/**************************************************************/

static void update_state_variables(
  State_Parameters * ptr_SP_p,
  double * Increment_E_plastic,
  double * Stress_k,
  double * Plastic_Flow, 
  double Lambda_k,
  double delta_lambda,
  double kappa_k1)
{

  ptr_SP_p->Equiv_Plast_Str = Lambda_k;
  ptr_SP_p->Kappa = kappa_k1;
  ptr_SP_p->Stress = Stress_k;
  ptr_SP_p->Increment_E_plastic = Increment_E_plastic;
  ptr_SP_p->Increment_E_plastic[0] = delta_lambda*Plastic_Flow[0];
  ptr_SP_p->Increment_E_plastic[1] = delta_lambda*Plastic_Flow[1];
  ptr_SP_p->Increment_E_plastic[2] = delta_lambda*Plastic_Flow[2];

}

/**************************************************************/

