#include "nl-partsol.h"

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
  double CC[3][3];

} Model_Parameters;

/*
  Auxiliar functions 
*/
static double eval_K1(double,double,double,double,double);
static double eval_K2(double,double,double,double,double);
static void eval_d_K2_d_stress(double *,double,double,double,double);
static double eval_b1(double,double,double,double,double);
static double eval_b2(double,double,double,double,double);
static void eval_d_b2_d_stress(double *,double *,double,double,double,double,double);
static void eval_kappa(double *,double,double,double,double,double,double);
static void eval_d_kappa1_d_stress(double *,double,double,double,double,double);
static double eval_d_kappa1_d_lambda(double,double,double,double,double);
static double eval_f_Matsuoka_Nakai(double,double);
static void eval_d_f_Matsuoka_Nakai_d_stress(double *,double *,double,double);
static double eval_g_Matsuoka_Nakai(double,double);
static void eval_d_g_Matsuoka_Nakai_d_stress(double *,double *,double,double);
static void eval_dd_g_Matsuoka_Nakai_dd_stress(double **,double *,double,double);
static double eval_F(double,double,double,double,double,double,double);
static void eval_d_F_d_stress(double *,double *,double,double,double,double,double,double,double);
static double eval_d_F_d_kappa1(double,double,double,double,double,double);
static double eval_G(double,double,double,double,double,double,double);
static void eval_d_G_d_stress(double *,double *,double,double,double,double,double,double,double);
static void eval_dd_G_dd_stress(double **,double *,double,double,double,double,double,double,double);
static void eval_dd_G_d_stress_d_kappa2(double *, double *,double,double,double,double,double,double);
static void eval_strain(double *,double *,double **);
static void assemble_residual(double *,double *,double *,double *,double,double,Model_Parameters);
static bool check_convergence(double *,double,int,int,int);
static void assemble_tangent_matrix(double *,double *,double *,double,double,Model_Parameters);
static void update_variables(double *,double *,double *,double *,double *,double *,double,double);

/**************************************************************/

State_Parameters Smooth_Mohr_Coulomb_Monolithic(
  State_Parameters Inputs_SP,
  Material MatProp)
/*	
	Radial returning algorithm for the Von-Mises plastic criterium
*/
{

  /*
    Define output state parameters
  */
  State_Parameters Outputs_VarCons;

  /*
    Initialise solver parameters
  */
  double EPS_k = Inputs_SP.EPS;
  double delta_Gamma_k = 0;
  double TOL = TOL_Radial_Returning;
  int MaxIter = Max_Iterations_Radial_Returning;
  int Iter = 0;
  bool Convergence = false;

  /*
    Get material variables
  */
  Model_Parameters Params = fill_model_paramters(Material MatProp);

  /*
    Compute invariants
  */
  double I1 = Stress_k[0] + Stress_k[1] + Stress_k[2];
  double I2 = Stress_k[0]*Stress_k[1] + Stress_k[1]*Stress_k[2] + Stress_k[0]*Stress_k[2];
  double I3 = Stress_k[0]*Stress_k[1]*Stress_k[2];

  /*
    Check yield condition
  */
  if(eval_F(c0,kappa1,pa,I1,I2,I3,m) > 0.0)
  {
    
    /*
      Newton-Rapson solver
    */
    while(Convergence == false)
    {
      
      assemble_residual(Residual,Stress_k,Strain_e_tri,kappa,delta_lambda,Lambda_k,Params);

      Convergence = check_convergence(residual,Error_0,TOL,iter,MaxIter)

      Convergence = check_convergence(Residual,TOL,Iter,MaxIter);

      if(Convergence == false)
      {
        
        Iter++;
        
        assemble_tangent_matrix(Tangent_Matrix,Stress_k,kappa_k,delta_lambda,Lambda_k,Params);
        
        update_variables(Tangent_Matrix,Residual,Stress_k,kappa_k,&delta_lambda,&lambda_k,lambda_n,alpha);

      }
    
    }

    /*
      Update plastic deformation gradient
    */
    update_increment_plastic_strain_tensor(Inputs_SP.Increment_E_plastic,plastic_flow_direction, delta_Gamma_k);

    /*
      Update equivalent plastic strain and increment of plastic deformation
    */
    Outputs_VarCons.EPS = EPS_k;
    Outputs_VarCons.Stress = Inputs_SP.Stress;
    Outputs_VarCons.Increment_E_plastic = Inputs_SP.Increment_E_plastic;
  }
  else
  {
    Outputs_VarCons.EPS = Inputs_SP.EPS;
    Outputs_VarCons.Stress = Inputs_SP.Stress;
    Outputs_VarCons.Increment_E_plastic = Inputs_SP.Increment_E_plastic;
  }
  
  return Outputs_VarCons;
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

    d_b2_d_stress[0] = (b2/I1)*(I1/(3*stress[0]) - m - 1.0);
    d_b2_d_stress[1] = (b2/I1)*(I1/(3*stress[1]) - m - 1.0);
    d_b2_d_stress[2] = (b2/I1)*(I1/(3*stress[2]) - m - 1.0);
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

static double eval_f_Matsuoka_Nakai(
  double I1,
  double I2)
  {
    return cbrt(I1*I2);
  } 

/**************************************************************/

static void eval_d_f_Matsuoka_Nakai_d_stress(
  double * d_f_Matsuoka_Nakai_d_stress,
  double * Stress,
  double I1,
  double I2)
  {
    d_f_Matsuoka_Nakai_d_stress[0] = (I1*(I1 - Stress[0]) + I2)/(3*pow(cbrt(I1*I2),2));
    d_f_Matsuoka_Nakai_d_stress[1] = (I1*(I1 - Stress[1]) + I2)/(3*pow(cbrt(I1*I2),2));
    d_f_Matsuoka_Nakai_d_stress[2] = (I1*(I1 - Stress[2]) + I2)/(3*pow(cbrt(I1*I2),2));
  } 

/**************************************************************/

static double eval_g_Matsuoka_Nakai(
  double I1,
  double I2)
  {
    return cbrt(I1*I2)**2;
  } 

/**************************************************************/

static void eval_d_g_Matsuoka_Nakai_d_stress(
  double * d_g_d_stress,
  double * Stress,
  double I1,
  double I2)
  {
    d_g_d_stress[0] = (I1*(I1 - Stress[0]) + I2)/(3*pow(cbrt(I1*I2),2));
    d_g_d_stress[1] = (I1*(I1 - Stress[1]) + I2)/(3*pow(cbrt(I1*I2),2));
    d_g_d_stress[2] = (I1*(I1 - Stress[2]) + I2)/(3*pow(cbrt(I1*I2),2));
  } 


/**************************************************************/

static void eval_dd_g_Matsuoka_Nakai_dd_stress(
  double ** dd_g_dd_stress,
  double * Stress,
  double I1,
  double I2)
  {
    for (int A = 0; A<3; A++)
    {
      for (int B = 0; B<3; B++)
      {
        dd_g_dd_stress[A][B] = (3.0*I1 - Stress[A] - Stress[B] - I1*(A==B))/(3.0*pow(cbrt(I1*I2),2)) 
        - (2.0/9.0)*(I1*(I1 - Stress[A]) + I2)*(I1*(I1 - Stress[B]) + I2)/(pow(cbrt(I1*I2),5));
      }
    }
  }


/**************************************************************/

static double eval_F(
  double c0,
  double kappa1,
  double pa,
  double I1,
  double I2,
  double I3,
  double m)
{

  double K1 = eval_K1(c0,kappa1,pa,I1,m);
  double f = eval_f_Matsuoka_Nakai(I1,I2);

  return cbrt(K1*I3) - f;
}

/**************************************************************/

static void eval_d_F_d_stress(
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
    double Grad_f[3], eval_d_f_Matsuoka_Nakai_d_stress(Grad_f,Stress,I1,I2);
    double K1 = eval_K1(kappa1,I1,c0,m,pa);
    double b1 = eval_b1(kappa1,I1,I3,m,pa);

    d_F_d_stress[0] = cbrt(K1*I3)/(3*Stress[0]) - b1/(3*pow(cbrt(K1),2)) - Grad_f[0];
    d_F_d_stress[0] = cbrt(K1*I3)/(3*Stress[0]) - b1/(3*pow(cbrt(K1),2)) - Grad_f[0];
    d_F_d_stress[0] = cbrt(K1*I3)/(3*Stress[0]) - b1/(3*pow(cbrt(K1),2)) - Grad_f[0];
  } 

/**************************************************************/

static double eval_d_F_d_kappa1(
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

static double eval_G(
  double c0,
  double kappa2,
  double pa,
  double I1,
  double I2,
  double I3,
  double m)
  {
    double K2 = eval_K2(kappa2,I1,c0,m,pa);
    double g = eval_g_Matsuoka_Nakai(I1,I2);
    return cbrt(K2*I3) - g;
  } 

/**************************************************************/

static void eval_d_G_d_stress(
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
    double Grad_g[3], eval_d_g_Matsuoka_Nakai_d_stress(Grad_g,Stress,I1,I2);
    double K2 = eval_K2(kappa2,I1,c0,m,pa);
    double b2 = eval_b2(kappa2,I1,I3,m,pa);

    d_G_d_stress[0] = cbrt(K2*I3)/(3*Stress[0]) - b2/(3*pow(cbrt(K2),2)) - Grad_g[0];
    d_G_d_stress[1] = cbrt(K2*I3)/(3*Stress[1]) - b2/(3*pow(cbrt(K2),2)) - Grad_g[1];
    d_G_d_stress[2] = cbrt(K2*I3)/(3*Stress[2]) - b2/(3*pow(cbrt(K2),2)) - Grad_g[2];
  } 

/**************************************************************/

static void eval_dd_G_dd_stress(
  double ** Hess_G,
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
    double Grad_K2[3], eval_d_K2_d_stress(Grad_K2,kappa2,I1,m,pa);
    double Grad_b2[3], eval_d_b2_d_stress(Grad_b2,Stress,kappa2,I1,I3,m,pa);
    double Hess_g[3][3], eval_dd_g_Matsuoka_Nakai_dd_stress(Hess_g,Stress,I1,I2);

    for (int A = 0; A<3; A++)
    {
      for (int B = 0; B<3; B++)
      {
        Hess_G[A][B] = (1.0/3.0)*cbrt(K2*I3)*(1.0/(3*Stress[A]*Stress[B]) - 1.0*(A==B)/pow(Stress[A],2)) 
        + (cbrt(I3)/Stress[A] + 2.0*b2/K2)*Grad_K2[B]/(9.0*pow(cbrt(K2),2)) 
        - Grad_b2[B]/(3.0*pow(cbrt(K2),2)) 
        - Hess_g[A][B];
      }
    }

  }

/**************************************************************/

static void eval_dd_G_d_stress_d_kappa2(
  double * d_G_d_stress_d_kappa2, 
  double * Stress,
  double I1,
  double I3,
  double m,
  double pa,
  double c0,
  double kappa2)
{

  double K2 = eval_K2(c0,kappa2,pa,I1,m);
  double b2 = eval_b2(kappa2,I1,I3,m,pa);

  d_G_d_stress_d_kappa2[0] = pow((pa/I1),m)*(cbrt(I3)/(3.0*Stress[0]) 
  + 2.0*b2/(3.0*K2) 
  - m*cbrt(I3)/I1)/(3.0*pow(cbrt(K2),2));
  d_G_d_stress_d_kappa2[1] = pow((pa/I1),m)*(cbrt(I3)/(3.0*Stress[1]) 
  + 2.0*b2/(3.0*K2) 
  - m*cbrt(I3)/I1)/(3.0*pow(cbrt(K2),2));
  d_G_d_stress_d_kappa2[2] = pow((pa/I1),m)*(cbrt(I3)/(3.0*Stress[2]) 
  + 2.0*b2/(3.0*K2) 
  - m*cbrt(I3)/I1)/(3.0*pow(cbrt(K2),2));

}

/**************************************************************/

static void eval_strain(
  double * Strain,
  double * Stress,
  double ** CC)
{
  Strain[0] = CC[0][0]*Stress[0] + CC[0][1]*Stress[1] + CC[0][2]*Stress[2]);
  Strain[1] = CC[1][0]*Stress[0] + CC[1][1]*Stress[1] + CC[1][2]*Stress[2]);
  Strain[2] = CC[2][0]*Stress[0] + CC[2][1]*Stress[1] + CC[2][2]*Stress[2]);
}

/**************************************************************/

static void assemble_residual(
  double * Residual,
  double * Stress_k,
  double * Strain_e_tri,
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
    double ** CC = Params.CC;
    double I1 = Stress_k[0] + Stress_k[1] + Stress_k[2];
    double I2 = Stress_k[0]*Stress_k[1] + Stress_k[1]*Stress_k[2] + Stress_k[0]*Stress_k[2];
    double I3 = Stress_k[0]*Stress_k[1]*Stress_k[2];
    double Strain_e_k[3], eval_strain(Strain_e_k,Stress_k,CC);
    double kappa_hat[3], eval_kappa(kappa_hat,Lambda_k,I1,a1,a2,a3,alpha);
    double Grad_G[3], eval_d_G_d_stress(Grad_G,Stress_k,I1,I2,I3,c0,kappa[1],pa,m);

    double F = eval_F(c0,kappa[0],pa,I1,I2,I3,m);

    Residual[0] = Strain_e_k[0] - Strain_e_tri[0] + delta_lambda*Grad_G[0];
    Residual[1] = Strain_e_k[1] - Strain_e_tri[1] + delta_lambda*Grad_G[1];
    Residual[2] = Strain_e_k[2] - Strain_e_tri[2] + delta_lambda*Grad_G[2];
    Residual[3] = kappa[0] - kappa_hat[0];
    Residual[4] = F;
  }

/**************************************************************/

static bool check_convergence(
  double * Residual,
  double TOL,
  int Iter,
  int MaxIter,
  int Step)
{
  bool convergence;
  int Ndim = NumberDimensions;
  int Ndof = NumberDOF;
  int Nnodes_mask = Residual.N_cols;
  double Error = 0.0;
  double Error_relative = 0.0;

  
  /*
    Compute absolute error 
  */
  for(int A = 0 ; A<5 ; A++)
  {
    Error += DSQR(Residual[A]);
  }

  Error = pow(Error,0.5);


  /*
    Compute relative error
  */
  if(Iter == 0)
  {
    Error0 = Error;
    Error_relative = Error/Error0;      
  }
  else if(Iter > MaxIter)
  {
    Error_relative = Error/Error0;
    fprintf(stderr,"\t Convergence no reached after : %i interations \n",Iter);
    fprintf(stderr,"\t Initial error : %1.4e \n",Error0);
    fprintf(stderr,"\t Total error : %1.4e \n",Error);
    fprintf(stderr,"\t Relative error : %1.4e \n",Error_relative);
    return true;
  }
  else
  {
      Error_relative = Error/Error0;
  }
      
  /*
    Check convergence using the relative error
  */
  if(Error_relative > TOL)
  {
    return false;
  }
  else
  {
    print_convergence_stats(Step, Iter, Error0, Error, Error_relative);
    return true;
  }

}

/**************************************************************/

static void assemble_tangent_matrix(
  double * Tangent_Matrix,
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
    double ** CC = Params.CC;
    double I1 = Stress_k[0] + Stress_k[1] + Stress_k[2];
    double I2 = Stress_k[0]*Stress_k[1] + Stress_k[1]*Stress_k[2] + Stress_k[0]*Stress_k[2];
    double I3 = Stress_k[0]*Stress_k[1]*Stress_k[2];
    double d_G_d_stress[3], eval_d_G_d_stress(d_G_d_stress,Stress_k,I1,I2,I3,c0,kappa_k[1],pa,m);
    double dd_G_dd_stress[3][3], eval_dd_G_dd_stress(dd_G_dd_stress,Stress_k,kappa_k[1],I1,I2,I3,m,pa,c0);
    double dd_G_d_stress_d_kappa2[3], eval_dd_G_d_stress_d_kappa2(d_G_d_stress_d_kappa2,Stress_k,I1,I3,m,pa,c0,kappa_k[1]);
    double d_kappa1_d_stress[3], eval_d_kappa1_d_stress(d_kappa1_d_stress,Lambda_k,I1,a1,a2,a3);
    double d_kappa1_d_lambda = eval_d_kappa1_d_lambda(Lambda_k,I1,a1,a2,a3);
    double d_F_d_stress[3], eval_d_F_d_stress(d_F_d_stress,Stress_k,I1,I2,I3,c0,kappa_k[0],pa,m);
    double d_F_d_kappa1 = eval_d_F_d_kappa1(I1,I3,c0,m,pa,kappa_k[0];
    
    /* First row */
    Tangent_Matrix[0] = CC[0][0] + delta_lambda*dd_G_dd_stress[0][0];
    Tangent_Matrix[1] = CC[0][1] + delta_lambda*dd_G_dd_stress[0][1];
    Tangent_Matrix[2] = CC[0][2] + delta_lambda*dd_G_dd_stress[0][2];
    Tangent_Matrix[3] = alpha*delta_lambda*dd_G_d_stress_d_kappa2[0];
    Tangent_Matrix[4] = d_G_d_stress[0];

    /* Second row */
    Tangent_Matrix[5] = CC[1][0] + delta_lambda*dd_G_dd_stress[1][0];
    Tangent_Matrix[6] = CC[1][1] + delta_lambda*dd_G_dd_stress[1][1];
    Tangent_Matrix[7] = CC[1][2] + delta_lambda*dd_G_dd_stress[1][2];
    Tangent_Matrix[8] = alpha*delta_lambda*dd_G_d_stress_d_kappa2[1];
    Tangent_Matrix[9] = d_G_d_stress[1];

    /* Third row */    
    Tangent_Matrix[10] = CC[2][0] + delta_lambda*dd_G_dd_stress[2][0];
    Tangent_Matrix[11] = CC[2][1] + delta_lambda*dd_G_dd_stress[2][1];
    Tangent_Matrix[12] = CC[2][2] + delta_lambda*dd_G_dd_stress[2][2];
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

static void update_variables(
  double * Tangent_Matrix,
  double * Residual
  double * Stress_k,
  double * kappa_k,
  double * delta_lambda,
  double * lambda_k,
  double lambda_n,
  double alpha)
{
  int Order = 5;
  int LDA   = 5;
  int LDB   = 5;
  char  TRANS = 'T'; /* (Transpose) */
  int   INFO = 3;
  int * IPIV = (int *)Allocate_Array(Order,sizeof(int));
  int NRHS = 1;

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
      "Error in solve_system",
      "The function","dgetrf_","returned an error message !!!" );
    exit(EXIT_FAILURE);
  }

  /*
    Solve
  */
  dgetrs_(&TRANS,&Order,&NRHS,Tangent_Matrix,&LDA,IPIV,Residual,&LDB,&INFO);
  free(IPIV);

  /*
    Check error messages in the LAPACK solver  
  */
  if(INFO)
  {
    fprintf(stderr,"%s : %s %s %s \n","Error in solve_system",
      "The function","dgetrs_","returned an error message !!!" );
    exit(EXIT_FAILURE);
  }

  /*
    Update 
  */
  Stress_k[0] -= Residual[0];
  Stress_k[1] -= Residual[1];
  Stress_k[2] -= Residual[2];
  kappa_k[0] -= Residual[3];
  kappa_k[1] = alpha*kappa_k[0];
  *delta_lambda -= Residual[4];
  *lambda_k = *delta_lambda + lambda_n;
}

/**************************************************************/