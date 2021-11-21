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

#define TOL_Error_0 1E-7
#define TOL_Error_i 1E-10
#define MAX_Iter 10

/*
  Auxiliar variables
*/
typedef struct
{
  double friction_angle;
  double dilatancy_angle;
  double CC[9];

} Model_Parameters;

/*
  Define local global variables
*/
int Particle_Idx;
double Error0;

/*
  Auxiliar functions 
*/
static Model_Parameters fill_model_paramters(Material Material);
static int eval_K_MN(double *,const double);
static int eval_Yield_Function(double *,const double,const double,const double,const double);
static int eval_d_Yield_Function_d_stress(double *,double *,const double,const double,const double,const double);
static int eval_Plastic_Potential(double*,const double,const double,const double,const double);
static int eval_d_Plastic_Potential_d_stress(double *, double *,const double,const double,const double,const double);
static int eval_dd_Plastic_Potential_dd_stress(double *,double *,const double,const double,const double,const double);
static void eval_strain(double *,double *,double *);
static int assemble_residual(double *,double *,double *,double *,double *,double,double,Model_Parameters);
static int compute_condition_number(double *,double *);
static bool check_convergence(double,int);
static int assemble_tangent_matrix(double *,double *,double *,double,double,double,Model_Parameters);
static int solver(double *,double *);
static int update_variables(double *,double *,double *,double *,double);
static int update_state_variables(State_Parameters *,double *,double *,double,double);

/**************************************************************/

int Matsuoka_Nakai__Constitutive__(
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
  double * Stress_k  = ptr_SP_p->Stress;
  double Plastic_Flow_k[3] = {0.0, 0.0, 0.0};
  double Lambda_n = ptr_SP_p->Equiv_Plast_Str;
  double Lambda_k = Lambda_n;
  double delta_lambda_k = 0.0;

  /*
    Initialize Newton-Raphson solver
  */
  double F = 0;
  double Strain_e_tri[3] = {0.0,0.0,0.0};
  double Residual[4]= {0.0,0.0,0.0,0.0};
  double Tangent_Matrix[16]; 
  double Norm_Residual = 1;
  int Iter = 0;
  bool Convergence = false;

  /*
    Compute invariants
  */
  double I1 = Stress_k[0] + Stress_k[1] + Stress_k[2];
  double I2 = Stress_k[0]*Stress_k[1] + Stress_k[1]*Stress_k[2] + Stress_k[0]*Stress_k[2];
  double I3 = Stress_k[0]*Stress_k[1]*Stress_k[2];

  /* Check yield condition */
  status = eval_Yield_Function(&F,I1,I2,I3,Params.friction_angle);
  if(status)
  {
    fprintf(stderr,"%s %s \n%s %s\n",
    "Error in the function",__func__,
    "File",__FILE__);
    return EXIT_FAILURE;
  }

  if(F > TOL_Error_i)
  {

    /* Compute elastic trial with the inverse elastic relation */
    eval_strain(Strain_e_tri,Stress_k,Params.CC);
    
    status = assemble_residual(&Norm_Residual,Residual,Stress_k,Strain_e_tri,Plastic_Flow_k,delta_lambda_k,Lambda_k,Params); 

    if(status)
    {
      fprintf(stderr,"%s %s \n%s %s\n",
      "Error in the function",__func__,
      "File",__FILE__);
      return EXIT_FAILURE;
    }

    Convergence = check_convergence(Norm_Residual,Iter);

    /* Newton-Rapson */
    while(Convergence == false)
    {
      Iter++;
             
      status = assemble_tangent_matrix(Tangent_Matrix,Plastic_Flow_k,Stress_k,delta_lambda_k,Lambda_k,Norm_Residual,Params);
      if(status)
      {
        fprintf(stderr,"%s %s \n%s %s\n",
        "Error in the function",__func__,
        "File",__FILE__);
        return EXIT_FAILURE;
      }

      status = solver(Tangent_Matrix,Residual);
      if(status)
      {
        fprintf(stderr,"%s %s \n%s %s\n",
        "Error in the function",__func__,
        "File",__FILE__);
        return EXIT_FAILURE;
      }

      status = update_variables(Residual,Stress_k,&delta_lambda_k,&Lambda_k,Lambda_n);
      if(status)
      {
        fprintf(stderr,"%s %s \n%s %s\n",
        "Error in the function",__func__,
        "File",__FILE__);
        return EXIT_FAILURE;
      }

      status = assemble_residual(&Norm_Residual,Residual,Stress_k,Strain_e_tri,Plastic_Flow_k,delta_lambda_k,Lambda_k,Params); 
      if(status)
      {
        fprintf(stderr,"%s %s \n%s %s\n",
        "Error in the function",__func__,
        "File",__FILE__);
        return EXIT_FAILURE;
      }

      Convergence = check_convergence(Norm_Residual,Iter);
    }

    /* Update equivalent plastic strain and increment of plastic deformation */
    status = update_state_variables(ptr_SP_p,Stress_k,Plastic_Flow_k,Lambda_k,delta_lambda_k);
    if(status)
    {
      fprintf(stderr,"%s %s \n%s %s\n",
      "Error in the function",__func__,
      "File",__FILE__);
      return EXIT_FAILURE;
    }

  }
  else
  {
    status = update_state_variables(ptr_SP_p,Stress_k,Plastic_Flow_k,Lambda_k,delta_lambda_k);
    if(status)
    {
      fprintf(stderr,"%s %s \n%s %s\n",
      "Error in the function",__func__,
      "File",__FILE__);
      return EXIT_FAILURE;
    }
  }

  
  return EXIT_SUCCESS;
}

/**************************************************************/

static Model_Parameters fill_model_paramters(Material MatProp)
{
  Model_Parameters Params;

  Params.friction_angle = MatProp.phi_Frictional;
  Params.dilatancy_angle = MatProp.psi_Frictional;

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

static int eval_K_MN(
  double * K_MN,
  const double phi)
{
    double DSQ_sin_phi = sin(phi)*sin(phi);
    
    (*K_MN) = (9.0 - DSQ_sin_phi)/(1.0 - DSQ_sin_phi);
 
  return EXIT_SUCCESS;
}
  
/**************************************************************/

static int eval_Yield_Function(
  double * Yield_Function,
  const double I1,
  const double I2,
  const double I3,
  const double friction_angle)
{
  int status = 0;
  double K_MN = 0;

  status = eval_K_MN(&K_MN,friction_angle);
  if(status)
  {
    fprintf(stderr,"%s %s \n%s %s\n",
    "Error in the function",__func__,
    "File",__FILE__);
    return EXIT_FAILURE;
  }

  (*Yield_Function) = K_MN*I3/I2 - I1;

  return EXIT_SUCCESS;
}

/**************************************************************/

static int eval_d_Yield_Function_d_stress(
  double * d_F_d_stress,
  double * Stress,
  const double I1,
  const double I2,
  const double I3,
  const double friction_angle)
  {
    int status = 0;
    double K_MN = 0.0;

    status = eval_K_MN(&K_MN,friction_angle);
    if(status)
    {
        fprintf(stderr,"%s %s \n%s %s\n",
        "Error in the function",__func__,
        "File",__FILE__);
        return EXIT_FAILURE;
    }

    d_F_d_stress[0] = K_MN*DSQR(Stress[1]*Stress[2]/I2) - 1;
    d_F_d_stress[1] = K_MN*DSQR(Stress[0]*Stress[2]/I2) - 1;
    d_F_d_stress[2] = K_MN*DSQR(Stress[0]*Stress[1]/I2) - 1;

    return EXIT_SUCCESS;
  } 


/**************************************************************/

static int eval_Plastic_Potential(
  double * Plastic_Potential,
  const double I1,
  const double I2,
  const double I3,
  const double dilatancy_angle)
  {
    int status = 0;
    double K_MN = 0.0; 

    status = eval_K_MN(&K_MN,dilatancy_angle);
    if(status)
    {
      fprintf(stderr,"%s %s \n%s %s\n",
      "Error in the function",__func__,
      "File",__FILE__);
      return EXIT_FAILURE;
    }

    (*Plastic_Potential) = K_MN*I3/I2 - I1;
    
    return EXIT_SUCCESS;
  } 

/**************************************************************/

static int eval_d_Plastic_Potential_d_stress(
  double * d_G_d_stress,
  double * Stress,
  const double I1,
  const double I2,
  const double I3,
  const double dilatancy_angle)
  {
    int status = 0;
    double K_MN = 0.0;

    status = eval_K_MN(&K_MN,dilatancy_angle);
    if(status)
    {
        fprintf(stderr,"%s %s \n%s %s\n",
        "Error in the function",__func__,
        "File",__FILE__);
        return EXIT_FAILURE;
    }

    d_G_d_stress[0] = K_MN*DSQR(Stress[1]*Stress[2]/I2) - 1;
    d_G_d_stress[1] = K_MN*DSQR(Stress[0]*Stress[2]/I2) - 1;
    d_G_d_stress[2] = K_MN*DSQR(Stress[0]*Stress[1]/I2) - 1;

    return EXIT_SUCCESS;
  } 

/**************************************************************/

static int eval_dd_Plastic_Potential_dd_stress(
  double * Hess_G,
  double * Stress,
  const double I1,
  const double I2,
  const double I3,
  const double dilatancy_angle)
  {

    int status = 0;
    double K_MN = 0.0;

    status = eval_K_MN(&K_MN,dilatancy_angle);
    if(status)
    {
        fprintf(stderr,"%s %s \n%s %s\n",
        "Error in the function",__func__,
        "File",__FILE__);
        return EXIT_FAILURE;
    }

    Hess_G[0] = - 2*K_MN*DSQR(Stress[1]*Stress[2]/I2)*((Stress[1] + Stress[2])/I2);
    Hess_G[4] = - 2*K_MN*DSQR(Stress[0]*Stress[2]/I2)*((Stress[0] + Stress[2])/I2);
    Hess_G[8] = - 2*K_MN*DSQR(Stress[0]*Stress[1]/I2)*((Stress[0] + Stress[1])/I2);

    Hess_G[1] = Hess_G[3] = 2*K_MN*(Stress[0]*Stress[1]*cbrt(Stress[2]))/cbrt(I2);
    Hess_G[2] = Hess_G[6] = 2*K_MN*(Stress[0]*cbrt(Stress[1])*Stress[2])/cbrt(I2);
    Hess_G[5] = Hess_G[7] = 2*K_MN*(cbrt(Stress[0])*Stress[1]*Stress[2])/cbrt(I2);

    return EXIT_SUCCESS;
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
  double delta_lambda,
  double Lambda_k,
  Model_Parameters Params)
  {
    int status = 0;
    double friction_angle = Params.friction_angle;
    double dilatancy_angle = Params.dilatancy_angle;
    double I1 = Stress_k[0] + Stress_k[1] + Stress_k[2];
    double I2 = Stress_k[0]*Stress_k[1] + Stress_k[1]*Stress_k[2] + Stress_k[0]*Stress_k[2];
    double I3 = Stress_k[0]*Stress_k[1]*Stress_k[2];
    double Strain_e_k[3];
    double F = 0;

    if((fabs(I1)/3 < 0.0))
    {
      fprintf(stderr,"%s: %s %i %s \n",
      "Error in Frictional_Monolithic()",
      "The stress state of particle",Particle_Idx,"has reach the apex");
      return EXIT_FAILURE;
    }

    eval_strain(Strain_e_k,Stress_k,Params.CC);
    

    status = eval_d_Plastic_Potential_d_stress(Grad_G,Stress_k,I1,I2,I3,dilatancy_angle);
    if(status)
    {
      fprintf(stderr,"%s %s \n%s %s\n",
      "Error in the function",__func__,
      "File",__FILE__);
      return EXIT_FAILURE;
    }

    status = eval_Yield_Function(&F,I1,I2,I3,friction_angle);
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
    Residual[3] = F;

    /*
      Compute absolute error from the residual
    */
   (*Norm_Residual) = 0.0;
   
    for(int A = 0 ; A<4 ; A++)
    {
      (*Norm_Residual) += DSQR(Residual[A]);
    }

    (*Norm_Residual) = pow((*Norm_Residual),0.5);

    return EXIT_SUCCESS;
  }

/**************************************************************/

static bool check_convergence(
  double Error,
  int Iter)
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

    if(Error0 < TOL_Error_0)
    {
      return true;
    }
  }
  else
  {
    Error_relative = Error/Error0;
  }


  if((Error < TOL_Error_0) 
  || (Error_relative < TOL_Error_i) 
  || (Iter > MAX_Iter))
  {
    if(Iter >= MAX_Iter) 
    {
      fprintf(stderr,"%s %s: %s %i \n",
      "Error in",__func__,
      "Maximm number of iteration reached in particle",Particle_Idx);
      fprintf(stderr,"Iter: %i | Error0: %e | Error: %e | Erro_rel: %e \n",Iter,Error0,Error,Error_relative);
    }

    return true;
  }
  else
  {
    return false;
  }

}

/**************************************************************/

static int assemble_tangent_matrix(
  double * Tangent_Matrix,
  double * d_G_d_stress,
  double * Stress_k, 
  double delta_lambda,
  double Lambda_k, 
  double Norm_Residual,
  Model_Parameters Params)
  {
    int status = 0;
    double rcond = 0;
    double friction_angle = Params.friction_angle;
    double dilatancy_angle = Params.dilatancy_angle;
    double I1 = Stress_k[0] + Stress_k[1] + Stress_k[2];
    double I2 = Stress_k[0]*Stress_k[1] + Stress_k[1]*Stress_k[2] + Stress_k[0]*Stress_k[2];
    double I3 = Stress_k[0]*Stress_k[1]*Stress_k[2];

    double dd_G_dd_stress[9]; 
    double d_F_d_stress[3] = {0.0, 0.0, 0.0};

    status = eval_dd_Plastic_Potential_dd_stress(dd_G_dd_stress,Stress_k,I1,I2,I3,dilatancy_angle); 
    if(status)
    {
      fprintf(stderr,"%s %s \n%s %s\n",
      "Error in the function",__func__,
      "File",__FILE__);
      return EXIT_FAILURE;
    }

    /* Evaluate yield function derivatives */
    status = eval_d_Yield_Function_d_stress(d_F_d_stress,Stress_k,I1,I2,I3,friction_angle);
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
    Tangent_Matrix[3] = d_G_d_stress[0];

    /* Second row */
    Tangent_Matrix[4] = Params.CC[3] + delta_lambda*dd_G_dd_stress[3];
    Tangent_Matrix[5] = Params.CC[4] + delta_lambda*dd_G_dd_stress[4];
    Tangent_Matrix[6] = Params.CC[5] + delta_lambda*dd_G_dd_stress[5];
    Tangent_Matrix[7] = d_G_d_stress[1];

    /* Third row */    
    Tangent_Matrix[8] = Params.CC[6] + delta_lambda*dd_G_dd_stress[6];
    Tangent_Matrix[9] = Params.CC[7] + delta_lambda*dd_G_dd_stress[7];
    Tangent_Matrix[10] = Params.CC[8] + delta_lambda*dd_G_dd_stress[8];
    Tangent_Matrix[11] = d_G_d_stress[2];

    /* Fourth row */    
    Tangent_Matrix[12] = d_F_d_stress[0];
    Tangent_Matrix[13] = d_F_d_stress[1];
    Tangent_Matrix[14] = d_F_d_stress[2];
    Tangent_Matrix[15] = 0.0 + Norm_Residual;
    
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
   //   return EXIT_FAILURE;
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
  int N_rows = 4;
  int N_cols = 4;
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
  double * Residual)
{
  int status = 0;
  int Order = 4;
  int LDA   = 4;
  int LDB   = 4;
  char  TRANS = 'N';
  int   INFO = 3;
  int * IPIV = (int *)malloc(4*sizeof(int));
  int NRHS = 1;

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
  dgetrs_(&TRANS,&Order,&NRHS,Tangent_Matrix,&LDA,IPIV,Residual,&LDB,&INFO);
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

static int update_variables(
  double * Residual,
  double * Stress_k,
  double * delta_lambda,
  double * lambda_k,
  double lambda_n)
{

  Stress_k[0] -= Residual[0];
  Stress_k[1] -= Residual[1];
  Stress_k[2] -= Residual[2];
  (*delta_lambda) -= Residual[3];
  (*lambda_k) = (*delta_lambda) + lambda_n; 

  return EXIT_SUCCESS;
}


/**************************************************************/

static int update_state_variables(
  State_Parameters * ptr_SP_p,
  double * Stress_k,
  double * Plastic_Flow, 
  double Lambda_k,
  double delta_lambda)
{

  if(delta_lambda < 0)
  {
    fprintf(stderr,"%s %s: %s \n%s %s\n",
    "Error in the function",__func__,
    "The increment of the discrete plastic multiplier (delta_lambda) is less than 0",
    "File",__FILE__);
    return EXIT_FAILURE;
  }

  ptr_SP_p->Equiv_Plast_Str = Lambda_k;
  ptr_SP_p->Stress = Stress_k;
  ptr_SP_p->Increment_E_plastic[0] = delta_lambda*Plastic_Flow[0];
  ptr_SP_p->Increment_E_plastic[1] = delta_lambda*Plastic_Flow[1];
  ptr_SP_p->Increment_E_plastic[2] = delta_lambda*Plastic_Flow[2];


  return EXIT_SUCCESS;
}

/**************************************************************/