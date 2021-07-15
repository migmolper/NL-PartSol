#include "nl-partsol.h"

/* List of solver implemented :
   - Newton_Rapson
   - The Conjugate Gradient method without preconditioning (CG)
   - The Conjugate Gradient method with the Jacobi preconditioner (CGJ)
   - The Generalized Minimal Resituals (GMRes) -> Not yet implemented...
*/


/*********************************************************************/


Matrix Newton_Rapson(Matrix(* Function)(Matrix, Matrix),Matrix Parameter_F,
		     Matrix(* Jacobian)(Matrix, Matrix),Matrix Parameter_J,
		     Matrix Y,Matrix X)

/*!
 *  Newton-Rapson method to solve non-linear sistems of equations :
 * Y = Y(X) -> We solve -> F(X) =  Y - Y(X) = 0
 * F(X + DeltaX) = F(X) + J(X)*DeltaX = 0 ---> DeltaX = - J(X)^{-1}*F(X)
 * Inputs :
 * - Y : Value of the function
 * - Function(X,Parameter_F) : Pointer to function to solve
 * - Parameter_F : F function optional parameters
 * - Jacobian(X,Parameter_J) : Pointer to the jacobian of the function 
 * - Parameter_J : Jacobian optional parameters
 * - X : Initial value of the objetive 
 *
 */
{

  /* Auxiliar variables */
  Matrix F_x;
  Matrix Y_x;
  Matrix dY_dX;
  Matrix dY_dX_m1;
  Matrix DeltaX;
  double TOL_NormDeltaX = pow(10,-23);
  double NormDeltaX = pow(10,4);
  int Num_Iter = 20;  
  int Iter_i = 0;
  int Bool;

  /* Allocate matrix for the function */
  F_x = alloc__MatrixLib__(Y.N_rows,Y.N_cols);
  
  /* 0º Check the convergence criterium */
  while( (NormDeltaX > TOL_NormDeltaX )
	 && (Iter_i < Num_Iter) ){
    
    /* 1º Get F(x) = Y - Y(x) */
    Y_x = Function(X,Parameter_F);
    for(int i = 0 ; i<Y.N_rows*Y.N_cols ; i++){
      F_x.nV[i] = Y.nV[i] - Y_x.nV[i];
    }
    free__MatrixLib__(Y_x);
    
    /* 2º Get the jacobian matrix in X0 DY_dX */
    /* Implement the numerical solution of the Jacobian for cases where the 
       Jacobian is not easy to derive */
    dY_dX = Jacobian(X,Parameter_J);

    /* 3º Solve the sistem DY_dX(X0)*DeltaX = F(X0) -> DeltaX */
    Bool = dY_dX.N_cols>3;
    switch(Bool)
    {
    case 0 : /* If the size of the Jacobian is less than 4, use analitical */
      dY_dX_m1 = inverse__MatrixLib__(dY_dX);
      free__MatrixLib__(dY_dX);
      DeltaX = matrix_product__MatrixLib__(dY_dX_m1,F_x);      
      free__MatrixLib__(dY_dX_m1);
      break;
    case 1 : /* If the size of the Jacobian is great than 4, use numerical */
      DeltaX = Jacobi_Conjugate_Gradient_Method(dY_dX,F_x,DeltaX);
      free__MatrixLib__(dY_dX);
      break;
    default :
      exit(EXIT_FAILURE);
    }

    /* 4º Update the variables of the convergence criterium */
    NormDeltaX = norm__MatrixLib__(DeltaX,2);
    Iter_i++;

    /* 4º Update the solution and free memory */
    X = increment__MatrixLib__(X,DeltaX);
    free__MatrixLib__(DeltaX);
        
  }

  /* Free array with the function value */
  free__MatrixLib__(F_x);
  
  /* 6º Return X */
  return X;
  
}

/*********************************************************************/

Matrix Conjugate_Gradient_Method(Matrix K, Matrix F, Matrix U0)
/*!
 *
 *  Practical aspects of the finite element method.
 * M.Pastor, P.Mira, J.A.Fernández Merodo
 *
 *  DOI : 10.1080/12795119.2002.9692737
 *  
 *  Section 6.1 : Conjugate Gradient Method with Preconditioning
 *  One of the most effective and simple iterative methods 
 * (when used with preconditioning) for solving $Ax = b$ is
 * the conjugate gradient algorithm. The algorithm is based on the idea that
 * the solution $Ax = b$ minimizes the total potential 
 * $\Pi = \frac{1}{2} x^T A x - x^T b $. Hence, the task in the iteration is,
 * given an approximate $x^k$ to $x$ for wich the potential is $\Pi^k$,
 * to find an improved aproximation $x^{k+1}$ for wich $\Pi^{k+1}<\Pi^k$.
 * Howerver, not only do we want the total potential to decrease each iteration
 * but we also want the total potential to decrease in each iteration but we also 
 * want $x^{k+1}$ to be calculate efficiently and the decrease in the total
 * potential to occur rapidly. Then the iteration will converge fast.
 *
 *  In the conjugate gradient method, we use in the kth iteration the linearly
 * independent vectors $p^1,p^2,p^3,\ldots,p^k$ and calculate the minimum of 
 * the potential in the space of the potential in the space spanned by these vectors.
 * This gives $x^{k+1}$. Also, we stablish the additional basis vector $p^{k+1}$ 
 * used in the subsequent iteration.
 *
 *  The algorithm can be summarized as follows:
 * \begin{enumerate}
 * 
 * \item Choose the iteration vector $x^1$ (frequently $x^1$ is the null vector).
 * 
 * \item Calculate the residual $r^1 = b - Ax^1$. If $r^1 = 0$, quit.
 *
 *  \item Else : 
 * \begin{enumerate}
 * \item $p^1 = r^1$
 * \item Calculate for $k = 1,2, \ldots$
 * \begin{equation}
 * \alpha^k = \frac{r^{k^T} r^k}{p^{k^T} A p^k}
 * x^{k+1} = x^k + \alpha^k p^k
 * r^{k+1} = r^k - \alpha^k A p^k
 * \beta^k = \frac{r^{k+1^T} r^{k+1}}{r^{k^T} A r^k}
 * p^{k+1} = p^{k+1} + \beta^k p^k
 * \end{equation}
 * \end{enumerate}
 * 
 * \end{enumerate}
 * 
 * We continue iterating until $||r^k|| \leseq \epsilon$,
 * where $\epsilon$ in the convergence tolerance. 
 * A convergence criterion on $||x^k||$ could also be used.
 * 
 * The conjugate gradient algorithm satisfies two important
 * orthogonality properties regarding the direction vectors $p_i$ 
 * and the residual $r_i$
 *
 * See also : 
 * https://en.wikipedia.org/wiki/Conjugate_gradient_method#The_preconditioned_conjugate_gradient_method
 *
 */
{

   /* First we check if the input data */
  if((K.N_cols != K.N_rows) ||
     ( (K.N_rows != F.N_rows) || (F.N_cols != 1)) ||
     ( (K.N_rows != U0.N_rows) || (U0.N_cols != 1) ))
    {
      printf("%s : %s \n",
	     "Error in Conjugate_Gradient_Method()",
	     "Wrong input data !");
      exit(EXIT_FAILURE);
    }

  int N = K.N_rows;
  double aux;
  double Tol_r;
  double alpha_k,beta_k;
  double dividend,divisor;
  double Norm_r;
  Matrix r_k,r_k1;
  Matrix p;
  Matrix U;
  int STOP_CRITERIA = 0;
  int Num_Iter,Num_Iter_Max;

  Num_Iter = 0;
  Num_Iter_Max = 25;
  Tol_r = 0.0001;

  /* Allocate the residual and the p basis arrays */
  r_k = alloc__MatrixLib__(N,1);
  r_k1 = alloc__MatrixLib__(N,1);
  p = alloc__MatrixLib__(N,1);
  
  /* Set initial solution to the initial solution */
  U = U0;

  /* Do initial step */
  for(int i = 0 ; i<N ; i++){

    /* Calculate the residual $r^1 = F - A x^1$ */
    aux = 0;
    for(int j = 0 ; j<N ; j++){
      aux += K.nM[i][j]*U.nV[j];
    }
    r_k.nV[i] = F.nV[i] - aux;

    /* Set p^1 = r^1   : */
    p.nV[i] = r_k.nV[i];    
  }

  /* Calcule the stopping criteria in the first step */
  Norm_r = norm__MatrixLib__(r_k,2);
  if( (Norm_r < Tol_r) ){
    STOP_CRITERIA = 1;
    puts("The initial solution was correct");   
  }

    
  /* iterate */
  while(STOP_CRITERIA == 0){

    /* 1th step : Get alpha */
    alpha_k = 0;
    dividend = 0;
    divisor = 0;
    for(int i = 0 ; i<N ; i++){
      /* Scalar product -> dividend = r_k^T \cdot r_k */
      dividend +=  r_k.nV[i]*r_k.nV[i];
      /* Scalar product -> divisor = r_k^T \cdot K \cdot r_k */
      aux = 0;
      for(int j = 0 ; j<N ; j++){
	aux += K.nM[i][j]*p.nV[j];
      }
      divisor += p.nV[i]*aux;      
    }
    alpha_k = dividend/divisor;
    
    for(int i = 0 ; i<N ; i++){
      /* 2th Step : Get the solution array (Update) */
      U.nV[i] += alpha_k*p.nV[i];
      /* 3th Step : Calcule the residual */
      aux = 0;
      for(int j = 0 ; j<N ; j++){
	aux += K.nM[i][j]*p.nV[j];
      }      
      r_k1.nV[i] = r_k.nV[i] - alpha_k*aux;
    }

    /* Calcule the stopping criteria */
    Norm_r = norm__MatrixLib__(r_k1,2);
    Num_Iter++;
    if( (Norm_r < Tol_r) ||
	(Num_Iter >= Num_Iter_Max) ){
      STOP_CRITERIA = 1;

      if(Num_Iter < Num_Iter_Max){
	printf("%s % i %s \n",
	       "Convergence criteria reached after",
	       Num_Iter,"iterations");
      }
      else{
	printf("%s : %s \n \t %s : %f \n",
	       "Warning in Conjugate_Gradient_Method",
	       "Maximum number of iterations",
	       "Norm of r",Norm_r);
      }
      
    }
    
    /* 4th step : */
    beta_k = 0;
    dividend = 0;
    divisor = 0;
    for(int i = 0 ; i<N ; i++){
      /* Scalar product -> dividend = r_{k+1}^T \cdot r_{k+1} */
      dividend +=  r_k1.nV[i]*r_k1.nV[i];
      /* Scalar product -> divisor = r_k^T \cdot r_k */
      divisor += r_k.nV[i]*r_k.nV[i];      
    }
    beta_k = dividend/divisor;

    for(int i = 0 ; i<N ; i++){
      /* 5th step : Update the basis vector p */
      p.nV[i] = r_k1.nV[i] + beta_k*p.nV[i];

      /* 6th step : Update r_{k} */
      r_k.nV[i] = r_k1.nV[i];
    }

  }

  /* Free memory */
  free__MatrixLib__(r_k);
  free__MatrixLib__(r_k1);
  free__MatrixLib__(p);
  
  return U;
  
}

/*********************************************************************/

Matrix Jacobi_Conjugate_Gradient_Method(Matrix K, Matrix F, Matrix U0)
/*!
 *  To increase the rate of convergence of Conjugate_Gradient_Method(),
 * preconditioning is used. The basic idea is that instead of solving $K U = F$,
 * we solve :
 * \begin{equation}
 * \tilde{K}^{-1}K U =\tilde{K}^{-1}F
 * \end{equation}
 * where $\tilde{K}$ is called the preconditioner. 
 * The objective with this transformation is to obtain a matrix $\tilde{K}^{-1}K$
 * with a much improved conditioned number choosing an easy inverting 
 * matrix $tilde{A}$. Various preconditioners have been proposed, 
 * the the choose of the diagonal part of $K$ results in the
 * Jacoby Conjugate method (JCG).
 * The new algorithm introduces an additional set of vectors $z^k$ defined by:
 * \begin{equation}
 * z^k = \tilde{K}^{-1} r^k 
 * \end{equation}
 * who modifies the definition of $\alpha^k$, $\beta^k$, $p^k$ :
 * \begin{equation}
 * \alpha^k = \frac{z^k^T r^k}{p^k^T A p^k}
 * \beta^k = \frac{z^{k+1^T} r^{k+1}}{z^{k^T} r^k}
 * p^{k+1} = z^{k+1} + \beta^k p^k
 * \end{equation}
 *
 */
{

   /* First we check if the input data */
  if((K.N_cols != K.N_rows) ||
     ( (K.N_rows != F.N_rows) || (F.N_cols != 1)) ||
     ( (K.N_rows != U0.N_rows) || (U0.N_cols != 1) ))
    {
      printf("%s : %s \n",
	     "Error in Jacobi_Conjugate_Gradient_Method()",
	     "Wrong input data !");
      exit(EXIT_FAILURE);
    }
  
  int N = K.N_rows;
  double aux;
  double Tol_r;
  double alpha_k,beta_k;
  double dividend,divisor;
  Matrix K_l; /* Lumped matrix */
  double Norm_r;
  Matrix r_k,r_k1;
  Matrix z_k,z_k1;
  Matrix p;
  Matrix U;
  int STOP_CRITERIA = 0;
  int Num_Iter,Num_Iter_Max;

  Num_Iter = 0;
  Num_Iter_Max = 25;
  Tol_r = 0.0000001;

  /* Allocate the residual arrays and the basis p array */
  r_k = alloc__MatrixLib__(N,1);
  r_k1 = alloc__MatrixLib__(N,1);
  z_k = alloc__MatrixLib__(N,1);
  z_k1 = alloc__MatrixLib__(N,1);
  p = alloc__MatrixLib__(N,1);
  
  
  /* Set initial solution to the initial solution */
  U = U0;

  /* Do initial step */
  for(int i = 0 ; i<N ; i++){
  
    /* Calculate the residual $r^1 = F - A x^1$ */
    aux = 0;
    for(int j = 0 ; j<N ; j++){
      aux += K.nM[i][j]*U.nV[j];
    }
    r_k.nV[i] = F.nV[i] - aux;

    /* Calcule the z vector $z^k = \tilde{K}^{-1} r^k $ */
    z_k.nV[i] = (1/K.nM[i][i])*r_k.nV[i];

    /* Set p^1 = z^1   : */
    p.nV[i] = z_k.nV[i];    
  }

  /* Calcule the stopping criteria in the first step */
  Norm_r = norm__MatrixLib__(r_k,2);
  if( (Norm_r < Tol_r) ){
    STOP_CRITERIA = 1;
    puts("The initial solution was correct");   
  }

  /* Get the Lumped-Mass matrix */
  K_l = lumped__MatrixLib__(K);

  /* iterate */
  while(STOP_CRITERIA == 0){

    /* 1th step : Get alpha */
    alpha_k = 0;
    dividend = 0;
    divisor = 0;
    for(int i = 0 ; i<N ; i++){
      /* Scalar product -> dividend = r_k^T \cdot r_k */
      dividend +=  r_k.nV[i]*z_k.nV[i];
      /* Scalar product -> divisor = r_k^T \cdot K \cdot r_k */
      aux = 0;
      for(int j = 0 ; j<N ; j++){
	aux += K.nM[i][j]*p.nV[j];
      }
      divisor += p.nV[i]*aux;      
    }
    alpha_k = dividend/divisor;
    
    for(int i = 0 ; i<N ; i++){
      /* 2th Step : Get the solution array (Update) */
      U.nV[i] += alpha_k*p.nV[i];
      /* 3th Step : Calcule the residual */
      aux = 0;
      for(int j = 0 ; j<N ; j++){
	aux += K.nM[i][j]*p.nV[j];
      }      
      r_k1.nV[i] = r_k.nV[i] - alpha_k*aux;
    }

    /* Calcule the stopping criteria */
    Norm_r = norm__MatrixLib__(r_k1,2);
    Num_Iter++;
    if( (Norm_r < Tol_r) ||
	(Num_Iter >= Num_Iter_Max) ){
      STOP_CRITERIA = 1;

      if(Num_Iter < Num_Iter_Max){
	printf("%s %i %s \n",
	       "Convergence criteria reached after",
	       Num_Iter,
	       "iterations");
      }
      else{
	printf("%s : %s \n \t %s : %f \n",
	       "Warning in Jacobi_Conjugate_Gradient_Method",
	       "Maximum number of iterations",
	       "Norm of r",Norm_r);
      }      
    }
       
    beta_k = 0;
    dividend = 0;
    divisor = 0;
    for(int i = 0 ; i<N ; i++){
      /* 4th step : */
      z_k1.nV[i] = ((double)1/K_l.nV[i])*r_k1.nV[i];
      /* 5th step :*/  
      /* Scalar product -> dividend = r_{k+1}^T \cdot r_{k+1} */
      dividend +=  z_k1.nV[i]*r_k1.nV[i];
      /* Scalar product -> divisor = r_k^T \cdot r_k */
      divisor += z_k.nV[i]*r_k.nV[i];      
    }
    beta_k = dividend/divisor;

    for(int i = 0 ; i<N ; i++){
      /* 6th step : Update the basis vector $p$ */
      p.nV[i] = z_k1.nV[i] + beta_k*p.nV[i];

      /* 7th step : Update the residual $r_{k}$ and the $z_k$ vector  */
      r_k.nV[i] = r_k1.nV[i];
      z_k.nV[i] = z_k1.nV[i];
    }

  }

  /* Free memory */
  free__MatrixLib__(r_k);
  free__MatrixLib__(r_k1);
  free__MatrixLib__(z_k);
  free__MatrixLib__(z_k1);
  free__MatrixLib__(p);
  free__MatrixLib__(K_l);
   
  return U;
  
}

/*********************************************************************/

Matrix One_Iteration_Lumped(Matrix K_l, Matrix F, Matrix U0){

  /* 0º First we check if the input data */
  if( K_l.N_cols*K_l.N_rows != F.N_cols*F.N_rows ){
    printf("%s : %s \n",
	   "Error in One_Iteration_Lumped()",
	   "Wrong input data !");
    exit(EXIT_FAILURE);
  }

  /* 1º Initialice the solution array */
  Matrix U = U0;

  /* 2º Solve the sistem */
  for(int i = 0 ; i<F.N_rows ; i++){
    U.nV[i] = (double)1/K_l.nV[i]*F.nV[i];
  }

  /* 3º Return the solution */
  return U;

}

/*********************************************************************/
