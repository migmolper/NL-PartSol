#include <stdio.h>
#include <stdlib.h>
#include "ToolsLib/TypeDefinitions.h"
#include "ToolsLib/Utils.h"

/*********************************************************************/

Matrix Conjugate_Gradient_Method(Matrix K, Matrix F, Matrix U0)
/*

  Practical aspects of the finite element method.
  M.Pastor, P.Mira, J.A.Fernández Merodo

  DOI : 10.1080/12795119.2002.9692737
  
  Section 6.1 : Conjugate Gradient Method with Preconditioning
  One of the most effective and simple iterative methods (when used with preconditioning) for solving $Ax = b$ is the conjugate gradient algorithm. The algorithm is based on the idea that the solution $Ax = b$ minimizes the total potential $\Pi = \frac{1}{2} x^T A x - x^T b $. Hence, the task in the iteration is, given an approximate $x^k$ to $x$ for wich the potential is $\Pi^k$, to find an improved aproximation $x^{k+1}$ for wich $\Pi^{k+1}<\Pi^k$. Howerver, not only do we want the total potential to decrease each iteration but we also want the total potential to decrease in each iteration but we also want $x^{k+1}$ to be calculate efficiently and the decrease in the total potential to occur rapidly. Then the iteration will converge fast.

In the conjugate gradient method, we use in the kth iteration the linearly independent vectors $p^1,p^2,p^3,\ldots,p^k$ and calculate the minimum of the potential in the space of the potential in the space spanned by these vectors. This gives $x^{k+1}$. Also, we stablish the additional basis vector $p^{k+1}$ used in the subsequent iteration.

The algorithm can be summarized as follows:
\begin{enumerate}

\item Choose the starting iteration vector $x^1$ (frequently $x^1$ is the null vector).

\item Calculate the residual $r^1 = b - Ax^1$. If $r^1 = 0$, quit.

\item Else : 
\begin{enumerate}
\item $p^1 = r^1$
\item Calculate for $k = 1,2, \ldots$
\begin{equation}
\alpha^k = \frac{r^{k^T} r^k}{p^{k^T} A p^k}
x^{k+1} = x^k + \alpha^k p^k
r^{k+1} = r^k - \alpha^k A p^k
\beta^k = \frac{r^{k+1^T} r^{k+1}}{r^{k^T} A r^k}
p^{k+1} = p^{k+1} + \beta^k p^k
\end{equation}
\end{enumerate}

\end{enumerate}

We continue iterating until $||r^k|| \leseq \epsilon$, where $\epsilon$ in the convergence tolerance. A convergence criterion on $||x^k||$ could also be used.

The conjugate gradient algorithm satisfies two important orthogonality properties regarding the direction vectors $p_i$ and the residual $r_i$

See also : 
https://en.wikipedia.org/wiki/Conjugate_gradient_method#The_preconditioned_conjugate_gradient_method

*/
{
  /* First we check if the input data */
  if((K.N_cols != K.N_rows) ||
     ( (K.N_rows != F.N_rows) || (F.N_cols != 1))){
    puts("Error in Conjugate_Gradient_Precondicionated() : Wrong input data !");
    exit(0);
  }
  
  int N = K.N_rows;
  double aux;
  double Tol_U,Tol_r;
  double alpha_k,beta_k;
  double dividend,divisor;
  Matrix Norm_U,Norm_r;
  Matrix r_k,r_kT,r_k1,r_k1T;
  Matrix p;
  Matrix U;
  int STOP_CRITERIA = 0;
  int Num_Iter,Num_Iter_Max;

  Num_Iter = 0;
  Num_Iter_Max = 25;
  Tol_U = 0.0001;
  Tol_r = 0.0001;
  
  /* Set initial solution to the initial solution */
  U = U0;

  /* Calculate the residual r^1 : 
     $r^1 = F - A x^1$
   */
  r_k = MatAlloc(N,1);
  r_k1 = MatAlloc(N,1);
  p = MatAlloc(N,1);

  for(int i = 0 ; i<N ; i++){
    aux = 0;
    for(int j = 0 ; j<N ; j++){
      aux += K.nM[i][j]*U.nV[j];
    }
    r_k.nV[i] = F.nV[i] - aux;

    /* Set p^1 = r^1   : */
    p.nV[i] = r_k.nV[i];    
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
    Norm_U = Norm_Mat(U,2);
    Norm_r = Norm_Mat(r_k1,2);
    Num_Iter++;
    if( (Norm_r.n < Tol_r) ||
	(Norm_U.n < Tol_U) ||
	(Num_Iter >= Num_Iter_Max) ){
      STOP_CRITERIA = 1;

      if(Num_Iter <= Num_Iter_Max){
	printf("Convergence criteria reached after %i iterations \n",Num_Iter);
      }
      else{
	printf("Warning not convergence reached after the maximum number of iterations: \n");
	printf("\t Norm of r: %f \n",Norm_r.n);
	printf("\t Norm of U: %f \n",Norm_U.n);
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
   
  return U;
  
}

/*********************************************************************/

Matrix Jacobi_Conjugate_Gradient_Method(Matrix K, Matrix F, Matrix U0)
/*
  To increase the rate of convergence of Conjugate_Gradient_Method(), preconditioning is used. The basic idea is that instead of solving $K U = F$, we solve :
  \begin{equation}
  \tilde{K}^{-1}K U =\tilde{K}^{-1}F
  \end{equation}
  where $\tilde{K}$ is called the preconditioner. The objective with this transformation is to obtain a matrix $\tilde{K}^{-1}K$ with a much improved conditioned number choosing an easy inverting matrix $tilde{A}$. Various preconditioners have been proposed, the the choose of the diagonal part of $K$ results in the Jacoby Conjugate method (JCG).
  The new algorithm introduces an additional set of vectors $z^k$ defined by:
  \begin{equation}
  z^k = \tilde{K}^{-1} r^k 
  \end{equation}
  who modifies the definition of $\alpha^k$, $\beta^k$, $p^k$ :
  \begin{equation}
  \alpha^k = \frac{}{}
  \beta^k = \frac{}{}
  p^{k+1} = 
  \end{equation}
 */
{

  
  
}
