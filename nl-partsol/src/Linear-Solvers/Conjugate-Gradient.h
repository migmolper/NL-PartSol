
#ifndef _CONJUGATE_GRADIENT_H_
#define _CONJUGATE_GRADIENT_H_

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "Macros.h"

/*!
 
   Practical aspects of the finite element method.
  M.Pastor, P.Mira, J.A.Fernández Merodo
 
   DOI : 10.1080/12795119.2002.9692737
 
   Section 6.1 : Conjugate Gradient Method with Preconditioning
   One of the most effective and simple iterative methods
  (when used with preconditioning) for solving $Ax = b$ is
  the conjugate gradient algorithm. The algorithm is based on the idea that
  the solution $Ax = b$ minimizes the total potential
  $\Pi = \frac{1}{2} x^T A x - x^T b $. Hence, the task in the iteration is,
  given an approximate $x^k$ to $x$ for wich the potential is $\Pi^k$,
  to find an improved aproximation $x^{k+1}$ for wich $\Pi^{k+1}<\Pi^k$.
  Howerver, not only do we want the total potential to decrease each iteration
  but we also want the total potential to decrease in each iteration but we
  also want $x^{k+1}$ to be calculate efficiently and the decrease in the total
  potential to occur rapidly. Then the iteration will converge fast.
 
   In the conjugate gradient method, we use in the kth iteration the linearly
  independent vectors $p^1,p^2,p^3,\ldots,p^k$ and calculate the minimum of
  the potential in the space of the potential in the space spanned by these
  vectors. This gives $x^{k+1}$. Also, we stablish the additional basis vector
  $p^{k+1}$ used in the subsequent iteration.
 
   The algorithm can be summarized as follows:
  \begin{enumerate}
 
  \item Choose the iteration vector $x^1$ (frequently $x^1$ is the null
  vector).
 
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
 
  We continue iterating until $||r^k|| \leseq \epsilon$,
  where $\epsilon$ in the convergence tolerance.
  A convergence criterion on $||x^k||$ could also be used.
 
  The conjugate gradient algorithm satisfies two important
  orthogonality properties regarding the direction vectors $p_i$
  and the residual $r_i$
 
  See also :
  https://en.wikipedia.org/wiki/Conjugate_gradient_method#The_preconditioned_conjugate_gradient_method
 
*/
int Conjugate_Gradient_Method(
    double * Tangent_Stiffness, 
    double * Residual, 
    double * U,
    unsigned Nactivedofs);
/*********************************************************************/

#endif