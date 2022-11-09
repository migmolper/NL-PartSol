
/*! \file Jacobi-Conjugate-Gradient.h
    \brief  To increase the rate of convergence of Conjugate_Gradient_Method(),
    preconditioning is used. The basic idea is that instead of solving $K U = F$,
    we solve :
    \begin{equation}
    \tilde{K}^{-1}K U =\tilde{K}^{-1}F
    \end{equation}
    where $\tilde{K}$ is called the preconditioner.
    The objective with this transformation is to obtain a matrix
    $\tilde{K}^{-1}K$ with a much improved conditioned number choosing an easy
    inverting matrix $tilde{A}$. Various preconditioners have been proposed, the
    the choose of the diagonal part of $K$ results in the Jacoby Conjugate method
    (JCG). The new algorithm introduces an additional set of vectors $z^k$
    defined by: \begin{equation} z^k = \tilde{K}^{-1} r^k \end{equation} who
    modifies the definition of $\alpha^k$, $\beta^k$, $p^k$ : \begin{equation}
    \alpha^k = \frac{z^k^T r^k}{p^k^T A p^k}
    \beta^k = \frac{z^{k+1^T} r^{k+1}}{z^{k^T} r^k}
    p^{k+1} = z^{k+1} + \beta^k p^k
    \end{equation}
*/

#ifndef _JACOBI_CONJUGATE_GRADIENT_H_
#define _JACOBI_CONJUGATE_GRADIENT_H_

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "Macros.h"


int Jacobi_Conjugate_Gradient_Method(
    double * Tangent_Stiffness, 
    double * Residual, 
    double * U,
    unsigned Nactivedofs);

#endif