// clang-format off
#include <stdlib.h>
#include <stdbool.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "Macros.h"
#include "Types.h"
#include "Matlib.h"
// clang-format on

/* List of solver implemented :
   - Newton_Rapson
   - The Conjugate Gradient method without preconditioning (CG)
   - The Conjugate Gradient method with the Jacobi preconditioner (CGJ)
   - The Generalized Minimal Resituals (GMRes) -> Not yet implemented...
*/

/*********************************************************************/

Matrix Newton_Rapson(Matrix (*Function)(Matrix, Matrix), Matrix Parameter_F,
                     Matrix (*Jacobian)(Matrix, Matrix), Matrix Parameter_J,
                     Matrix Y, Matrix X)

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
  double TOL_NormDeltaX = pow(10, -23);
  double NormDeltaX = pow(10, 4);
  int Num_Iter = 20;
  int Iter_i = 0;
  int Bool;

  /* Allocate matrix for the function */
  F_x = alloc__MatrixLib__(Y.N_rows, Y.N_cols);

  /* 0º Check the convergence criterium */
  while ((NormDeltaX > TOL_NormDeltaX) && (Iter_i < Num_Iter)) {

    /* 1º Get F(x) = Y - Y(x) */
    Y_x = Function(X, Parameter_F);
    for (int i = 0; i < Y.N_rows * Y.N_cols; i++) {
      F_x.nV[i] = Y.nV[i] - Y_x.nV[i];
    }
    free__MatrixLib__(Y_x);

    /* 2º Get the jacobian matrix in X0 DY_dX */
    /* Implement the numerical solution of the Jacobian for cases where the
       Jacobian is not easy to derive */
    dY_dX = Jacobian(X, Parameter_J);

    /* 3º Solve the sistem DY_dX(X0)*DeltaX = F(X0) -> DeltaX */
    Bool = dY_dX.N_cols > 3;
    switch (Bool) {
    case 0: /* If the size of the Jacobian is less than 4, use analitical */
      dY_dX_m1 = inverse__MatrixLib__(dY_dX);
      free__MatrixLib__(dY_dX);
      DeltaX = matrix_product__MatrixLib__(dY_dX_m1, F_x);
      free__MatrixLib__(dY_dX_m1);
      break;
    default:
      exit(EXIT_FAILURE);
    }

    /* 4º Update the variables of the convergence criterium */
    NormDeltaX = norm__MatrixLib__(DeltaX, 2);
    Iter_i++;

    /* 4º Update the solution and free memory */
    X = increment__MatrixLib__(X, DeltaX);
    free__MatrixLib__(DeltaX);
  }

  /* Free array with the function value */
  free__MatrixLib__(F_x);

  /* 6º Return X */
  return X;
}


/*********************************************************************/
