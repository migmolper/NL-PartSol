
#ifndef _DGETRS_LAPACK_H_
#define _DGETRS_LAPACK_H_

#include <stdlib.h>
#include <stdio.h>
#include "Macros.h"

#ifdef __linux__
#include <lapacke.h>
#elif __APPLE__
#include <Accelerate/Accelerate.h>
#endif

/*!
 DGETRS solves a system of linear equations
    A * X = B  or  A**T * X = B
 with a general N-by-N matrix A using the LU factorization computed
 by DGETRF.
*/
int dgetrs_LAPACK(
    double * Tangent_Stiffness,
    double * Residual,
    unsigned Nactivedofs);

#endif