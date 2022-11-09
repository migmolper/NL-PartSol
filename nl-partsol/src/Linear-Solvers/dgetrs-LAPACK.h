
#ifndef _DGETRS_LAPACK_H_
#define _DGETRS_LAPACK_H_

/**
 * @brief solves a system of linear equations
 *  A * X = B  or  A**T * X = B
 * with a general N-by-N matrix A using the LU factorization computed 
 * by DGETRF.
 * 
 * @param Tangent_Stiffness 
 * @param Residual 
 * @param Nactivedofs 
 * @return int 
 */
int dgetrs_LAPACK(
    double * Tangent_Stiffness,
    double * Residual,
    unsigned Nactivedofs);

#endif