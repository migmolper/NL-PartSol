
#ifdef USE_PETSC
#ifndef _KSP_PETSC_H_
#define _KSP_PETSC_H_

#include <petscksp.h>

/*!
    \brief KSP solves a system of linear equations
    A * X = B  or  A**T * X = B
    it is an abstract PETSc object that manages all Krylov methods. 
    This is the object that manages the linear solves in PETSc
    (even those such as direct solvers that do no use Krylov accelerators). 

    \param[in] ptr_Tangent_Stiffness Pointer to the tangent matrix of the problem
    \param[in,out] ptr_Residual Pointer to the residual, returs the solution
    \param[in] Nactivedofs Number of active degrees of freedom

*/
int krylov_PETSC(
    Mat * ptr_Tangent_Stiffness,
    Vec * ptr_Residual,
    unsigned Nactivedofs);

#endif
#endif