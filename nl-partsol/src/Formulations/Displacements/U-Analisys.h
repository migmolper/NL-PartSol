#ifndef _U_ANALISYS_H_
#define _U_ANALISYS_H_

// clang-format off
#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>
#include "Macros.h"
#include "Types.h"
#include "Globals.h"
#include "Matlib.h"
// clang-format on

/*!
 * \fn Fields allocate_U_vars__Fields__(int NumParticles)
 * 
 * \brief This function is devoted to allocate memory to store the field data
 * 
 * \param NumParticles : Number of particles in the domain
 * */
Fields allocate_U_vars__Fields__(int);

/*!
 * \fn void free_U_vars__Fields__(Fields Phi)
 * 
 * \brief This function is devoted to free memory
 * 
 * \param Fields : Variable that stores all the field information
 * */
void free_U_vars__Fields__(Fields);
/*****************************************************************/
#endif