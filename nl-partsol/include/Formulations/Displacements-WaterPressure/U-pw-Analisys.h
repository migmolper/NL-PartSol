#ifndef _UPW_ANALISYS_H_
#define _UPW_ANALISYS_H_

#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>
#include "Macros.h"
#include "Types.h"
#include "Matlib.h"


/*!
 * \fn Fields allocate_upw_vars__Fields__(int NumParticles);
 * 
 * \brief This function is devoted to allocate memory to store the field data for a upw analisys
 *  
 * \param NumParticles : Number of particles in the domain 
 * */
Fields allocate_upw_vars__Fields__(int);

/*!
 * \fn void free_upw_vars__Fields__(Fields Phi)
 *  
 * \brief This function is devoted to free memory
 * 
 * \param Fields : Variable that stores all the field information
 * */
void free_upw_vars__Fields__(Fields);
/*****************************************************************/


#endif