
#ifndef _UP_ANALISYS_H_
#define _UP_ANALISYS_H_

#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>
#include "Macros.h"
#include "Types.h"
#include "Matlib.h"


/*!
 * \fn Fields allocate_Up_vars__Fields__(int NumParticles)
 * 
 * \brief This function is devoted to allocate memory to store the field data for U-p formulations
 * 
 * \param Fields : Variable that stores all the field information
 * */
Fields allocate_Up_vars__Fields__(int);

/*!
 * \fn void free_Up_vars__Fields__(Fields Phi)
 * 
 * \brief This function is devoted to free memory for U-p formulations
 * 
 * \param Fields : Variable that stores all the field information
 * */
void free_Up_vars__Fields__(Fields);
/*****************************************************************/


#endif