/*! \file Fields.h
    \brief File with the prototype with the function to free memory
*/

#ifndef _FIELDS_H_
#define _FIELDS_H_


/*!
  \fn Fields allocate_Fields(int NumParticles)

  \brief This function is devoted to allocate memory to store the field data

  \param NumParticles : Number of particles in the domain 
*/
Fields allocate_Fields(int);
/*****************************************************************/

/*!
	\fn Fields allocate_upw_vars__Fields__(int NumParticles);
	
	\brief This function is devoted to allocate memory to store the field data for a upw analisys

	\param NumParticles : Number of particles in the domain 
*/
Fields allocate_upw_vars__Fields__(int NumParticles);
/*****************************************************************/

/*!
  \fn void free_Fields(Fields Phi)

  \brief This function is devoted to free memory

  \param Fields : Variable that stores all the field information
*/
void free_Fields(Fields);
/*****************************************************************/

/*!
  \fn void free_upw_vars__Fields__(Fields Phi)

  \brief This function is devoted to free memory

  \param Fields : Variable that stores all the field information
*/
void free_upw_vars__Fields__(Fields);
/*****************************************************************/

#endif
