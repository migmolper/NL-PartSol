/*! \file Fields.h
    \brief File with the prototype with the function to free memory
*/

#ifndef _FIELDS_H_
#define _FIELDS_H_


/*!
  \fn Fields allocate_Fields(int NumParticles)

  \brief This function is devoted to allocate memory to store the field data

  Inputs :
  \param NumParticles : Number of particles in the domain 
*/
Fields allocate_Fields(int);


/*!
  \fn void free_Fields(Fields Phi)

  \brief This function is devoted to free memory

  Inputs :
  \param Fields : Variable that stores all the field information
*/
void free_Fields(Fields);

#endif
