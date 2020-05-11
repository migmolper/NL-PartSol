/*! 
  \file InOutFun.h
  \brief File with the prototype in/out functions
*/

#ifndef _INOUTFUN_H_
#define _INOUTFUN_H_

/*! 
  \fn int parse(char ** Words, char * String, char * Delimiters)
  
  \brief Parse string in to its components 
  
  \param String : Input string
  \param Delimiters : parser symbol o list of symbols
  \param Words : Output string with the words parsed

*/
int parse(char **, char *, char *);

/*****************************************************************/

/*
  \def void generate_route(char * Route, char * File)

  \brief Generate the relative route to a file

  \param : Route
  \param : File symbol o list of symbols

*/
void generate_route(char *, char *);

/*****************************************************************/

Matrix Read_CSV(char *, int);
Curve ReadCurve(char *);
void free_Curve(Curve);
Mesh ReadGidMesh(char *);
ChainPtr File2Chain(char *);

/*! 
  \fn void print_Status(char Message,int Time)
  
  \brief Print screen message in a safety way.
  
  \param Message : Message to print
  \param Time : Current step of the simulation 
*/
void print_Status(char *,int);

/*****************************************************************/

/*! 
  \fn void print_step(int Time,double DeltaTimeStep)
 
  \brief Print current step in a safety way.
 
  \param Time : Current step of the simulation 
  \param DeltaTimeStep : Current time step
 */
void print_step(int,double);

/*****************************************************************/

/* Read .gfd format */
Mesh GramsBox(char *);
Boundaries GramsBoundary(char *,int);
GaussPoint GramsSolid2D(char *,Mesh);
void GramsTime(char * );
Material * GramsMaterials(char *, GaussPoint, int);
void GramsInitials(char *, GaussPoint, int);
void GramsShapeFun(char * );
void GramsOutputs(char * );
Load * GramsNeumannBC(char *, int, int);
Load * GramsBodyForces(char *, int, int);

/* Print .vtk format */
void WriteVtk_MPM(char *, GaussPoint, Matrix, int);
void WriteVtk_FEM(char *, Mesh, Matrix, int);
void WriteVtk_Float_Scalar(char *, Matrix);
void WriteVtk_Float_Vector(char *, Matrix);
void WriteVtk_Float_Tensor(char *, Matrix);

/* Print Gnuplot */
void WriteGnuplot(Matrix, Matrix, double, double, int, int, char[20]);

#endif
