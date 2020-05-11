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
  \fn void generate_route(char * Route, char * File)

  \brief Generate the relative route to a file

  \param Route : Route to the file
  \param File : Route to the goberning file

*/
void generate_route(char *, char *);

/*****************************************************************/

/*! 
  \fn Matrix Read_CSV(char * File, int NumData)

  \brief Read text separated by commas

  \param File : Route to the file
  \param NumData : Number of elements in the file

*/
Matrix Read_CSV(char *, int);

/*****************************************************************/

/*! 
  \fn Curve ReadCurve(char * File)

  \brief Generate a user defined curve for loads, boundary conditions ...

  List of functions and interface:
  CONSTANT_CURVE SCALE#double
  RAMP CURVE SCALE#double 
  HEAVISIDE_CURVE SCALE#double Tc#integer
  DELTA_CURVE SCALE#double Tc#integer
  HAT_CURVE SCALE#double T0#integer T1#integer
  CUSTOM_CURVE

  Example :
  DAT_CURVE NUM#3
  HEAVISIDE_CURVE SCALE#2.5 Tc#integer
  CONSTANT_CURVE SCALE#0


  \param File : Route to the file with the instructions

*/
Curve ReadCurve(char *);

/*****************************************************************/

/*!

  \fn void fill_ConstantCurve(Curve DatCurve,char * Prop_1)

  \brief Generate a constant curve

  \param DatCurve : Curve to fill 
  \param Prop_1

*/
void fill_ConstantCurve(Curve,char *);

/*****************************************************************/

/*!

  \fn void fill_RampCurve(Curve DatCurve,char * Prop_1)

  \brief Generate the ramp curve

  \param DatCurve : Curve to fill 
  \param Prop_1

*/
void fill_RampCurve(Curve,char *);

/*****************************************************************/

/*!

  \fn void fill_HeavisideCurve(Curve DatCurve,char * Prop_1,char * Prop_2);

  \brief Generate the heaviside curve

  \param DatCurve : Curve to fill 
  \param Prop_1
  \param Prop_2

*/
void fill_HeavisideCurve(Curve,char *,char *);

/*****************************************************************/

/*!

  \fn void fill_DeltaCurve(Curve DatCurve,char * Prop_1,char * Prop_2);

  \brief Generate the delta curve

  \param DatCurve : Curve to fill 
  \param Prop_1
  \param Prop_2

*/
void fill_DeltaCurve(Curve,char *,char *);

/*****************************************************************/

/*!
  \fn void fill_HatCurve(Curve DatCurve,char * Prop_1,char * Prop_2,char * Prop_3);

  \brief Generate the hat curve

  \param DatCurve : Curve to fill 
  \param Prop_1
  \param Prop_2
  \param Prop_3

*/
void fill_HatCurve(Curve,char *,char *,char *);

/*****************************************************************/


/*! 
  \fn void free_Curve(Curve A)

  \brief Free the memory allocated by the curve function

  \param A : Curve to free
*/
void free_Curve(Curve);

/*****************************************************************/

/*
  \fn Mesh ReadGidMesh(char * MeshName)

  \brief This function is devoted to read mesh files generated by GID

  \param MeshName : Route to the mesh
*/
Mesh ReadGidMesh(char *);

/*****************************************************************/
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
