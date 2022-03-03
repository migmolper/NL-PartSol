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

int Parse_string__InOutFun__(char **, char *, char *, int *);
/*****************************************************************/

FILE * Open_and_Check_File__InOutFun__(char *, int *);
/*****************************************************************/

ChainPtr File_to_Chain__InOutFun__(char *, int *);
/*****************************************************************/

Tensor Read_Vector__InOutFun__(char *, int * );
/*****************************************************************/

/*
  \fn void generate_route(char * Route, char * File)

  \brief Generate the relative route to a file

  \param Route : Route to the file
  \param File : Route to the goberning file

*/
void generate_route(char *, char *);


/*****************************************************************/


/*
  \fn int get_ResultStep(char * File_Result);

  \brief Read the step of the file

  \param File_Result : Route to the file

*/
int get_ResultStep(char *);

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
  \fn Matrix Read_Delimited_File__InOutLib__(char * Name_File)

  \brief Read a file with a given format
  \param Name_File : file name
*/
Matrix Read_Delimited_File__InOutLib__(char *);

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

/*! 

  \fn void print_convergence_stats(int TimeStep,int Iter, double Error_total, double Error_relative)
 
  \brief Print number of iterations required.
 
  \param Time : Current step of the simulation 
  \param Iter : Number of iterations required
  \param Error_total 
  \param Error_relative

 */
void print_convergence_stats(int, int, double, double, double);

/*****************************************************************/

/*!

  \fn Mesh GramsBox(char * File,Time_Int_Params Parameters_Solver)
  
  \brief Function devoted to generate the computational grid of nodes
  
  \param File: Name of the file with the instructions

*/
Mesh GramsBox(char *,Time_Int_Params);
/*****************************************************************/

/*!

  \fn Boundaries Read_u_Dirichlet_Boundary_Conditions__InOutFun__(char * File,int NumBounds,int NumTimeSteps)

  \brief Function that reads the boundary conditions for the mesh

  \param File: Name of the file with the instructions
  \param NumBounds : Numerb of boundary defined in the domain
  \param NumTimeSteps

*/
Boundaries Read_u_Dirichlet_Boundary_Conditions__InOutFun__(char *,int,int);
/*****************************************************************/

/*!

  \fn Boundaries Read_u_Neumann_Boundary_Conditions__InOutFun__(char * Name_File,int NumBounds,int GPxElement, int NumTimeSteps)
  
  \brief Generate the load state

  \param File: Name of the file with the instructions
  \param NumNeumannBC : Number of loads
  \param GPxElement : As the particle discretization is performed thorouht 
  \param NumTimeSteps
*/
Boundaries Read_u_Neumann_Boundary_Conditions__InOutFun__(char *,int,int,int);

/*****************************************************************/

/*!
  \fn void Check_u_Neumann_Boundary_Conditions__InOutFun__(Boundaries Neumann_Contours,int NumParticles)
  
  \brief Check if the Neumann boundary condition is well applied

  \param Neumann_Contours : List with the properties of the boundary condition
  \param NumParticles : Number of particles in the all simulation
*/
void Check_u_Neumann_Boundary_Conditions__InOutFun__(Boundaries,int);

/*****************************************************************/


/*!

  \fn Boundaries Read_upw_Dirichlet_Boundary_Conditions__InOutFun__(char * Name_File,int NumBounds, int NumTimeSteps)

  \brief Function that reads the boundary conditions for the mesh

  \param File: Name of the file with the instructions
  \param NumBounds : Numerb of boundary defined in the domain
  \param NumTimeSteps : 

*/
Boundaries Read_upw_Dirichlet_Boundary_Conditions__InOutFun__(char *,int,int);
/*****************************************************************/

/*!
  \fn Boundaries Read_upw_Neumann_Boundary_Conditions__InOutFun__(char * Name_File,int NumBounds,int GPxElement, int NumTimeSteps)

  \brief Function that reads the neumann boundary conditions for the mesh

  \param File: Name of the file with the instructions
  \param NumBounds : Numerb of boundary defined in the domain
  \param GPxElement
  \param NumTimeSteps:
*/
Boundaries Read_upw_Neumann_Boundary_Conditions__InOutFun__(char *,int,int,int);
/*****************************************************************/


Particle Generate_Gauss_Point_Analysis__InOutFun__(char *);
/*****************************************************************/

/*!
  \fn Particle Generate_Soil_Water_Coupling_Analysis__InOutFun__(char * Name_File, Mesh FEM_Mesh, Time_Int_Params Parameters_Solver);
 
  \brief Function that generate a set of material points for soil-water coupling applications

  \param File: Name of the file with the instructions
  \param Nodes : Set of background nodes
  \param Parameters_Solver
*/
Particle Generate_Soil_Water_Coupling_Analysis__InOutFun__(char *, Mesh, Time_Int_Params);
/*****************************************************************/

/*!
  \fn void Solver_selector__InOutFun__(char * File)

  \brief Function to define the solver
  
  \param File: Name of the file with the instructions
}
*/
Time_Int_Params Solver_selector__InOutFun__(char *);
/*****************************************************************/

/*!
  \fn Material * GramsMaterials(char * File, Particle Particles, int GPxElement)

  \brief Generate the libreary of materials

  \param File: Name of the file with the instructions
  \param Particles : Particle discretization
  \param GPxElement : As the particle discretization is performed thorouht 
*/
Material * GramsMaterials(char *, Particle, int);
/*****************************************************************/

/*!
  \fn Material * Read_Materials__InOutFun__(char * SimulationFile, int NumberMaterials);
*/
Material * Read_Materials__InOutFun__(char *, int);
/*****************************************************************/

/*!
  \fn Particle Generate_One_Phase_Analysis__InOutFun__(char * Name_File, Mesh FEM_Mesh, Time_Int_Params Parameters_Solver)
*/
Particle Generate_One_Phase_Analysis__InOutFun__(char *, Mesh,Time_Int_Params);
/*****************************************************************/

/*!
\fn Mixture * Read_Soil_Water_Mixtures__InOutFun__(char * SimulationFile, int Number_Soil_Water_Mixtures);

*/
Mixture * Read_Soil_Water_Mixtures__InOutFun__(char *, int);
/*****************************************************************/

/*!
  
  \fn void Initial_condition_particles__InOutFun__(char * File, Particle Particles, int GPxElement)
  
  \brief Define initial conditions for the particles
  
  \param File: Name of the file with the instructions
  \param Particles : Particle discretization
  \param GPxElement : As the particle discretization is performed thorouht 

*/
void Initial_condition_particles__InOutFun__(char *, Particle, int);
/*****************************************************************/

/*!
  
  \fn void Initial_condition_nodes__InOutFun__(char * File, Particle Particles, Mesh FEM_Mesh)
  
  \brief Define initial conditions for particles using nodal values
  
  \param File: Name of the file with the instructions
  \param Particles : Particle discretization
  \param FEM_Mesh
*/
void Initial_condition_nodes__InOutFun__(char *, Particle, Mesh);
/*****************************************************************/

/*!

  \fn void GramsShapeFun(char * File)

  \brief Initialize the shape functions. 

  Example : 
  GramsShapeFun (Type=MPMQ4) {
  }
  GramsShapeFun (Type=uGIMP) {
  }
  GramsShapeFun (Type=LME) {
  gamma=2.3
  TOL_lambda=10e-6
  }

  \param File: Name of the file with the instructions

*/
void GramsShapeFun(char * );
/*****************************************************************/

/*!

  \fn void GramsOutputs(char * File)

  \brief Define the output directory and the interval

  \param File: Name of the file with the instructions

*/
void GramsOutputs(char * );
/*****************************************************************/

/*!
  \fn void OutputSimulation(Particle Set_Particles,int Numerical_T,double Physical_T,
  double DeltaT,Event Parameters)
 
  \brief Generate backup copy
  
  \param Set_Particles
  \param Numerical_T
  \param Physical_T
  \param DeltaT
  \param Parameters
*/
void OutputSimulation(Particle,int,double,double,Event);
/*****************************************************************/

/*!

  \fn Load * GramsBodyForces(char * File, int NumBodyForces, int GPxElement, int NumTimeSteps)
  
  \brief Generate the load state

  \param File: Name of the file with the instructions
  \param NumBodyForces : Number of body forces
  \param GPxElement : As the particle discretization is performed thorouht 
  \param NumTimeSteps:
*/
Load * GramsBodyForces(char *, int, int, int);

/*****************************************************************/

void particle_results_vtk__InOutFun__(Particle, int, int);

/*****************************************************************/

void particle_backup_vtk__InOutFun__(Particle,int,int);

/*****************************************************************/

void nodal_results_vtk__InOutFun__(Mesh, Mask, Matrix, int, int);

/*****************************************************************/

void WriteGnuplot(Matrix, Matrix, double, double, int, int, char[20]);

/*****************************************************************/

/*!
\fn Particle restart_Simulation(char * Parameters,char * Restart,Mesh Nodes)

\brief Generate simulation structure from previous results

\param Parameters
\param Restart
\param Nodes
*/
Particle restart_Simulation(char *,char *,Mesh);

/*****************************************************************/

void path_nodes_analysis_csv__InOutFun__(Matrix, Matrix, char *, Mask, Event, int, int, double);
void path_particles_analysis_csv__InOutFun__(Matrix, Matrix, char *, Event, int, int, double);
/*****************************************************************/


void NLPS_Out_nodal_path_csv__InOutFun__(char * Name_File);
void NLPS_Out_particles_path_csv__InOutFun__(char * Name_File);
/*****************************************************************/


void Gauss_Point_evolution__InOutFun__(Particle, Event, char *, int, int);
/*****************************************************************/


int Hidrostatic_condition_particles__InOutFun__(char *, Particle, int);
/*****************************************************************/

Material Define_Solid_Rigid(FILE *,char *,int);
Material Define_Linear_Elastic(FILE *,char *,int);
Material Define_Saint_Venant_Kirchhoff(FILE *,char *,int);
Material Define_Neo_Hookean_Wriggers(FILE *,char *,int);
Material Define_Von_Mises(FILE *,char *,int);

Material Define_Frictional(FILE *,char *,int);
void Initialize_Frictional(double *,double *,Material);


int Define_Drucker_Prager(
  Material *DP_Material /**< [out] List with the material properties */,
  FILE *Simulation_file /**< [in] Simulation file */,
  char *Material_Model /**< [in] Character identifier for the material */, 
  int Material_Idx /**< [in] Index identifier for the material */);
/*****************************************************************/

Material Define_Compressible_Newtonian_Fluid(FILE *,char *,int);
Material Define_Incompressible_Newtonian_Fluid(FILE *,char *,int);
/*****************************************************************/

#endif
