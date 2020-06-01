#include "nl-partsol.h"

/*
  List of auxiliar functions for the main function : restart_Simulation
*/
static Fields restart_Fields(char *);
static int    restart_ReadVtk_Points(char *);
static Matrix restart_ReadVtk_Scalars(char *,char *);
static int *  restart_ReadVtk_I0(char *,char *);
static int *  restart_ReadVtk_MatIdx(char *,char *);
static Matrix restart_ReadVtk_Vectors(char *,char *);
static Matrix restart_ReadVtk_Tensors(char *,char *);

/**********************************************************************/

GaussPoint restart_Simulation(char * File_Parameters,
			      char * File_Restart,
			      Mesh FEM_Mesh)
{

  int Ndim = NumberDimensions;
  int Np;

  /* Simulation file */
  FILE * Sim_dat;
 
  GaussPoint Set_Particles;

  /* Set to false check variables */
  bool Is_GramsSolid2D = false;
  bool Is_GramsShapeFun = false;
  bool Is_GramsMaterials = false;
  bool Is_GramsInitials = false;
  bool Is_GramsBodyForces = false;
  bool Is_GramsNeumannBC = false;

  /* Initialize counters */
  int Counter_Materials = 0;
  int Counter_BodyForces = 0;
  int Counter_GramsNeumannBC = 0;

  while(fgets(Line_GramsSolid2D,sizeof(Line_GramsSolid2D),Sim_dat) != NULL)
    {
    
      /* Read the line with the space as separators */
      Num_words_parse = parse(Parse_GramsSolid2D,Line_GramsSolid2D," \n\t");
      if (Num_words_parse < 0)
	{
	  fprintf(stderr,"%s : %s \n",
		  "Error in GramsSolid2D()",
		  "Parser failed");
	  exit(EXIT_FAILURE);
	}
    
      if ((Num_words_parse > 0) &&
	  (strcmp(Parse_GramsSolid2D[0],"GramsSolid2D") == 0))
	{
	  Is_GramsSolid2D = true;    
	}
      if ((Num_words_parse > 0) &&
	  (strcmp(Parse_GramsSolid2D[0],"GramsShapeFun") == 0))
	{
	  Is_GramsShapeFun = true;
	}
      if ((Num_words_parse > 0) &&
	  (strcmp(Parse_GramsSolid2D[0],"GramsMaterials") == 0))
	{
	  Is_GramsMaterials = true;
	  Counter_Materials++;
	}
      if ((Num_words_parse > 0) &&
	  (strcmp(Parse_GramsSolid2D[0],"GramsInitials") == 0))
	{
	  Is_GramsInitials = true;
	}
      if ((Num_words_parse > 0) &&
	  (strcmp(Parse_GramsSolid2D[0],"GramsBodyForces") == 0))
	{
	  Is_GramsBodyForces = true;
	  Counter_BodyForces++;
	}
      if ((Num_words_parse > 0) &&
	  (strcmp(Parse_GramsSolid2D[0],"GramsNeumannBC") == 0))
	{
	  Is_GramsNeumannBC = true;
	  Counter_GramsNeumannBC++;
	}    
    }
  
  /* Read particle mesh properties */
  if(Is_GramsSolid2D)
    {    
      /* Read the number of particles */
      Np = restart_ReadVtk_Points(File_Restart);
      Set_Particles.NumGP = Np;
      
      /* Closest node to the particle */
      Set_Particles.I0 = restart_ReadVtk_I0(File_Restart,Np);
      
      /* Number of tributary nodes for each particle */
      Set_Particles.NumberNodes = (int *)Allocate_ArrayZ(Np,sizeof(int));
      
      /* Tributary nodes for each particle */
      Set_Particles.ListNodes =  allocate_SetTable(Np);
      
      /* List of particles close to each particle */
      Set_Particles.Beps =  allocate_SetTable(Np);
      
      /* Read list of fields */
      Set_Particles.Phi = restart_Fields(File_Restart,Np);
    }
  
  /* Read Material parameters */
  if(Is_GramsMaterials)
    {
      Set_Particles.NumberMaterials = Counter_Materials;
      Set_Particles.MatIdx = (int *)malloc(NumParticles*sizeof(int));
      Set_Particles.Mat = GramsMaterials(Name_File,MPM_Mesh,GPxElement);
    }
  
  /* Read external forces */
  if(Is_GramsNeumannBC)
    {
      MPM_Mesh.NumNeumannBC = Counter_GramsNeumannBC;
      MPM_Mesh.F = GramsNeumannBC(Name_File,Counter_GramsNeumannBC,GPxElement);
    }
  
  /* Read body forces */    
  if(Is_GramsBodyForces)
    {
      MPM_Mesh.NumberBodyForces = Counter_BodyForces;
      MPM_Mesh.B = GramsBodyForces(Name_File,Counter_BodyForces,GPxElement); 
    }
  
  /* Initialize shape functions */
  if(strcmp(ShapeFunctionGP,"MPMQ4") == 0)
    {
      Q4_Initialize(MPM_Mesh, FEM_Mesh);
    }
  else if(strcmp(ShapeFunctionGP,"uGIMP") == 0)
    {
      Set_Particles.lp = MatAllocZ(NumParticles,Ndim);
      strcpy(Set_Particles.lp.Info,"Voxel lenght GP");
      uGIMP_Initialize(MPM_Mesh,FEM_Mesh);
    }
  else if(strcmp(ShapeFunctionGP,"LME") == 0)
    {
      Set_Particles.lambda = MatAllocZ(NumParticles,Ndim);
      strcpy(Set_Particles.lambda.Info,"Lagrange Multiplier");
      Set_Particles.Beta = MatAllocZ(NumParticles,Ndim);
      strcpy(Set_Particles.Beta.Info,"Beta parameter");
      LME_Initialize(MPM_Mesh,FEM_Mesh);
    } 
  
  return Set_Particles;
}

/**********************************************************************/

static Fields restart_Fields(char * File_Restart,int Np)
{
  Fields Phi;
  
  /* Density field */
  Phi.rho = restart_ReadVtk_Scalars(File_Restart,"DENSITY",Np);
  
  /* Mass field */
  Phi.mass = restart_ReadVtk_Scalars(File_Restart,"MASS",Np);
  
  /* Deformation Energy */
  Phi.W = restart_ReadVtk_Scalars(File_Restart,"W",Np);

  /* Strain during crack */
  Phi.Strain_If = restart_ReadVtk_Scalars(File_Restart,"STRAIN_IF",Np);
  
  /* Damage parameter (Fracture) */
  Phi.ji = restart_ReadVtk_Scalars(File_Restart,"Ji",Np);

  /* Position in global coordinates */
  Phi.x_GC = restart_ReadVtk_Vectors(File_Restart,"X_GC",Np); 
  
  /* Position in element coordiantes */
  Phi.x_EC = restart_ReadVtk_Vectors(File_Restart,"X_EC",Np); 
 
  /* Displacement field */
  Phi.dis = restart_ReadVtk_Vectors(File_Restart,"DISPLACEMENT",Np);
  
  /* Velocity field */
  Phi.vel = restart_ReadVtk_Vectors(File_Restart,"VELOCITY",Np);
    
  /* Acceleration field */
  Phi.acc = restart_ReadVtk_Vectors(File_Restart,"ACCELERATION",Np);
  
  /* Stress field */
  Phi.Stress = restart_ReadVtk_Tensors(File_Restart,"STRESS",Np);
  
  /* Strain field */
  Phi.Strain = restart_ReadVtk_Tensors(File_Restart,"STRAIN",Np);
  
  return Phi;
}

/**********************************************************************/

/*!
  \fn static Matrix restart_ReadVtk_Points(char * File_Parameters)

  \brief Read the number of particles from a .vtk file

*/
static int restart_ReadVtk_Points(char * Simulation_File)
{
  int Np = 0;
  bool Is_POINTS = false;
  
  /* Open and check file */
  FILE * Simulation_Ptr = fopen(Simulation_File,"r");  
  if (Sim_dat==NULL){
    fprintf(stderr,"%s : \n\t %s %s",
	    "Error in GramsInitials()","Incorrect lecture of",
	    Name_File);
    exit(EXIT_FAILURE);
  }

  /* Variable to read a file */
  char Line[MAXC] = {0}; 
  char * Parse_Line[MAXW] = {NULL};
  int Num_Words;

  /* Search "POINTS" keyword in the .vtk file */
  while(fgets(Line,sizeof(Line),Simulation_Ptr) != NULL)
    {
      Num_Words = parse(Parse_Line,Line," \n\t");

      if ((Num_words == 3) &&
	  (strcmp(Parse_Line[0],"POINTS") == 0) &&
	  (strcmp(Parse_Line[2],"float") == 0))
	{
	  Np = atoi(Parse_Line[1]);
	  Is_POINTS = true;
	  break;
	}
      
    }
  
  /* If there is not POINTS keyword in the file */
  if(!Is_POINTS)
    {
      fprintf(stderr,"Error in restart_ReadVtk_Points(%s) : %s \n",
	      Simulation_File,"Not keyword POINTS");
      exit(EXIT_FAILURE);
    }
  
  return Np;

}

/**********************************************************************/

static Matrix restart_ReadVtk_Scalars(char * Simulation_File,char * Parameter,int Np)
{

  bool Is_SCALAR = false;
  Matrix Scalar;

  /* Open and check file */
  FILE * Simulation_Ptr = fopen(Simulation_File,"r");  
  if (Sim_dat==NULL){
    fprintf(stderr,"%s : \n\t %s %s",
	    "Error in GramsInitials()","Incorrect lecture of",
	    Name_File);
    exit(EXIT_FAILURE);
  }

  /* Variable to read a file */
  char Line[MAXC] = {0}; 
  char * Parse_Line[MAXW] = {NULL};
  int Num_Words;
  
  /* Search "Parameter" keyword in the .vtk file */
  while(fgets(Line,sizeof(Line),Simulation_Ptr) != NULL)
    {
      Num_Words = parse(Parse_Line,Line," \n\t");

      if ((Num_words == 3) &&
	  (strcmp(Parse_Line[0],"SCALARS") == 0) &&
	  (strcmp(Parse_Line[1],Parameter) == 0) &&
	  (strcmp(Parse_Line[2],"float") == 0))
	{
	  Is_SCALAR = true;
	  break;
	}
      
    }
  
  /* Read scalar variable */
  if(Is_SCALAR)
    {
      Scalar = MatAllocZ(Np,1);

      for(int p = -1 ; p<Np ; p++)
	{

	  fgets(Line,sizeof(Line),Simulation_Ptr);
	  Num_Words = parse(Parse_Line,Line," \n\t");
	
	  if(p >= 0)
	    {
	      Scalar.nV[p] = atof(Parse_Line[0]);
	    }
	
	}
      
    }

  return Scalar; 
}

/**********************************************************************/

static int * restart_ReadVtk_I0(char * Simulation_File,int Np)
{

  bool Is_I0 = false;
  Matrix I0;
  
  /* Open and check file */
  FILE * Simulation_Ptr = fopen(Simulation_File,"r");  
  if (Sim_dat==NULL){
    fprintf(stderr,"%s : \n\t %s %s",
	    "Error in GramsInitials()","Incorrect lecture of",
	    Name_File);
    exit(EXIT_FAILURE);
  }

  /* Variable to read a file */
  char Line[MAXC] = {0}; 
  char * Parse_Line[MAXW] = {NULL};
  int Num_Words;
  
  /* Search "ELEM_i" keyword in the .vtk file */
  while(fgets(Line,sizeof(Line),Simulation_Ptr) != NULL)
    {
      Num_Words = parse(Parse_Line,Line," \n\t");

      if ((Num_words == 3) &&
	  (strcmp(Parse_Line[0],"SCALARS") == 0) &&
	  (strcmp(Parse_Line[1],"ELEM_i") == 0) &&
	  (strcmp(Parse_Line[2],"integer") == 0))
	{
	  Is_I0 = true;
	  break;
	}
      
    }
  
  /* Read scalar variable */
  if(Is_I0)
    {
      Scalar = (int *)Allocate_ArrayZ(Np,sizeof(int));

      for(int p = -1 ; p<Np ; p++)
	{

	  fgets(Line,sizeof(Line),Simulation_Ptr);
	  Num_Words = parse(Parse_Line,Line," \n\t");
	
	  if(p >= 0)
	    {
	      I0[p] = atoi(Parse_Line[0]);
	    }
	
	}
      
    }

  return I0; 
}

/**********************************************************************/

static int * restart_ReadVtk_MatIdx(char * File,char * Parameter,int Np)
{

  bool Is_MatIdx = false;
  Matrix MatIdx;

  /* Open and check file */
  FILE * Simulation_Ptr = fopen(Simulation_File,"r");  
  if (Sim_dat==NULL){
    fprintf(stderr,"%s : \n\t %s %s",
	    "Error in GramsInitials()","Incorrect lecture of",
	    Name_File);
    exit(EXIT_FAILURE);
  }

  /* Variable to read a file */
  char Line[MAXC] = {0}; 
  char * Parse_Line[MAXW] = {NULL};
  int Num_Words;
  
  /* Search "MatIdx" keyword in the .vtk file */
  while(fgets(Line,sizeof(Line),Simulation_Ptr) != NULL)
    {
      Num_Words = parse(Parse_Line,Line," \n\t");

      if ((Num_words == 3) &&
	  (strcmp(Parse_Line[0],"SCALARS") == 0) &&
	  (strcmp(Parse_Line[1],"MatIdx") == 0) &&
	  (strcmp(Parse_Line[2],"integer") == 0))
	{
	  Is_MatIdx = true;
	  break;
	}
      
    }
  
  /* Read scalar variable */
  if(Is_MatIdx)
    {
      Scalar = (int *)Allocate_ArrayZ(Np,sizeof(int));

      for(int p = -1 ; p<Np ; p++)
	{

	  fgets(Line,sizeof(Line),Simulation_Ptr);
	  Num_Words = parse(Parse_Line,Line," \n\t");
	
	  if(p >= 0)
	    {
	      MatIdx[p] = atoi(Parse_Line[0]);
	    }
	
	}
      
    }

  return MatIdx; 
}

/**********************************************************************/

static Matrix restart_ReadVtk_Vectors(char * File,char * Parameter,int Np)
{
  
  int Ndim = NumberDimensions;
  bool Is_Vector = false;
  Matrix Vector;

  /* Open and check file */
  FILE * Simulation_Ptr = fopen(Simulation_File,"r");  
  if (Sim_dat==NULL){
    fprintf(stderr,"%s : \n\t %s %s",
	    "Error in GramsInitials()","Incorrect lecture of",
	    Name_File);
    exit(EXIT_FAILURE);
  }

  /* Variable to read a file */
  char Line[MAXC] = {0}; 
  char * Parse_Line[MAXW] = {NULL};
  int Num_Words;

  /* Search "MatIdx" keyword in the .vtk file */
  while(fgets(Line,sizeof(Line),Simulation_Ptr) != NULL)
    {
      Num_Words = parse(Parse_Line,Line," \n\t");

      if ((Num_words == 3) &&
	  (strcmp(Parse_Line[0],"VECTORS") == 0) &&
	  (strcmp(Parse_Line[1],Parameter) == 0) &&
	  (strcmp(Parse_Line[2],"float") == 0))
	{
	  Is_Vector = true;
	  break;
	}
      
    }
    
  /* Generate the matrix and fill it */
  if(Is_Vector)
    {
      /* Allocate coordinates */
      Vector = MatAllocZ(Np,Ndim);

      /* Fill coordinates matrix */
      for(int p = 0 ; p<Np ; p++)
	{
	  /* Read each line and parse it */
	  fgets(Line,sizeof(Line),Simulation_Ptr);
	  Num_Words = parse(Parse_Line,Line," \n\t");

	  /* Fill coordinates */
	  for(int i = 0 ; i<Ndim ; i++)
	    {
	      Vector.nM[p][i] = atof(Parse_Line[i]);
	    }
	}

      /* Close the file */
      fclose(Simulation_Ptr);

    }
}

/**********************************************************************/

static Matrix restart_ReadVtk_Tensors(char * File,char * Parameter)
{

  int Ndim = NumberDimensions;
  bool Is_Tensor = false;
  Matrix Tensor;

  /* Open and check file */
  FILE * Simulation_Ptr = fopen(Simulation_File,"r");  
  if (Sim_dat==NULL){
    fprintf(stderr,"%s : \n\t %s %s",
	    "Error in GramsInitials()","Incorrect lecture of",
	    Name_File);
    exit(EXIT_FAILURE);
  }

  /* Variable to read a file */
  char Line[MAXC] = {0}; 
  char * Parse_Line[MAXW] = {NULL};
  int Num_Words;

  /* Search "MatIdx" keyword in the .vtk file */
  while(fgets(Line,sizeof(Line),Simulation_Ptr) != NULL)
    {
      Num_Words = parse(Parse_Line,Line," \n\t");

      if ((Num_words == 3) &&
	  (strcmp(Parse_Line[0],"TENSORS") == 0) &&
	  (strcmp(Parse_Line[1],Parameter) == 0) &&
	  (strcmp(Parse_Line[2],"float") == 0))
	{
	  Is_Tensor = true;
	  break;
	}
      
    }

  /* Generate the matrix and fill it */
  if(Is_Tensor)
    {
      /* Allocate coordinates */
      Tensor = MatAllocZ(Np,Ndim*Ndim);

      /* Fill coordinates matrix */
      for(int p = 0 ; p<Np ; p++)
	{
	  /* Fill coordinates */
	  for(int i = 0 ; i<Ndim ; i++)
	    {	      
	      /* Read each line and parse it */
	      fgets(Line,sizeof(Line),Simulation_Ptr);
	      Num_Words = parse(Parse_Line,Line," \n\t");

	      for(int j = 0 ; j<Ndim ; j++){	      
		Tensor.nM[p][i*j] = atof(Parse_Line[j]);
	      }
	    }
	  
	  /* Read white line at the end of each tensor */
	  fgets(Line,sizeof(Line),Simulation_Ptr);
	}

      /* Close the file */
      fclose(Simulation_Ptr);

    }
  
}

/**********************************************************************/
