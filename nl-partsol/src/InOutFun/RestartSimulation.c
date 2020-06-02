#include "nl-partsol.h"

/*
  List of auxiliar functions for the main function : restart_Simulation
*/
static Fields restart_Fields(char *,int);
static int    restart_ReadVtk_Points(char *);
static Matrix restart_ReadVtk_Scalars(char *,char *,int);
static int *  restart_ReadVtk_I0(char *,int);
static Matrix restart_ReadVtk_Vectors(char *,char *,int);
static Matrix restart_ReadVtk_Tensors(char *,char *,int);

/**********************************************************************/

GaussPoint restart_Simulation(char * File_Parameters,
			      char * File_Restart,
			      Mesh FEM_Mesh)
{

  int Ndim = NumberDimensions;
  int Np;

  GaussPoint Set_Particles;

  /* Set to false check variables */
  bool Is_GramsSolid2D = false;
  bool Is_GramsShapeFun = false;
  bool Is_GramsMaterials = false;
  bool Is_GramsInitials = false;
  bool Is_GramsBodyForces = false;
  bool Is_GramsNeumannBC = false;

  int GPxElement = 1;

  /* Initialize counters */
  int Counter_Materials = 0;
  int Counter_BodyForces = 0;
  int Counter_GramsNeumannBC = 0;

  /* Open and check file */
  FILE * Simulation_Ptr = fopen(File_Parameters,"r");  
  if (Simulation_Ptr == NULL)
    {
      fprintf(stderr,"%s : \n\t %s %s",
	      "Error in GramsInitials()",
	      "Incorrect lecture of",File_Parameters);
      exit(EXIT_FAILURE);
    }

  /* Variable to read a file */
  char Line[MAXC] = {0}; 
  char * Parse_Line[MAXW] = {NULL};
  int Num_Words;

  while(fgets(Line,sizeof(Line),Simulation_Ptr) != NULL)
    {
    
      /* Read the line with the space as separators */
      Num_Words = parse(Parse_Line,Line," \n\t");
      if (Num_Words < 0)
	{
	  fprintf(stderr,"%s : %s \n",
		  "Error in GramsSolid2D()",
		  "Parser failed");
	  exit(EXIT_FAILURE);
	}
    
      if ((Num_Words > 0) &&
	  (strcmp(Parse_Line[0],"GramsSolid2D") == 0))
	{
	  Is_GramsSolid2D = true;    
	}
      if ((Num_Words > 0) &&
	  (strcmp(Parse_Line[0],"GramsShapeFun") == 0))
	{
	  Is_GramsShapeFun = true;
	}
      if ((Num_Words > 0) &&
	  (strcmp(Parse_Line[0],"GramsMaterials") == 0))
	{
	  Is_GramsMaterials = true;
	  Counter_Materials++;
	}
      if ((Num_Words > 0) &&
	  (strcmp(Parse_Line[0],"GramsInitials") == 0))
	{
	  Is_GramsInitials = true;
	}
      if ((Num_Words > 0) &&
	  (strcmp(Parse_Line[0],"GramsBodyForces") == 0))
	{
	  Is_GramsBodyForces = true;
	  Counter_BodyForces++;
	}
      if ((Num_Words > 0) &&
	  (strcmp(Parse_Line[0],"GramsNeumannBC") == 0))
	{
	  Is_GramsNeumannBC = true;
	  Counter_GramsNeumannBC++;
	}    
    }

  /* Close the file */
  fclose(Simulation_Ptr);
  
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
      Set_Particles.ListNodes = allocate_SetTable(Np);
      
      /* List of particles close to each particle */
      Set_Particles.Beps =  allocate_SetTable(Np);
      
      /* Read list of fields */
      Set_Particles.Phi = restart_Fields(File_Restart,Np);
    }

  
  if(Is_GramsMaterials)
    {
      puts("*************************************************");
      printf(" \t %s \n","* Read materials properties");
      Set_Particles.NumberMaterials = Counter_Materials;
      Set_Particles.MatIdx = (int *)malloc(Np*sizeof(int));
      Set_Particles.Mat = GramsMaterials(File_Parameters,Set_Particles,GPxElement);
    }

  /* Read external forces */
  if(Is_GramsNeumannBC)
    {
      puts("*************************************************");
      printf(" \t %s \n","* Read Neumann boundary conditions");
      Set_Particles.NumNeumannBC = Counter_GramsNeumannBC;
      Set_Particles.F = GramsNeumannBC(File_Parameters,
				       Counter_GramsNeumannBC,
				       GPxElement);
    }
  
  /* Read body forces */    
  if(Is_GramsBodyForces)
    {
      puts("*************************************************");
      printf(" \t %s \n","* Read body forces");
      Set_Particles.NumberBodyForces = Counter_BodyForces;
      Set_Particles.B = GramsBodyForces(File_Parameters,
					Counter_BodyForces,
					GPxElement); 
    }
  
  /* Initialize shape functions */
  puts("*************************************************");
  printf(" \t %s \n","* Read shape functions");
  GramsShapeFun(File_Parameters);
  
  if(strcmp(ShapeFunctionGP,"MPMQ4") == 0)
    {
      Q4_Initialize(Set_Particles, FEM_Mesh);
    }
  else if(strcmp(ShapeFunctionGP,"uGIMP") == 0)
    {
      Set_Particles.lp = MatAllocZ(Np,Ndim);
      strcpy(Set_Particles.lp.Info,"Voxel lenght GP");
      uGIMP_Initialize(Set_Particles,FEM_Mesh);
    }
  else if(strcmp(ShapeFunctionGP,"LME") == 0)
    {
      Set_Particles.lambda = MatAllocZ(Np,Ndim);
      strcpy(Set_Particles.lambda.Info,"Lagrange Multiplier");
      Set_Particles.Beta = MatAllocZ(Np,Ndim);
      strcpy(Set_Particles.Beta.Info,"Beta parameter");
      LME_Initialize(Set_Particles,FEM_Mesh);
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
static int restart_ReadVtk_Points(char * File)
{
  int Np = 0;
  bool Is_POINTS = false;
  
  /* Open and check file */
  FILE * Simulation_Ptr = fopen(File,"r");  
  if (Simulation_Ptr == NULL)
    {
      fprintf(stderr,"%s : \n\t %s %s",
	      "Error in GramsInitials()","Incorrect lecture of",
	      File);
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

      if ((Num_Words == 3) &&
	  (strcmp(Parse_Line[0],"POINTS") == 0) &&
	  (strcmp(Parse_Line[2],"float") == 0))
	{
	  Np = atoi(Parse_Line[1]);
	  Is_POINTS = true;
	  break;
	}
      
    }
  
  /* If there is not POINTS keyword in the file */
  if(Is_POINTS)
    {
      /* Close the file */
      fclose(Simulation_Ptr);
      /* Number of particles */
      return Np;
    }
  else
    {
      /* Close the file */
      fclose(Simulation_Ptr);
      /* Error message */      
      fprintf(stderr,"Error in restart_ReadVtk_Points(%s) : %s \n",
	      File,"Not keyword POINTS");
      exit(EXIT_FAILURE);
    }
}

/**********************************************************************/

static Matrix restart_ReadVtk_Scalars(char * File,char * Parameter,int Np)
{

  bool Is_SCALAR = false;
  Matrix Scalar;

  /* Open and check file */
  FILE * Simulation_Ptr = fopen(File,"r");  
  if (Simulation_Ptr == NULL)
    {
      fprintf(stderr,"%s : \n\t %s %s",
	      "Error in GramsInitials()","Incorrect lecture of",
	      File);
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

      if ((Num_Words == 3) &&
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

  /* Close the file */
  fclose(Simulation_Ptr);


  return Scalar; 
}

/**********************************************************************/

static int * restart_ReadVtk_I0(char * File,int Np)
{

  bool Is_I0 = false;
  int * I0;
  
  /* Open and check file */
  FILE * Simulation_Ptr = fopen(File,"r");  
  if (Simulation_Ptr == NULL)
    {
      fprintf(stderr,"%s : \n\t %s %s",
	    "Error in GramsInitials()","Incorrect lecture of",
	      File);
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

      if ((Num_Words == 3) &&
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
      I0 = (int *)Allocate_ArrayZ(Np,sizeof(int));

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

  /* Close the file */
  fclose(Simulation_Ptr);

  return I0; 
}

/**********************************************************************/

static Matrix restart_ReadVtk_Vectors(char * File,char * Parameter,int Np)
{
  
  int Ndim = NumberDimensions;
  bool Is_Vector = false;
  Matrix Vector_Field;

  /* Open and check file */
  FILE * Simulation_Ptr = fopen(File,"r");  
  if (Simulation_Ptr == NULL)
    {
      fprintf(stderr,"%s : \n\t %s %s",
	      "Error in GramsInitials()","Incorrect lecture of",File);
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

      if ((Num_Words == 3) &&
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
      Vector_Field = MatAllocZ(Np,Ndim);

      /* Fill coordinates matrix */
      for(int p = 0 ; p<Np ; p++)
	{
	  /* Read each line and parse it */
	  fgets(Line,sizeof(Line),Simulation_Ptr);
	  Num_Words = parse(Parse_Line,Line," \n\t");

	  /* Fill coordinates */
	  for(int i = 0 ; i<Ndim ; i++)
	    {
	      Vector_Field.nM[p][i] = atof(Parse_Line[i]);
	    }
	}
    }

  /* Close the file */
  fclose(Simulation_Ptr);
  
  return Vector_Field;  
}

/**********************************************************************/

static Matrix restart_ReadVtk_Tensors(char * File,char * Parameter,int Np)
{

  int Ndim = NumberDimensions;
  bool Is_Tensor = false;
  Matrix Tensor_Field;

  /* Open and check file */
  FILE * Simulation_Ptr = fopen(File,"r");  
  if (Simulation_Ptr == NULL)
    {
      fprintf(stderr,"%s : \n\t %s %s",
	      "Error in GramsInitials()","Incorrect lecture of",File);
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

      if ((Num_Words == 3) &&
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
      Tensor_Field = MatAllocZ(Np,Ndim*Ndim);

      /* Fill coordinates matrix */
      for(int p = 0 ; p<Np ; p++)
	{
	  /* Fill coordinates */
	  for(int i = 0 ; i<Ndim ; i++)
	    {	      
	      /* Read each line and parse it */
	      fgets(Line,sizeof(Line),Simulation_Ptr);
	      Num_Words = parse(Parse_Line,Line," \n\t");

	      for(int j = 0 ; j<Ndim ; j++)
		{	      
		  Tensor_Field.nM[p][i*j] = atof(Parse_Line[j]);
		}
	    }
	  
	  /* Read white line at the end of each tensor */
	  fgets(Line,sizeof(Line),Simulation_Ptr);
	}
    }

  /* Close the file */
  fclose(Simulation_Ptr);
  
  return Tensor_Field;  
}

/**********************************************************************/
