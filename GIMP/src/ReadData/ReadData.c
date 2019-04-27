#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "../ToolsLib/GlobalVariables.h"
#include "../ToolsLib/Utils.h"

/***************************************************************************/


void ReadData(char * Name_File)
/*
  Read data from the .DAT file and initialize the variables
  
  Inputs
  - Name_file : Name of the file
  - 
  Outputs
  - Conectivity matrix
  - Coordenates of the nodes
  - phi_n : Values of the variables in the n step, 
        initializerd with the Initial conditions
  - DeltaT : Time-step
  - A_el : Area of the element
  - type_elem : Type of the element (1)
  - 
  - N_nodes : Number of nodes
  - N_elem : Number of elements
  - N_steps : Number of time steps
*/
{
  /* Number of lines with coments in the header */
  int num_coment;

  /* Number of element in the line , just for check */
  int aux_read;

  /* Array with the possible index fields */
  int AuxiliarArrayFields[14];

  /* Variable for the index in nodes and conectivitie */
  int ix;

  /* Auxiliar variable for reading the lines in the files */
  char line[200];
  
  /* Allocate a string with enough space for the extensions. */
  char *Name_Simulation = malloc(strlen(Name_File+5));
  
  /* Copy the name with extensions into fn. */
  sprintf(Name_Simulation, "%s.dat", Name_File); 
  
  /* Open and check .dat file */
  Sim_dat = fopen(Name_Simulation,"r");  
  if (Sim_dat==NULL){
    puts("Error during the lecture of .dat file");
    exit(0);
  }

  /* Read number of lines with a coment */
  fscanf(Sim_dat,"%d",
	 &num_coment);

  /* Skip commented lines */
  for(int i = 0;i<=num_coment;i++){
    fgets(line, sizeof(line), Sim_dat);
  }

  /* Read :
     - Number of nodes
     - Number of elements  
     - Number of spatial dimensions
     - Kind of elements
  */
  fgets(line, sizeof(line), Sim_dat); /* Comented line */
  fgets(line, sizeof(line), Sim_dat); /* Line with usefull data */
  aux_read = sscanf(line,"%d%d%d%d",
	 &NumberNodes,
	 &NumberElements,
	 &SpatialDimensions,
	 &KindOfElement);
  if(aux_read < 4){
    printf("Insuficient parameters !!! \n");
    printf("Check if you miss some of this parameters : \n");
    printf("\t Number of nodes \n");
    printf("\t Number of elements \n");
    printf("\t Number of spatial dimensions \n");
    printf("\t Kind of element for the mesh \n");    
    exit(0);
  }
  if(aux_read > 4){
    printf("Too many parameters !!! \n");
    printf("Check if you add some aditional parameter : \n");
    printf("\t Number of nodes \n");
    printf("\t Number of elements \n");
    printf("\t Number of spatial dimensions \n");
    printf("\t Kind of element for the mesh \n");    
    exit(0);
  }
  if((SpatialDimensions == 1) || (SpatialDimensions == 3)){
    printf("For the moment I am only able to deal with 2D problems \n");
    exit(0);
  }
  if( (SpatialDimensions != 1) &&
      (SpatialDimensions != 2) &&
      (SpatialDimensions != 3) ){
    printf(" For the moment our Universe it is only \n");
    printf(" able to deal with 1, 2 or 3 spatial dimensions. \n ");
    printf(" Wait for future updtes ... \n");
    exit(0);
  }

  /* Read Number of fields and chek for dummie mistakes of the users */
  fgets(line, sizeof(line), Sim_dat); /* Comented line */
  fgets(line, sizeof(line), Sim_dat); /* Line with usefull data */
  aux_read = sscanf(line,"%d",&NumberFields);
  if(aux_read > 1){
    printf("Too many parameters !!! \n");
    printf("Check if you add some aditional parameter : \n");
    printf("\t Number of fields for the simulation \n");
    exit(0);
  }
  if(aux_read < 1){
    printf("Insuficient parameters !!! \n");
    printf("Check if you miss some of this parameters : \n");
    printf("\t Number of fields for the simulation \n");
    exit(0);
  }
  if(NumberFields > 14){
    printf(" This problem is only able to deal with a maximum \n");
    printf(" of 14 different parameter simoultanesly \n");
    exit(0);
  }
  if(NumberFields <= 0){
    printf(" The minimum number of fields to simulate is 1 \nx");
    exit(0);
  }

  /* Allocate an array with the list of the fields */
  Fields = (double *) AllocateArray(NumberNodes*NumberFields,sizeof(double));

  fgets(line, sizeof(line), Sim_dat); /* Usefull data */
  aux_read = sscanf(line,"%d%d%d%d%d%d%d%d%d%d%d%d%d%f",
	 &AuxiliarArrayFields[0],
	 &AuxiliarArrayFields[1], 
	 &AuxiliarArrayFields[2], 
	 &AuxiliarArrayFields[3], 
	 &AuxiliarArrayFields[4], 
	 &AuxiliarArrayFields[5], 
	 &AuxiliarArrayFields[6], 
	 &AuxiliarArrayFields[7],
	 &AuxiliarArrayFields[8], 
	 &AuxiliarArrayFields[9], 
	 &AuxiliarArrayFields[10], 
	 &AuxiliarArrayFields[11], 
	 &AuxiliarArrayFields[12], 
	 &AuxiliarArrayFields[13]); 

  
  // Populate the array with pointers.
  for (int i=0; i<numPointers; i++) {
    pointerLUT[i] = NewFoo(); // NewFoo() should return (Foo *)
  }  
  // Access the LUT.
  OneField * double = Fields[i];
  
  
  /* Read the list of fields to simulate */
  fgets(line, sizeof(line), Sim_dat); /* Comented line */
  switch(NumberFields){
  case 1  :
    fgets(line, sizeof(line), Sim_dat); /* Usefull data */
    sscanf(line,"%lf",
	   &Fields[0]); 
    break;
    
  case 2  :
    fgets(line, sizeof(line), Sim_dat); /* Usefull data */
    sscanf(line,"%lf%lf",
	   &Fields[i + 0],
	   &Fields[i + 1]); 
    break;
    
  case 3  :
    fgets(line, sizeof(line), Sim_dat); /* Usefull data */
    sscanf(line,"%lf%lf%lf",
	   &Fields[i + 0],
	   &Fields[i + 1],
	   &Fields[i + 2]); 
    break;
    
  case 4  :
    fgets(line, sizeof(line), Sim_dat); /* Usefull data */
    sscanf(line,"%lf%lf%lf%lf",
	   &Fields[i + 0],
	   &Fields[i + 1],
	   &Fields[i + 2],
	   &Fields[i + 3]); 
    break;
    
  case 5  :
    fgets(line, sizeof(line), Sim_dat); /* Usefull data */
    sscanf(line,"%lf%lf%lf%lf%f",
	   &Fields[i + 0],
	   &Fields[i + 1],
	   &Fields[i + 2],
	   &Fields[i + 3],
	   &Fields[i + 4]); 
    break;
  case 6  :
    break;
  case 7  :
    break;
  case 8  :
    break;
  case 9  :
    break;
  case 10 :
    break;
  case 11 :
    break;
  case 12 :
    break;
  case 13 :
    break;
  case 14 :
    break;
  default :
    puts("Too much number of fiels \n");
    exit(0);
  }

  
  for(int i = 0;i<NumberFields;i++){
    fgets(line, sizeof(line), Sim_dat);
    sscanf(line,"%d",&Fields[i]);    
  }


  /* Read delta_t and number of time-steps */
  fgets(line, sizeof(line), Sim_dat);
  fgets(line, sizeof(line), Sim_dat);
  sscanf(line,"%lf%d",
	 &DeltaT,
	 &N_steps);

  /* Type of problem, type of element and type of time integration */
  fgets(line, sizeof(line), Sim_dat);
  fgets(line, sizeof(line), Sim_dat);
  sscanf(line,"%d%d%d",
	 &type_problem,
	 &type_elem,
	 &temp_integration);
  
  /* Allocate data ( Depends of the type of simulation ) */
  coordinates_1D = (double * ) AllocateArray(N_nodes,sizeof(double));
  conectivity_1D = (int * ) AllocateArray(N_elem*2,sizeof(int));
  
  /* Read coordinates of the node of the mesh and assemble it */
  fgets(line, sizeof(line), Sim_dat);
  for(int i = 0;i<N_nodes;i++){
    
    /* Read each nodal position */
    fgets(line, sizeof(line), Sim_dat);
    sscanf(line,"%d %lf",
	   &ix, /* Index of the node */
	   &coordinates_1D[i]); /* Coordinate of the node */
  }


  
  /* Read conectivity of the mesh */
  fgets(line, sizeof(line), Sim_dat);
  for(int i = 0;i<N_elem;i++){

    /* Read each element conectivity */
    fgets(line, sizeof(line), Sim_dat);
    sscanf(line,"%d%d%d",
	     &ix,
	     &conectivity_1D[i*2],
	     &conectivity_1D[i*2+1]);
  }
  
  /* Read physical parameters */
  fgets(line, sizeof(line), Sim_dat);
  fgets(line, sizeof(line), Sim_dat);
  sscanf(line,"%lf%lf",
	 &A_el,
	 &g);

  /* Allocate array of variables in the t = n */
  phi_n = (double * ) AllocateArray(N_nodes,sizeof(double));
  
  /* Read initial conditions and save them in the phi_n array */
  fgets(line, sizeof(line), Sim_dat);
  for(int i = 0;i<N_nodes;i++){
    fgets(line, sizeof(line), Sim_dat);
    sscanf(line,"%lf",
	   &phi_n[i]);
  }
  
  /* /\* Read boundary conditions *\/ */
  /* fgets(line, sizeof(line), Sim_dat); */
  /* if(type_problem == 1){ */

  /*   /\* Boundary conditions for 1D cases *\/ */
  /*   BoundaryConditions  = (int **) AllocateMatrix(2,3,sizeof(int)); */
    
  /*   /\* X = 0 *\/ */
  /*   fgets(line, sizeof(line), Sim_dat); */
  /*   sscanf(line,"%d%d%d", */
  /* 	   &BoundaryConditions[0][0], /\* Node *\/ */
  /* 	   &BoundaryConditions[0][1], /\* Type of BCC *\/ */
  /* 	   &BoundaryConditions[0][2]); /\* Label of BCC *\/ */

  /*   /\* X = L *\/ */
  /*   fgets(line, sizeof(line), Sim_dat); */
  /*   sscanf(line,"%d%d%d", */
  /* 	   &BoundaryConditions[1][0], /\* Node *\/ */
  /* 	   &BoundaryConditions[1][1], /\* Type of BCC *\/ */
  /* 	   &BoundaryConditions[1][2]); /\* Label of BCC *\/ */
    
  /* } /\* if(type_problem == 1) *\/ */
     
    
  // Close .dat file
  free(Name_Simulation); // Free the memory gained from malloc. 
  fclose(Sim_dat);
  

} /* void read_dat(char * Name_File) */



/***************************************************************************/
