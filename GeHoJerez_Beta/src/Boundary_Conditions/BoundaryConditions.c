#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "../ToolsLib/TypeDefinitions.h"
#include "../ToolsLib/GlobalVariables.h"
#include "../ToolsLib/Utils.h"
#include "../InOutFun/InOutFun.h"


/***************************************************************************/

void Read_FEM_BCC(char * Name_File, Mesh FEM_Mesh)
/*
  Read the boundary conditions file :
  Inputs
  - Name_file : Name of the file
  BCC_BOTTOM DIR#[0,1] CURVE#{curve.txt}
  Note : Only read those lines with a BCC in the init 
*/
{

  /* Create a array with the number of steps and fill it with  V#[1]={-1} V#[6]={0} */
  BoundaryConditions BCC;
  
  /* Simulation file */
  FILE * Sim_dat;

  /* Auxiliar variable for reading the lines in the files */
  char line[MAXC] = {0};
  char line_nodes[MAXC] = {0};
  
  /* Number of element in the line , just for check */
  int nkwords,nparam;
  char * kwords[MAXW] = {NULL};
  char * param[MAXW] = {NULL};

  /* Variable for the load parser */
  int aux_DIR;
  char * DIR[MAXW] = {NULL};

  int aux_CURVE;
  char * CURVE[MAXW] = {NULL};

  printf("************************************************* \n");
  printf("Begin of set boundary conditions !!! \n");
  printf(" * Begin of read boundary files : %s \n",Name_File);
  
  /* Open and check .bcc file */
  Sim_dat = fopen(Name_File,"r");  
  if (Sim_dat==NULL){
    puts("Error during the lecture of .bcc file");
    exit(0);
  }

  /* Read the file line by line */
  while( fgets(line, sizeof line, Sim_dat) != NULL ){

    /* Read the line with the space as separators */
    nkwords = parse (kwords, line," \n");

    /* Set BCC in those node in the left of the domain */
    if ( ( strcmp(kwords[0],"BCC_TOP") == 0 ) ||
	 ( strcmp(kwords[0],"BCC_BOTTOM") == 0 ) ||
	 ( strcmp(kwords[0],"BCC_RIGHT") == 0 ) ||
	 ( strcmp(kwords[0],"BCC_LEFT") == 0) ){

        
      for(int i  = 0 ; i<nkwords ; i++){ /* Loop over the words */
	/* Parse the keywords parameters */
	nparam = parse (param,kwords[i],"#\n");
	if(nparam>1){ /* Read only keywords with an asignement */
	  /* Parse the direction of the BCC */
	  if(strcmp(param[0],"DIR") == 0){ 

	    /* Read the DOFS to impose the BCC :*/
	    aux_DIR = parse(DIR,param[1],"[,]\n");
	    if(aux_DIR != NumberDimensions){
	      printf("Error in ReadBCC() : Invalid direction in BCC !!! ");
	      exit(0);
	    }
	    BCC.Dir = (int *)Allocate_ArrayZ(aux_DIR,sizeof(int));
	    for(int j = 0 ; j<aux_DIR ; j++){
	      BCC.Dir[j] = atoi(DIR[j]);
	    }
	    
	  }
	  /* Parse the curve assigned to this BC */
	  if(strcmp(param[0],"CURVE") == 0){
	    /* Read file of the curve */
	    aux_CURVE = parse(CURVE,param[1],"{}\n");
	    if(aux_CURVE != 1){
	      printf("Error in ReadBCC() : Wrong format for the CURVE#{File}");
	      exit(0);
	    }
	    BCC.Value = ReadCurve(CURVE[0]); 
	  }
	  	  
	} 
      } /* Read word by word */

      /* Apply this boundary conditions and copy information of the BCC */
      if(strcmp(kwords[0],"BCC_TOP") == 0){
	BCC.Nodes = FEM_Mesh.TOP;
	BCC.NumNodes = FEM_Mesh.NumTOP;
	strcpy(BCC.Info,"TOP");
      }
      else if(strcmp(kwords[0],"BCC_BOTTOM") == 0){
	BCC.Nodes = FEM_Mesh.BOTTOM;
	BCC.NumNodes = FEM_Mesh.NumBOTTOM;
	strcpy(BCC.Info,"BOTTOM");
      }
      else if(strcmp(kwords[0],"BCC_RIGHT") == 0){
	BCC.Nodes = FEM_Mesh.RIGHT;
	BCC.NumNodes = FEM_Mesh.NumRIGHT;
	strcpy(BCC.Info,"RIGHT");
      }
      else if(strcmp(kwords[0],"BCC_LEFT") == 0){
	BCC.Nodes = FEM_Mesh.LEFT;
	BCC.NumNodes = FEM_Mesh.NumLEFT;
	strcpy(BCC.Info,"LEFT");
      }


      
    } /* End of read BC */
    
  }  /* End of read file */
  printf("End of read boundary conditions file !!! \n");
  fclose(Sim_dat);

  /* Return output */
  return BCC;  
  
} /* BoundayConditions ReadBCC(char * Name_File) */

/**********************************************************************/

/* void Read_MPM_BCC(char * Name_File, GaussPoint GP_Mesh) */
/*   /\* */

/*     - Load format :  */
/*     LOAD_GP DIR#[0,1] CURVE#{curve.txt} NUM_NODES#integer */
/*     LOAD_NODES  */
/*     . */
/*     . integer (NLIST) */
/*     . */
/*     END LIST_NODES */

/*   *\/ */

/* { */


/* /\* Precribed loads in some material point *\/ */

/* if (strcmp(kwords[0],"LOAD_GP") == 0 ){ */

/*   for(int i  = 1 ; i<nkwords ; i++){ */
/* 	nparam = parse (param,kwords[i],"#\n");	 */
/* 	if(nparam == 2){ */

/* 	  if(strcmp(param[0],"DIR") == 0){  */

/* 	    /\* Read the DOFS to impose the BCC :*\/ */
/* 	    aux_DIR = parse(DIR,param[1],"[,]\n"); */
/* 	    if(aux_DIR != NumberDimensions){ */
/* 	      printf("Error in ReadBCC() : Invalid direction in LOAD_GP !!! "); */
/* 	      exit(0); */
/* 	    } */
/* 	    BCC.Dir = (int *)Allocate_ArrayZ(aux_DIR,sizeof(int)); */
/* 	    for(int j = 0 ; j<aux_DIR ; j++){ */
/* 	      BCC.Dir[j] = atoi(DIR[j]); */
/* 	    } */
	    
/* 	  }	   */
/* 	  /\* Parse the curve assigned to this BC *\/ */
/* 	  if(strcmp(param[0],"CURVE") == 0){ */
/* 	    /\* Read file of the curve *\/ */
/* 	    aux_CURVE = parse(CURVE,param[1],"{}\n"); */
/* 	    if(aux_CURVE != 1){ */
/* 	      printf("Error in ReadBCC() : Wrong format for the CURVE#{File}"); */
/* 	      exit(0); */
/* 	    } */
/* 	    BCC.Value = ReadCurve(CURVE[0]);  */
/* 	  } */
/* 	  if(strcmp(param[0],"NUM_NODES") == 0){ // integer */
/* 	    BCC.NumNodes = atoi(param[1]); */
/* 	    BCC.Nodes = */
/* 	      (int *)Allocate_Array(BCC.NumNodes,sizeof(int)); */
/* 	  } */
	  
/* 	}	 */
/*   } */
/* } */
/* if ( strcmp(kwords[0],"LIST_NODES") == 0 ){ */

/*   /\* Fill the list of nodes *\/ */
/*   for(int i = 0 ; i<BCC.NumNodes ; i++){ */
/* 	fgets(line_nodes, sizeof line_nodes, Sim_dat); */
/* 	nparam = parse(param,line_nodes," \n"); */
/* 	if(nparam == 1){ */
/* 	  BCC.Nodes[i] = atoi(param[0]); */
/* 	} */
/* 	else{ */
/* 	  puts("Error in ReadLoads_GP() : Check the list of nodes "); */
/* 	  exit(0); */
/* 	} */
/*   } */
/* } */

/* } */

/* void BCC_GP_Forces(GaussPoint GP_Mesh, BoundaryConditions Loads, int TimeStep) */
/* /\* */
/*   Forces defined in the Gauss Points : */
/*   Inputs */
/* *\/ */
/* { */
/*   /\* 0º  Check the time step *\/ */
/*   if( (TimeStep < 0) || */
/*       (TimeStep > Loads.Value.Num)){ */
/*     puts("Error in BCC_GP_Forces() : The time step is out of the curve !!"); */
/*     exit(0); */
/*   } */
  
/*   /\* 1º Fill the matrix with the local forces *\/ */
/*   for(int i = 0 ; i<Loads.NumNodes ; i++){ */
/*     /\* 2º Check if this GP has a force applied *\/ */
/*     if( (Loads.Nodes[i] > GP_Mesh.NumGP) || */
/* 	(Loads.Nodes[i] < 0)){ */
/*       puts("Error in BCC_GP_Forces() : This GP does not exist !!"); */
/*       exit(0); */
/*     } */
/*     /\* 3º Loop over the dimensions *\/ */
/*     for(int j = 0 ; j<NumberDimensions ; j++){ */
/*       /\* 5º Assign the BC to the nodes *\/ */
/*       if((Loads.Dir[j] == 1) || */
/* 	 (Loads.Dir[j] == -1)){ */
/* 	/\* 6º Apply the force in the node *\/ */
/* 	GP_Mesh.Phi.F.nM[Loads.Nodes[i]][j] = */
/* 	  (double)Loads.Dir[j]*Loads.Value.Fx[TimeStep]; */
/*       }       */
/*     }    */
/*   } */

/* } */

/**********************************************************************/

BoundaryConditions SetBCC(int * Dir, int NumNodes, int * Nodes,
			  char * Info, char * Curve_File){
  
  /* Define boundary conditions */
  BoundaryConditions BCC;

  /* Apply this boundary conditions */
  BCC.Dir = Dir;
  BCC.Nodes = Nodes;
  BCC.NumNodes = NumNodes;

  /* Asign curve of values */
  BCC.Value = ReadCurve(Curve_File);

  /* Copy information of the BCC */
  strcpy(BCC.Info,Info);

  /* Return boundary condition */
  return BCC;
}


/**********************************************************************/

void BCC_Nod_Momentum(Mesh FEM_Mesh, Matrix Nodal_MOMENTUM, int TimeStep)
/*
  Apply the momentum boundary conditions over the nodes 
*/
{

  /* 1º Define auxilar variables */
  int Dim_BCC; /* Dimension where is applied the boundary condition */
  int Id_BCC; /* Index of the node where we apply the BCC */
  

  /* 2º Check the time range */
  if( (TimeStep < 0) ||
      (TimeStep > FEM_Mesh.TOP.Value.Num)){
    puts("Error in BCC_Nodal_MOMENTUM() : The time step is out of the curve !!");
    exit(0);
  }
  /* 2º Loop over the nodes of the TOP boundary */
  for(int i = 0 ; i<FEM_Mesh.TOP.NumNodes ; i++){
    /* 3º Get the index of the element */
    Id_BCC = FEM_Mesh.TOP.Nodes[i];
    /* 4º Loop over the dimensions */
    for(int j = 0 ; j<NumberDimensions ; j++){
      /* 5º Assign the BC to the nodes */
      if((FEM_Mesh.TOP.Dir[j] == 1) ||
	 (FEM_Mesh.TOP.Dir[j] == -1)){
	Nodal_MOMENTUM.nM[j][Id_BCC] =
	  (double)FEM_Mesh.TOP.Dir[j]*FEM_Mesh.TOP.Value.Fx[TimeStep];
      }
    }          
  }

  /* 2º Check the time range */
  if( (TimeStep < 0) ||
      (TimeStep > FEM_Mesh.BOTTOM.Value.Num)){
    puts("Error in BCC_Nodal_MOMENTUM() : The time step is out of the curve !!");
    exit(0);
  }
  /* 2º Loop over the nodes of the BOTTOM boundary */
  for(int i = 0 ; i<FEM_Mesh.BOTTOM.NumNodes ; i++){
    /* 3º Get the index of the element */
    Id_BCC = FEM_Mesh.BOTTOM.Nodes[i];
    /* 4º Loop over the dimensions */
    for(int j = 0 ; j<NumberDimensions ; j++){
      /* 5º Assign the BC to the nodes */
      if((FEM_Mesh.BOTTOM.Dir[j] == 1) ||
	 (FEM_Mesh.BOTTOM.Dir[j] == -1)){
	Nodal_MOMENTUM.nM[j][Id_BCC] =
	  (double)FEM_Mesh.BOTTOM.Dir[j]*FEM_Mesh.BOTTOM.Value.Fx[TimeStep];
      }
    }          
  }

  /* 2º Check the time range */
  if( (TimeStep < 0) ||
      (TimeStep > FEM_Mesh.RIGHT.Value.Num)){
    puts("Error in BCC_Nodal_MOMENTUM() : The time step is out of the curve !!");
    exit(0);
  }
  /* 2º Loop over the nodes of the RIGHT boundary */
  for(int i = 0 ; i<BCC_Momentum.NumNodes ; i++){
    /* 3º Get the index of the element */
    Id_BCC = BCC_Momentum.Nodes[i];
    /* 4º Loop over the dimensions */
    for(int j = 0 ; j<NumberDimensions ; j++){
      /* 5º Assign the BC to the nodes */
      if((BCC_Momentum.Dir[j] == 1) ||
	 (BCC_Momentum.Dir[j] == -1)){
	Nodal_MOMENTUM.nM[j][Id_BCC] =
	  (double)BCC_Momentum.Dir[j]*BCC_Momentum.Value.Fx[TimeStep];
      }
    }          
  }

  
  /* 2º Check the time range */
  if( (TimeStep < 0) ||
      (TimeStep > FEM_Mesh.LEFT.Value.Num)){
    puts("Error in BCC_Nodal_MOMENTUM() : The time step is out of the curve !!");
    exit(0);
  }  
  /* 2º Loop over the nodes of the LEFT boundary */
  for(int i = 0 ; i<BCC_Momentum.NumNodes ; i++){
    /* 3º Get the index of the element */
    Id_BCC = BCC_Momentum.Nodes[i];
    /* 4º Loop over the dimensions */
    for(int j = 0 ; j<NumberDimensions ; j++){
      /* 5º Assign the BC to the nodes */
      if((BCC_Momentum.Dir[j] == 1) ||
	 (BCC_Momentum.Dir[j] == -1)){
	Nodal_MOMENTUM.nM[j][Id_BCC] =
	  (double)BCC_Momentum.Dir[j]*BCC_Momentum.Value.Fx[TimeStep];
      }
    }          
  } 
  

}

/*********************************************************************/



/*********************************************************************/

/* void BC_Nod_Springs(GaussPoint GP_Mesh, BoundaryConditions Loads, int TimeStep) */
/* /\* */
/*   Add to a surface springs to simulate the winkler */
/* *\/ */
/* { */
  
/* } */
