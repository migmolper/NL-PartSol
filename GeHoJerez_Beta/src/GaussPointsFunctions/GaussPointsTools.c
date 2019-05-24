#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "../ToolsLib/TypeDefinitions.h"
#include "../ToolsLib/Utils.h"
#include "../InOutFun/InOutFun.h"
#include "../ElementsFunctions/ShapeFunctions.h"
#include "../ElementsFunctions/ElementTools.h"
#include "../Constitutive/Constitutive.h"
#include "GaussPointsTools.h"

/*********************************************************************/

GaussPoint Initialize_GP_Mesh(char * MPM_GID_MeshName,
			      Matrix InputFields,
			      double Density0,
			      Matrix D)
/*
  
*/
{
  /* Material point mesh (Gauss-Points) */
  Mesh MPM_GID_Mesh;
  GaussPoint MPM_Mesh;
  int Init_Num_GP_Elem;
  int Size_MPM_Mesh;
  int NumFields;
  char * Field[MAXW] = {NULL};
  /* Initialize parser to read files */
  ParserDictionary Dict = InitParserDictionary();
  char * delims = Dict.sep[4];

  /* Screen message */
  printf("************************************************* \n");
  printf("Begin of initialize the Gauss-Points mesh !!! \n");

  /* Read GP mesh */
  MPM_GID_Mesh = ReadGidMesh(MPM_GID_MeshName);

  /* The number of Gauss-Points is the same as the number of elements
   in the input mesh, because we set a GP in the middle of each element */
  MPM_Mesh.NumGP = MPM_GID_Mesh.NumElemMesh;

  /* Allocate fields */
  
  /* Index of the Element */
  MPM_Mesh.Element_id = (int *)Allocate_ArrayZ(MPM_Mesh.NumGP,sizeof(int));

  /* Coordinates of the GP (Global/Local)*/
  MPM_Mesh.Phi.x_GC = MatAllocZ(MPM_Mesh.NumGP,3);
  MPM_Mesh.Phi.x_EC = MatAllocZ(MPM_Mesh.NumGP,2);

  /* Displacement field (Vectorial) */
  MPM_Mesh.Phi.dis = MatAllocZ(MPM_Mesh.NumGP,2);

  /* Velocity field (Vectorial) */
  MPM_Mesh.Phi.vel = MatAllocZ(MPM_Mesh.NumGP,2);

  /* Acceleration field (Vectorial) */
  MPM_Mesh.Phi.acc = MatAllocZ(MPM_Mesh.NumGP,2);
  
  /* Stress field (Tensor) */
  MPM_Mesh.Phi.Stress = MatAllocZ(MPM_Mesh.NumGP,3);

  /* Strain field (Tensor) */
  MPM_Mesh.Phi.Strain = MatAllocZ(MPM_Mesh.NumGP,3);

  /* Mass field (Scalar) */
  MPM_Mesh.Phi.mass = MatAllocZ(MPM_Mesh.NumGP,1);

  /* Density field (Scalar) */
  MPM_Mesh.Phi.rho = MatAllocZ(MPM_Mesh.NumGP,1);

  if(strcmp(MPM_GID_Mesh.TypeElem,"Triangle") == 0 ){
    Matrix a = MatAllocZ(3,1);
    Matrix b = MatAllocZ(3,1);
    Matrix c;
    Matrix A_el;
    int NOD_e[3];

    for(int i = 0 ; i<MPM_Mesh.NumGP ; i++){
      
      /* Element connectivity */
      NOD_e[0] = MPM_GID_Mesh.Connectivity[i][0];
      NOD_e[1] = MPM_GID_Mesh.Connectivity[i][1];
      NOD_e[2] = MPM_GID_Mesh.Connectivity[i][2];
      /* a array */
      a.nV[0] = MPM_GID_Mesh.Coordinates.nM[NOD_e[0]][0] -
	MPM_GID_Mesh.Coordinates.nM[NOD_e[1]][0];
      a.nV[1] = MPM_GID_Mesh.Coordinates.nM[NOD_e[0]][1] -
	MPM_GID_Mesh.Coordinates.nM[NOD_e[1]][1];
      a.nV[2] = 0.0;
      /* b array */
      b.nV[0] = MPM_GID_Mesh.Coordinates.nM[NOD_e[0]][0] -
	MPM_GID_Mesh.Coordinates.nM[NOD_e[2]][0];
      b.nV[1] = MPM_GID_Mesh.Coordinates.nM[NOD_e[0]][1] -
	MPM_GID_Mesh.Coordinates.nM[NOD_e[2]][1];
      b.nV[2] = 0.0;
      /* c array */
      c = Vectorial_prod(a,b);
      /* Get the area of the element */
      A_el = Norm_Mat(c,2);
      A_el.n *= 0.5;
      free(c.nV);
      
      /* Assign the mass parameter */
      MPM_Mesh.Phi.mass.nV[i] = fabs(A_el.n)*Density0;
      /* Set the initial density */
      MPM_Mesh.Phi.rho.nV[i] = Density0;
      /* Get the coordinates of the centre */
      MPM_Mesh.Phi.x_GC.nM[i][0] = (double)1/3*(MPM_GID_Mesh.Coordinates.nM[NOD_e[0]][0] +
						MPM_GID_Mesh.Coordinates.nM[NOD_e[1]][0] +
						MPM_GID_Mesh.Coordinates.nM[NOD_e[2]][0]);
      MPM_Mesh.Phi.x_GC.nM[i][1] = (double)1/3*(MPM_GID_Mesh.Coordinates.nM[NOD_e[0]][1] +
						MPM_GID_Mesh.Coordinates.nM[NOD_e[1]][1] +
						MPM_GID_Mesh.Coordinates.nM[NOD_e[2]][1]);
      MPM_Mesh.Phi.x_GC.nM[i][2] = 0.0;
      /* Local coordinates of the element */
      MPM_Mesh.Element_id[i] = -999;
      /* Location in the natural coordinates
	 of the element (Init to zero) */
      MPM_Mesh.Phi.x_EC.nM[i][0] = 0.0;
      MPM_Mesh.Phi.x_EC.nM[i][1] = 0.0;
      /* Initial displacement */
      MPM_Mesh.Phi.dis.nM[i][0] = 0.0;
      MPM_Mesh.Phi.dis.nM[i][1] = 0.0;
      /* Initial Velocity */
      MPM_Mesh.Phi.vel.nM[i][0] = 0.0;
      MPM_Mesh.Phi.vel.nM[i][1] = 0.0;
      /* Initial Acceleration */
      MPM_Mesh.Phi.acc.nM[i][0] = 0.0;
      MPM_Mesh.Phi.acc.nM[i][1] = 0.0;
      /* Initial Stress */
      MPM_Mesh.Phi.Stress.nM[i][0] = 0.0;
      MPM_Mesh.Phi.Stress.nM[i][1] = 0.0;
      MPM_Mesh.Phi.Stress.nM[i][2] = 0.0;
      /* Initial Strain */
      MPM_Mesh.Phi.Strain.nM[i][0] = 0.0;
      MPM_Mesh.Phi.Strain.nM[i][1] = 0.0;
      MPM_Mesh.Phi.Strain.nM[i][2] = 0.0;      

      
    }

    /* Free auxiliar arrays */
    free(a.nV);
    free(b.nV);

    /* Initialize all the fields :
       - Mass of the material point
       - Position field (Vectorial) in global coordiantes and in element coordinates :
       - Displacement, Velocity and acceleration field (Vectorial)
       - Stress and Strain fields (Tensorial)
       Note : It is not necessary to allocate memory...think about it ;)
    */
    if(InputFields.nM != NULL){ /* Only if we want to do a hot-start */

      /* NumFields = parse (Field, InputFields.Info, ";\n"); */
    
      /* for(int i = 0; i<NumFields ; i++){ */
      /* 	if(strcmp(Field[i],"X_GP")==0) */
      /* 	  MPM_Mesh.Phi.x_GC.nM[0] = InputFields.nM[i]; */

      /* 	if(strcmp(Field[i],"Y_GP")==0) */
      /* 	  MPM_Mesh.Phi.x_GC.nM[1] = InputFields.nM[i]; */

      /* 	if(strcmp(Field[i],"V_X")==0) */
      /* 	  MPM_Mesh.Phi.vel.nM[0] = InputFields.nM[i]; */

      /* 	if(strcmp(Field[i],"V_Y")==0) */
      /* 	  MPM_Mesh.Phi.vel.nM[1] = InputFields.nM[i]; */

      /* 	if(strcmp(Field[i],"SIGMA_X")==0) */
      /* 	  MPM_Mesh.Phi.Stress.nM[0] = InputFields.nM[i]; */

      /* 	if(strcmp(Field[i],"SIGMA_Y")==0) */
      /* 	  MPM_Mesh.Phi.Stress.nM[1] = InputFields.nM[i]; */
      
      /* 	if(strcmp(Field[i],"TAUB_XY")==0) */
      /* 	  MPM_Mesh.Phi.Stress.nM[2] = InputFields.nM[i]; */

      /* 	if(strcmp(Field[i],"MASS")==0) */
      /* 	  MPM_Mesh.Phi.mass.nV = InputFields.nM[i]; */
      /* } */

    }

    /* Allocate and Initialize the constitutive response */
    MPM_Mesh.D = D;
  }


  /* Free the input data */
  free(MPM_GID_Mesh.Coordinates.nM);
  free(MPM_GID_Mesh.Connectivity);
  free(MPM_GID_Mesh.ActiveElem);

  /* Final messages */
  printf("End of initialize the Gauss-Points mesh !!! \n");
  
  return MPM_Mesh;
}

/*********************************************************************/

void LocateGP(GaussPoint MPM_Mesh, Mesh FEM_Mesh, int TimeStep){

  Matrix X_GP;
  X_GP.N_cols = 3;
  X_GP.N_rows = 1;
  X_GP.n = NAN;
  Matrix Poligon = MatAllocZ(FEM_Mesh.NumNodesElem,3);
  int * Poligon_Connectivity;
  
  for(int i = 0 ; i<MPM_Mesh.NumGP ; i++){

    X_GP.nV = MPM_Mesh.Phi.x_GC.nM[i];

    for(int j = 0 ; j<FEM_Mesh.NumElemMesh ; j++){

      /* Connectivity of the Poligon */
      Poligon_Connectivity = FEM_Mesh.Connectivity[j];

      /* Fill the poligon Matrix */
      for(int k = 0 ; k<FEM_Mesh.NumNodesElem ; k++){
	Poligon.nM[k][0] = FEM_Mesh.Coordinates.nM[Poligon_Connectivity[k]][0];
	Poligon.nM[k][1] = FEM_Mesh.Coordinates.nM[Poligon_Connectivity[k]][1];
	Poligon.nM[k][2] = FEM_Mesh.Coordinates.nM[Poligon_Connectivity[k]][2];
      }

      /* Check out if the GP is in the Element */
      if(InOut_Poligon(X_GP,Poligon) == 1){
	/* If the GP is in the element, set the index of the position */
	MPM_Mesh.Element_id[i] = j;
	/* If the GP is in the element, get its natural coordinates */
	MPM_Mesh.Phi.x_EC.nM[i] = GetNaturalCoordinates(Poligon,
							MPM_Mesh.Phi.x_EC.nM[i],
							FEM_Mesh.dNdX_ref);
	break;
      }
      
    } /* Loop over the elements */

  } /* Loop over the GP */

  /* Free memory */
  free(Poligon.nM);
  
}

/*********************************************************************/

Matrix GetMassMatrix_L(Mesh ElementMesh,
		       GaussPoint GP_Mesh){
  
  Matrix M_l = MatAllocZ(ElementMesh.NumNodesMesh,1);
  int Elem_GP_i; /* Element where the GP is placed */
  double M_GP; /* Mass of the Gauss-Point */
  int * Nodes_Elem_GP_i; /* Nodes of the element where the GP is placed */
  Matrix GP_i_XE; /* Element coordinates of the Gauss Points */
  Matrix N_ref_XG; /* Element function evaluated in the Gauss-Points */

  for(int i=0 ; i<GP_Mesh.NumGP ; i++){
     
    /* 4º Index of the element where the G-P is located */
    Elem_GP_i = GP_Mesh.Element_id[i];
    
    /* 5º List of nodes of the element */
    Nodes_Elem_GP_i = ElementMesh.Connectivity[Elem_GP_i];
    
    /* 6º Element coordinates of the G-P */
    GP_i_XE.n = GP_Mesh.Phi.x_EC.nV[i];
    
    /* 7º Evaluate the shape functions of the element in the 
       coordinates of the G-P */
    N_ref_XG = ElementMesh.N_ref(GP_i_XE);

    /* 8º Get the mass of the material point */
    M_GP = GP_Mesh.Phi.mass.nV[i];
    
    /* 8º Iterate over the nodes of the element */
    for(int j=0 ; j<ElementMesh.NumNodesElem ; j++){
      M_l.nV[Nodes_Elem_GP_i[j]-1] += M_GP*N_ref_XG.nV[j];
    }
    
  }

  return M_l;
}


/*********************************************************************/

void GaussPointsToMesh(Mesh ElementMesh,
		       GaussPoint GP_Mesh,
		       Matrix Phi_n_GP,
		       Matrix Phi_n_Nod,
		       Matrix M_l){

  printf(" * Transfer the Gauss-Points information to the mesh \n");
  /************* Transfer the GP information to the mesh ****************/
  /**********************************************************************/
  /*  [N0]-------(e0)-------[N1]-------(e1)-------[N2]-- ;  t = n       */
  /* Phi_n_GP0 <-- Phi_e0 --> Phi_n_GP1 <-- Phi_e1 --> Phi_n_GP         */
  /**********************************************************************/
  
  /* 0º Define auxiliar variable to the transference from
     the nodes to the Gauss Points */
  int NumDOF; /* Degree of freedom of the GP */
  int Elem_GP_i; /* Element where the GP is placed */
  double M_GP; /* Mass of the Gauss-Point */
  int * Nodes_Elem_GP_i; /* Nodes of the element where the GP is placed */
  Matrix GP_i_XE; /* Element coordinates of the Gauss Points */
  Matrix N_ref_XG; /* Element function evaluated in the Gauss-Points */
   
  /* 1º Set the number of degree of freedom of the problem */
  NumDOF = Phi_n_GP.N_rows;

  /* 2º Check if it is necessary to do a resize of Phi_n_Nod  */
  if(ElementMesh.NumNodesMesh != Phi_n_Nod.N_cols){
    if(NumDOF >= 2){ /* 2aº Vectorial */
      free(Phi_n_Nod.nM); /* Free previous data field */
      Phi_n_Nod = MatAllocZ(NumDOF,ElementMesh.NumNodesMesh); /* Resize */
    }
    else{/* 2bº Scalar */
      free(Phi_n_Nod.nV); /* Free previous data field */
      Phi_n_Nod = MatAllocZ(NumDOF,ElementMesh.NumNodesMesh); /* Resize */
    }
  }
  else{ /* 3º If is not necessary to resize, initialize it to zero */
    for(int i = 0 ; i<ElementMesh.NumNodesMesh ; i++){
      if(NumDOF >= 2){ /* 3aº Vectorial */
	for(int j = 0 ; j<NumDOF ; j++){
	  Phi_n_Nod.nM[j][i] = 0; 
	}
      }
      else{ /* 3bº Scalar */
	Phi_n_Nod.nV[i] = 0; 
      }      
    }
  }
  
  for(int i=0 ; i<GP_Mesh.NumGP ; i++){
     
    /* 4º Index of the element where the G-P is located */
    Elem_GP_i = GP_Mesh.Element_id[i];
    
    /* 5º List of nodes of the element */
    Nodes_Elem_GP_i = ElementMesh.Connectivity[Elem_GP_i];
    
    /* 6º Element coordinates of the G-P */
    GP_i_XE.n = GP_Mesh.Phi.x_EC.nV[i];
    
    /* 7º Evaluate the shape functions of the element in the 
       coordinates of the G-P */
    N_ref_XG = ElementMesh.N_ref(GP_i_XE);

    /* 8º Get the mass of the material point */
    M_GP = GP_Mesh.Phi.mass.nV[i];
    
    /* 8º Iterate over the nodes of the element */
    for(int j=0 ; j<ElementMesh.NumNodesElem ; j++){
      /* 9º Iterate over the fields */
      for(int k=0 ; k<NumDOF ; k++){
	/* 10º Update the field value */;
	Phi_n_Nod.nM[k][Nodes_Elem_GP_i[j]-1] += Phi_n_GP.nM[k][i]*M_GP*N_ref_XG.nV[j];
      }
    }
    
  }

  /* 11º Fix the nodal value with the lumped mass matrix */
  for(int i = 0 ; i<ElementMesh.NumNodesMesh ; i++){
    for(int j = 0 ; j<NumDOF ; j++){
      Phi_n_Nod.nM[j][i] /= M_l.nV[i]; 
    }
  }


}

/*********************************************************************/

void MeshToGaussPoints(Mesh ElementMesh,
		       GaussPoint GP_Mesh,
		       Matrix Phi_n_Nod,
		       Matrix Phi_n_GP,
		       Matrix M_l){

  printf(" * Transfer the nodal values to the Gauss-Points \n");
  /************ Transfer the nodal values of U_n1 to the G-P *************/
  /***********************************************************************/
  /*   Phi_n ---> Phi_e <--- Phi_n ---> Phi_e <--- Phi_n                 */
  /*   [N0]-------(e0)-------[N1]-------(e1)-------[N2]-- ;  t = n + 1   */
  /***********************************************************************/

  /* 0º Define auxiliar variable to the transference from
     the Gauss Points to the nodes */
  int NumDOF; /* Degree of freedom of the GP */
  int Elem_GP_i; /* Element where the GP is placed */
  int * Nodes_Elem_GP_i; /* Nodes of the element where the GP is placed */
  Matrix GP_i_XE; /* Element coordinates of the Gauss Points */
  Matrix N_ref_XG; /* Element function evaluated in the Gauss-Points */

  NumDOF = Phi_n_Nod.N_rows;
  
  /* 1º Iterate over the G-P */
  for(int i=0 ; i<GP_Mesh.NumGP ; i++){
    
    /* 2º Index of the element where the G-P is located */
    Elem_GP_i = GP_Mesh.Element_id[i];
    
    /* 3º List of nodes of the element */
    Nodes_Elem_GP_i = ElementMesh.Connectivity[Elem_GP_i];
    
    /* 4º Element coordinates of the G-P */
    GP_i_XE.n = GP_Mesh.Phi.x_EC.nV[i];

    /* 5º Evaluate the shape functions of the element in the coordinates of the G-P */
    N_ref_XG = ElementMesh.N_ref(GP_i_XE);

    /* 6º Set to zero the field value */
    for(int k=0 ; k<NumDOF ; k++){
      /* 7º Update the field value */
      Phi_n_GP.nM[k][i] = 0;
    }    

    /* 8º Iterate over the nodes of the element */
    for(int j=0 ; j<ElementMesh.NumNodesElem ; j++){
      /* 9º Iterate over the fields */
      for(int k=0 ; k<NumDOF ; k++){
  	/* 10º Update the field value */
  	Phi_n_GP.nM[k][i] += N_ref_XG.nV[j]*Phi_n_Nod.nM[k][Nodes_Elem_GP_i[j]-1];
      }
    }
  }
}

/*********************************************************************/

/* double GetJacobian(){ */
/* } */

/*********************************************************************/


/* Matrix Get_Lagrangian_CG_Tensor(GaussPoint * GP) */
/* /\*  */
/*    The Right Cauchy Green Deformation Tensor : */
/*    C = F^{T} \cdot F  */

/*    Inputs : Gauss point */
/*    Output : C */
/* *\/ */
/* { */
/*   /\* Check if we have a null matrix *\/ */
/*   if (GP->F.nM == NULL){ */
/*     puts("Error in Get_Lagrangian_CG_Tensor : GP->F tensor = NULL"); */
/*     exit(0); */
/*   } */

/*   /\* Output, we dont need to allocate it because the memory reserve is done  */
/*      in the function Scalar_prod() *\/ */
/*   Matrix C; */
/*   /\* Make a temporal copy of F because Scalar_prod() destroy the input *\/ */
/*   Matrix F = CopyMat(GP->F); */
/*   /\* Get C *\/ */
/*   C = Scalar_prod(Transpose_Mat(F),F); */

/*   return C; */
  
/* } */

/*********************************************************************/

/* Matrix Get_Eulerian_CG_Tensor(GaussPoint * GP)  */
/* /\* */
/*  The Left Cauchy Green Deformation Tensor : */
/*  B = F \cdot F^{T}  */

/*  Inputs : Gauss point  */
/*  Output : B */
/* *\/ */
/* {   */
/*   /\* Check if we have a null matrix *\/ */
/*   if (GP->F.nM == NULL){ */
/*     puts("Error in Get_Eulerian_CG_Tensor : GP->F tensor = NULL"); */
/*     exit(0); */
/*   } */

/*   /\* Output, we dont need to allocate it because the memory reserve is done  */
/*    in the function Scalar_prod() *\/ */
/*   Matrix B; */
/*   /\* Make a temporal copy of F because Scalar_prod() destroy the input *\/ */
/*   Matrix F = CopyMat(GP->F); */
  
/*   /\* Get B *\/ */
/*   B = Scalar_prod(F,Transpose_Mat(F)); */

/*   return B; */
  
/* } */
