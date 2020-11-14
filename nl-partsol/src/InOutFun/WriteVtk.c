#include "nl-partsol.h"

void WriteVtk_MPM(char * Name_File, GaussPoint MPM_Mesh,
		  char * List_Fields, int TimeStep_i, int ResultsTimeStep)
{

  /* Number of dimensions */
  int Ndim = NumberDimensions;

  FILE * Vtk_file;
  char Name_file_t[80];
  double P_GP; /* Trace of the stress tensor (volumetric) */
  Tensor Stress_p, EV_Stress_p; 
  Tensor Strain_p, EV_Strain_p; 

  sprintf(Name_file_t,"%s/MPM_%s_%i.vtk",OutputDir,Name_File,TimeStep_i);
  
  Vtk_file = fopen(Name_file_t,"w");

  /* Header */
  fprintf(Vtk_file,"# vtk DataFile Version 3.0 \n");
  fprintf(Vtk_file,"Results time step %i \n",ResultsTimeStep);
  fprintf(Vtk_file,"ASCII \n");
  fprintf(Vtk_file,"DATASET UNSTRUCTURED_GRID \n");

  /* Coordinates */
  fprintf(Vtk_file,"POINTS %i double \n",MPM_Mesh.NumGP);
  for(int i = 0 ; i<MPM_Mesh.NumGP ; i++){
    for(int j = 0 ; j<3 ; j++){
      if(j<Ndim){
	fprintf(Vtk_file,"%lf ",MPM_Mesh.Phi.x_GC.nM[i][j]);
      }
      else{
	fprintf(Vtk_file,"%lf ",0.0);
      }
    }
    fprintf(Vtk_file,"\n");
  }

  /* Connectivity */
  fprintf(Vtk_file,"CELLS %i %i \n",MPM_Mesh.NumGP,Ndim*MPM_Mesh.NumGP);
  for(int i = 0 ; i<MPM_Mesh.NumGP ; i++){
    fprintf(Vtk_file,"%i %i \n",1,i);
  }

  /* Type of element */
  fprintf(Vtk_file,"CELL_TYPES %i \n",MPM_Mesh.NumGP);
  for(int i = 0 ; i<MPM_Mesh.NumGP ; i++){
    fprintf(Vtk_file,"%i \n",1);
  }

  /* Point data */
  fprintf(Vtk_file,"POINT_DATA %i \n",MPM_Mesh.NumGP);
  fprintf(Vtk_file,"VECTORS %s double \n","X_GC");
  for(int i =  0 ; i<MPM_Mesh.NumGP ; i++){
    for(int j = 0 ; j<3 ; j++){
      if(j<Ndim){
	fprintf(Vtk_file,"%lf ",MPM_Mesh.Phi.x_GC.nM[i][j]);
      }
      else{
	fprintf(Vtk_file,"%lf ",0.0);
      }
    }
    fprintf(Vtk_file,"\n");
  }
  fprintf(Vtk_file,"VECTORS %s double \n","X_EC");
  for(int i =  0 ; i<MPM_Mesh.NumGP ; i++){
    for(int j = 0 ; j<3 ; j++){
      if(j<Ndim){
	fprintf(Vtk_file,"%lf ",MPM_Mesh.Phi.x_EC.nM[i][j]);
      }
      else{
	fprintf(Vtk_file,"%lf ",0.0);
      }
    }
    fprintf(Vtk_file,"\n");
  }

  /* Cell data */  
  fprintf(Vtk_file,"CELL_DATA %i \n",MPM_Mesh.NumGP);

  fprintf(Vtk_file,"SCALARS MASS double \n");
  fprintf(Vtk_file,"LOOKUP_TABLE default \n");
  for(int i =  0 ; i<MPM_Mesh.NumGP ; i++){
    fprintf(Vtk_file,"%lf \n",MPM_Mesh.Phi.mass.nV[i]);
  }

  fprintf(Vtk_file,"SCALARS DENSITY double \n");
  fprintf(Vtk_file,"LOOKUP_TABLE default \n");
  for(int i =  0 ; i<MPM_Mesh.NumGP ; i++){
    fprintf(Vtk_file,"%lf \n",MPM_Mesh.Phi.rho.nV[i]); 
  }

  fprintf(Vtk_file,"SCALARS P double \n");
  fprintf(Vtk_file,"LOOKUP_TABLE default \n");
  for(int i =  0 ; i<MPM_Mesh.NumGP ; i++){
    P_GP = 0;
    for(int j = 0 ; j<3 ; j++ ){
      if(j < Ndim){
	P_GP += MPM_Mesh.Phi.Stress.nM[i][j];
      }
      else{
	P_GP = 0;
      }
    }    
    fprintf(Vtk_file,"%lf \n",P_GP/Ndim);
  }

  fprintf(Vtk_file,"SCALARS W double \n");
  fprintf(Vtk_file,"LOOKUP_TABLE default \n");
  for(int i =  0 ; i<MPM_Mesh.NumGP ; i++){
    fprintf(Vtk_file,"%lf \n",MPM_Mesh.Phi.W.nV[i]); 
  }

  fprintf(Vtk_file,"SCALARS Chi double \n");
  fprintf(Vtk_file,"LOOKUP_TABLE default \n");
  for(int i =  0 ; i<MPM_Mesh.NumGP ; i++){
    fprintf(Vtk_file,"%lf \n",MPM_Mesh.Phi.chi.nV[i]); 
  }

  fprintf(Vtk_file,"SCALARS EPS double \n");
  fprintf(Vtk_file,"LOOKUP_TABLE default \n");
  for(int i =  0 ; i<MPM_Mesh.NumGP ; i++){
    fprintf(Vtk_file,"%lf \n",MPM_Mesh.Phi.EPS.nV[i]); 
  }

  /* double -> integer */
  fprintf(Vtk_file,"SCALARS ELEM_i integer \n");
  fprintf(Vtk_file,"LOOKUP_TABLE default \n");
  for(int i =  0 ; i<MPM_Mesh.NumGP ; i++){
    fprintf(Vtk_file,"%i \n",MPM_Mesh.I0[i]);
  }

  fprintf(Vtk_file,"SCALARS MatIdx integer \n");
  fprintf(Vtk_file,"LOOKUP_TABLE default \n");
  for(int i =  0 ; i<MPM_Mesh.NumGP ; i++){
    fprintf(Vtk_file,"%i \n",MPM_Mesh.MatIdx[i]);
  }

  fprintf(Vtk_file,"VECTORS VELOCITY double \n");
  for(int i =  0 ; i<MPM_Mesh.NumGP ; i++){
    for(int j = 0 ; j<3 ; j++){
      if(j<Ndim){
	fprintf(Vtk_file,"%lf ",MPM_Mesh.Phi.vel.nM[i][j]);
      }
      else{
	fprintf(Vtk_file,"%lf ",0.0);
      }
    }
    fprintf(Vtk_file,"\n");
  }

  fprintf(Vtk_file,"VECTORS ACCELERATION double \n");
  for(int i =  0 ; i<MPM_Mesh.NumGP ; i++){
    for(int j = 0 ; j<3 ; j++){
      if(j<Ndim){
	fprintf(Vtk_file,"%lf ",MPM_Mesh.Phi.acc.nM[i][j]);
      }
      else{
	fprintf(Vtk_file,"%lf ",0.0);
      }
    }
    fprintf(Vtk_file,"\n");
  }

  fprintf(Vtk_file,"VECTORS DISPLACEMENT double \n");
  for(int i =  0 ; i<MPM_Mesh.NumGP ; i++){
    for(int j = 0 ; j<3 ; j++){
      if(j<Ndim){
	fprintf(Vtk_file,"%lf ",MPM_Mesh.Phi.dis.nM[i][j]);
      }
      else{
	fprintf(Vtk_file,"%lf ",0.0);
      }
    }
    fprintf(Vtk_file,"\n");
  }

  fprintf(Vtk_file,"TENSORS STRESS double \n");
  for(int i =  0 ; i<MPM_Mesh.NumGP ; i++){
    for(int j = 0 ; j<3 ; j++){
      for(int k = 0 ; k<3 ; k++){
	if((j<Ndim) && (k<Ndim)){
	  fprintf(Vtk_file,"%lf ",MPM_Mesh.Phi.Stress.nM[i][j*Ndim+k]);
	}
	else{
	  fprintf(Vtk_file,"%lf ",0.0);
	}
      }
      fprintf(Vtk_file,"\n");
    }
    fprintf(Vtk_file,"\n");
  }

  fprintf(Vtk_file,"VECTORS STRESS_EV double \n");
  for(int i =  0 ; i<MPM_Mesh.NumGP ; i++){
    Stress_p = memory_to_tensor__TensorLib__(MPM_Mesh.Phi.Stress.nM[i], 2);
    EV_Stress_p = Eigenvalues__TensorLib__(Stress_p);
    for(int j = 0 ; j<3 ; j++){
      if(j<Ndim){
	fprintf(Vtk_file,"%lf ",EV_Stress_p.n[j]);
      }
      else{
	fprintf(Vtk_file,"%lf ",0.0);
      }
    }
    fprintf(Vtk_file,"\n");
    free__TensorLib__(EV_Stress_p);
  }

  fprintf(Vtk_file,"TENSORS STRAIN double \n");
  for(int i =  0 ; i<MPM_Mesh.NumGP ; i++){
    for(int j = 0 ; j<3 ; j++){
      for(int k = 0 ; k<3 ; k++){
	if((j<Ndim) && (k<Ndim)){
	  fprintf(Vtk_file,"%lf ",MPM_Mesh.Phi.Strain.nM[i][j*Ndim+k]);
	}
	else{
	  fprintf(Vtk_file,"%lf ",0.0);
	}
      }
      fprintf(Vtk_file,"\n");
    }
    fprintf(Vtk_file,"\n");
  }

  fprintf(Vtk_file,"VECTORS STRAIN_EV double \n");
  for(int i =  0 ; i<MPM_Mesh.NumGP ; i++){
    Strain_p = memory_to_tensor__TensorLib__(MPM_Mesh.Phi.Strain.nM[i], 2);
    EV_Strain_p = Eigenvalues__TensorLib__(Strain_p);
    for(int j = 0 ; j<3 ; j++){
      if(j<Ndim){
	fprintf(Vtk_file,"%lf ",EV_Strain_p.n[j]);
      }
      else{
	fprintf(Vtk_file,"%lf ",0.0);
      }
    }
    fprintf(Vtk_file,"\n");
    free__TensorLib__(EV_Strain_p);
  }

  fprintf(Vtk_file,"TENSORS DEFORMATION-GRADIENT double \n");
  for(int i =  0 ; i<MPM_Mesh.NumGP ; i++){
    for(int j = 0 ; j<3 ; j++){
      for(int k = 0 ; k<3 ; k++){
	if((j<Ndim) && (k<Ndim)){
	  fprintf(Vtk_file,"%lf ",MPM_Mesh.Phi.F_n.nM[i][j*Ndim+k]);
	}
	else{
	  fprintf(Vtk_file,"%lf ",0.0);
	}
      }
      fprintf(Vtk_file,"\n");
    }
    fprintf(Vtk_file,"\n");
  }

  /* Close the file */
  fclose(Vtk_file);
  
}

/*********************************************************************/

void WriteVtk_FEM(char * Name_File, Mesh ElementMesh,
		  Matrix List_Nod_Fields, int TimeStep_i){

  /* Number of dimensions */
  int Ndim = NumberDimensions;
  FILE * Vtk_file;
  char Name_file_t[80];
  int NumberFields;
  char * FieldsList[MAXW] = {NULL};
  int i_Field;
  int NumNodesElem;
  ChainPtr Elem_Conn;

  sprintf(Name_file_t,"%s/FEM_%s_%i.vtk",OutputDir,Name_File,TimeStep_i);
  
  Vtk_file = fopen(Name_file_t,"w");

  /* Header */
  fprintf(Vtk_file,"# vtk DataFile Version 3.0 \n");
  fprintf(Vtk_file,"vtk output \n");
  fprintf(Vtk_file,"ASCII \n");
  fprintf(Vtk_file,"DATASET UNSTRUCTURED_GRID \n");

  /* Coordinates */
  fprintf(Vtk_file,"POINTS %i double \n",ElementMesh.NumNodesMesh);
  for(int i = 0 ; i<ElementMesh.NumNodesMesh ; i++){
    for(int j = 0 ; j<3 ; j++){
      if(j<Ndim){
	fprintf(Vtk_file,"%lf ",ElementMesh.Coordinates.nM[i][j]);
      }
      else{
	fprintf(Vtk_file,"%lf ",0.0);
      }
    }
    fprintf(Vtk_file,"\n");
  }


  /* TEMPORARY --> FIX IT*/
  NumNodesElem = ElementMesh.NumNodesElem[0];
  /* Connectivity */
  fprintf(Vtk_file,"CELLS %i %i \n",
	  ElementMesh.NumElemMesh,
	  ElementMesh.NumElemMesh*(1+NumNodesElem));
  for(int i = 0 ; i<ElementMesh.NumElemMesh ; i++){
    fprintf(Vtk_file,"%i",NumNodesElem);
    Elem_Conn = ElementMesh.Connectivity[i];
    while(Elem_Conn != NULL){
      fprintf(Vtk_file," %i ",Elem_Conn->I);
      Elem_Conn = Elem_Conn->next;
    }
    fprintf(Vtk_file,"\n");
  }

  /* Type of element */
  fprintf(Vtk_file,"CELL_TYPES %i \n",ElementMesh.NumElemMesh);
  for(int i = 0 ; i<ElementMesh.NumElemMesh ; i++){
    fprintf(Vtk_file,"%i \n",9);
  }

  /* Point data */
  fprintf(Vtk_file,"POINT_DATA %i \n",ElementMesh.NumNodesMesh);
  fprintf(Vtk_file,"VECTORS %s double \n","X_GC");
  for(int i =  0 ; i<ElementMesh.NumNodesMesh ; i++){
    for(int j = 0 ; j<3 ; j++){
      if(j<Ndim){
	fprintf(Vtk_file,"%lf ",ElementMesh.Coordinates.nM[i][j]);
      }
      else{
	fprintf(Vtk_file,"%lf ",0.0);
      }
    }
    fprintf(Vtk_file,"\n");
  }

  /* Active nodes */
  fprintf(Vtk_file,"SCALARS Active_Nod_Active int \n");
  fprintf(Vtk_file,"LOOKUP_TABLE default \n");
  for(int i = 0 ; i<ElementMesh.NumNodesMesh ; i++){
    fprintf(Vtk_file,"%i\n",
	    ElementMesh.NumParticles[i]);
  } 

  if(List_Nod_Fields.nM != NULL){
    
    NumberFields = parse (FieldsList,List_Nod_Fields.Info,";\n");

    /* Auxiliar index for the input fields */
    i_Field = 0;

    /* Loop over the fields */
    for(int i = 0 ; i<NumberFields ; i++){
      if(strcmp(FieldsList[i],"MASS") == 0){
	fprintf(Vtk_file,"SCALARS Nod_Mass double \n");
	fprintf(Vtk_file,"LOOKUP_TABLE default \n");
	for(int j =  0 ; j<ElementMesh.NumNodesMesh ; j++){
	  fprintf(Vtk_file,"%lf \n",List_Nod_Fields.nM[j][i_Field]);
	}
	/* Update the index of the field */
	i_Field += 1;
      }

      if(strcmp(FieldsList[i],"MOMENTUM") == 0){
	fprintf(Vtk_file,"VECTORS %s double \n","MOMENTUM");
	for(int j =  0 ; j<ElementMesh.NumNodesMesh ; j++){  
	  /* Print the dimensions of the array */
	  for(int k = 0 ; k<3 ; k++){
	    if(k<Ndim){
	      fprintf(Vtk_file,"%lf ",List_Nod_Fields.nM[j][i_Field+k]);
	    }
	    else{
	      fprintf(Vtk_file,"%lf ",0.0);
	    }

	  }
	  fprintf(Vtk_file,"\n");	
	}
	/* Update the index of the field */
	i_Field += Ndim;
      }

      if(strcmp(FieldsList[i],"ACCELERATION_t0") == 0){
	fprintf(Vtk_file,"VECTORS %s double \n","ACCELERATION_t0");
	for(int j =  0 ; j<ElementMesh.NumNodesMesh ; j++){
	  /* Print the dimensions of the array */
	  for(int k = 0 ; k<3 ; k++){
	    if(k<Ndim){
	      fprintf(Vtk_file,"%lf ",List_Nod_Fields.nM[j][i_Field+k]);
	    }
	    else{
	      fprintf(Vtk_file,"%lf ",0.0);
	    }    
	  }
	  fprintf(Vtk_file,"\n");	
	}
	/* Update the index of the field */
	i_Field += Ndim;
      }

      if(strcmp(FieldsList[i],"ACCELERATION_t1") == 0){
	fprintf(Vtk_file,"VECTORS %s double \n","ACCELERATION_t1");
	for(int j =  0 ; j<ElementMesh.NumNodesMesh ; j++){
	  /* Print the dimensions of the array */
	  for(int k = 0 ; k<3 ; k++){
	    if(k<Ndim){
	      fprintf(Vtk_file,"%lf ",List_Nod_Fields.nM[j][i_Field+k]);
	    }
	    else{
	      fprintf(Vtk_file,"%lf ",0.0);
	    }    
	  }
	  fprintf(Vtk_file,"\n");	
	}
	/* Update the index of the field */
	i_Field += Ndim;
      }

      if(strcmp(FieldsList[i],"VELOCITY") == 0){
	fprintf(Vtk_file,"VECTORS %s double \n","VELOCITY");
	for(int j =  0 ; j<ElementMesh.NumNodesMesh ; j++){
	  /* Print the dimensions of the array */
	  for(int k = 0 ; k<3 ; k++){
	    if(k<Ndim){
	      fprintf(Vtk_file,"%lf ",List_Nod_Fields.nM[j][i_Field+k]);
	    }
	    else{
	      fprintf(Vtk_file,"%lf ",0.0);
	    }    
	  }
	  fprintf(Vtk_file,"\n");	
	}
	/* Update the index of the field */
	i_Field += Ndim;
      }


      if(strcmp(FieldsList[i],"F_INT") == 0){
	fprintf(Vtk_file,"VECTORS %s double \n","F_INT");
	for(int j =  0 ; j<ElementMesh.NumNodesMesh ; j++){
	  /* Print the dimensions of the array */
	  for(int k = 0 ; k<3 ; k++){
	    if(k<Ndim){
	      fprintf(Vtk_file,"%lf ",List_Nod_Fields.nM[j][i_Field+k]);
	    }
	    else{
	      fprintf(Vtk_file,"%lf ",0.0);
	    }    
	  }
	  fprintf(Vtk_file,"\n");	
	}
	/* Update the index of the field */
	i_Field += Ndim;
      }


      if(strcmp(FieldsList[i],"F_GRAV") == 0){
	fprintf(Vtk_file,"VECTORS %s double \n","F_GRAV");
	for(int j =  0 ; j<ElementMesh.NumNodesMesh ; j++){
	  /* Print the dimensions of the array */
	  for(int k = 0 ; k<3 ; k++){
	    if(k<Ndim){
	      fprintf(Vtk_file,"%lf ",List_Nod_Fields.nM[j][i_Field+k]);
	    }
	    else{
	      fprintf(Vtk_file,"%lf ",0.0);
	    }    
	  }
	  fprintf(Vtk_file,"\n");	
	}
	/* Update the index of the field */
	i_Field += Ndim;
      }

      if(strcmp(FieldsList[i],"F_TOT") == 0){
	fprintf(Vtk_file,"VECTORS %s double \n","F_TOT");
	for(int j =  0 ; j<ElementMesh.NumNodesMesh ; j++){
	  /* Print the dimensions of the array */
	  for(int k = 0 ; k<3 ; k++){
	    if(k<Ndim){
	      fprintf(Vtk_file,"%lf ",List_Nod_Fields.nM[j][i_Field+k]);
	    }
	    else{
	      fprintf(Vtk_file,"%lf ",0.0);
	    }    
	  }
	  fprintf(Vtk_file,"\n");	
	}
	/* Update the index of the field */
	i_Field += Ndim;
      }       

      
      if(strcmp(FieldsList[i],"REACTIONS") == 0){
	fprintf(Vtk_file,"VECTORS %s double \n","REACTIONS");
	for(int j =  0 ; j<ElementMesh.NumNodesMesh ; j++){
	  /* Print the dimensions of the array */
	  for(int k = 0 ; k<3 ; k++){
	    if(k<Ndim){
	      fprintf(Vtk_file,"%lf ",List_Nod_Fields.nM[j][i_Field+k]);
	    }
	    else{
	      fprintf(Vtk_file,"%lf ",0.0);
	    }    
	  }
	  fprintf(Vtk_file,"\n");	
	}
	/* Update the index of the field */
	i_Field += NumberDimensions ;
      }

    }

  }
  
  /* Cell data */  
  fprintf(Vtk_file,"CELL_DATA %i \n",ElementMesh.NumElemMesh);

  /* Close the file */
  fclose(Vtk_file);

}

/*********************************************************************/

void WriteVtk_Float_Scalar(char * Name_File, Matrix Field){

  FILE * Vtk_file; 
  Vtk_file = fopen(Name_File,"a");

  fprintf(Vtk_file,"SCALARS %s float \n",Field.Info);
  fprintf(Vtk_file,"LOOKUP_TABLE default \n");
  for(int j =  0 ; j<Field.N_cols ; j++){
    fprintf(Vtk_file,"%lf \n",Field.nV[j]);
  }

  /* Close the file */
  fclose(Vtk_file);

}

/*********************************************************************/

void WriteVtk_Float_Vector(char * Name_File, Matrix Field){  

  int Ndim = NumberDimensions;
  FILE * Vtk_file; 
  Vtk_file = fopen(Name_File,"a");
  
  fprintf(Vtk_file,"VECTORS %s float \n",Field.Info);
  for(int i =  0 ; i<Field.N_cols ; i++){
    /* Print the dimensions of the array */
    for(int j = 0 ; j<Ndim ; j++){
      fprintf(Vtk_file,"%lf ",Field.nM[j][i]);
    }
    fprintf(Vtk_file,"\n");	
  }

  /* Close the file */
  fclose(Vtk_file); 
  
}

/*********************************************************************/

void WriteVtk_Float_Tensor(char * Name_File, Matrix Field){  

  int Ndim = NumberDimensions;
  FILE * Vtk_file; 
  Vtk_file = fopen(Name_File,"a");
  
  fprintf(Vtk_file,"TENSORS %s float \n",Field.Info);
  for(int i =  0 ; i<Field.N_cols ; i++){
    for(int j = 0 ; j<Ndim ; j++){
      for(int k = 0 ; k<Ndim ; k++){
	fprintf(Vtk_file,"%lf ",Field.nM[i][j*Ndim+k]);
      }
      fprintf(Vtk_file,"\n");
    }
    fprintf(Vtk_file,"\n");
  }
 
  /* Close the file */
  fclose(Vtk_file); 
  
}

/*********************************************************************/
