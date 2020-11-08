#include "nl-partsol.h"


/*
  Call global variables
*/
bool Out_global_coordinates;
bool Out_element_coordinates;
bool Out_mass;
bool Out_density;
bool Out_damage;
bool Out_nodal_idx;
bool Out_material_idx;
bool Out_velocity;
bool Out_acceleration;
bool Out_displacement;
bool Out_stress;
bool Out_eigenvalues_stress;
bool Out_volumetric_stress;
bool Out_strain;
bool Out_eigenvalues_strain;
bool Out_deformation_gradient;
bool Out_energy;

/*
  Auxiliar functions 
*/
static void vtk_Out_X_GC(FILE *, Matrix, int);
static void vtk_Out_X_EC(FILE *, Matrix, int);
static void vtk_Out_mass(FILE *, Matrix, int);
static void vtk_Out_density(FILE *, Matrix, int);
static void vtk_Out_damage(FILE *, Matrix, int);
static void vtk_Out_nodal_idx(FILE *, int *, int);
static void vtk_Out_material_idx(FILE *, int *, int);
static void vtk_Out_vel(FILE *, Matrix, int);
static void vtk_Out_acc(FILE *, Matrix, int);
static void vtk_Out_dis(FILE *, Matrix, int);
static void vtk_Out_Stress(FILE *, Matrix, int);
static void vtk_Out_Stress_EV(FILE *, Matrix, int);
static void vtk_Out_Stress_P(FILE *, Matrix, int);
static void vtk_Out_Stain(FILE *, Matrix, int);
static void vtk_Out_Strain_EV(FILE *, Matrix, int);
static void vtk_Out_Deformation_Gradient(FILE *, Matrix, int);
static void vtk_Out_Deformation_Energy(FILE *, Matrix, int);
static void vtk_Out_Kinetic_Energy(FILE *, Matrix, Matrix, int);

/*****************************************************************/

void particle_results_vtk__InOutFun__(char * Name_File, GaussPoint MPM_Mesh, char * List_Fields, int TimeStep_i, int ResultsTimeStep)
{

  /* Number of dimensions */
  int Ndim = NumberDimensions;

  int NumParticles = MPM_Mesh.NumGP;

  FILE * Vtk_file;
  char Name_file_t[80];

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

  /* Print particle global coordiantes */
  if(Out_global_coordinates)
  {
    vtk_Out_X_GC(Vtk_file, MPM_Mesh.Phi.x_GC, NumParticles);
  }

  /* Print particle local coordinates */
  if(Out_element_coordinates)
  {
    vtk_Out_X_EC(Vtk_file, MPM_Mesh.Phi.x_EC, NumParticles);
  }
  
  /* Cell data */  
  fprintf(Vtk_file,"CELL_DATA %i \n",MPM_Mesh.NumGP);

  /* Print particle mass */
  if(Out_mass)
  {
    vtk_Out_mass(Vtk_file, MPM_Mesh.Phi.mass, NumParticles);
  }
  
  /* Print particle density */
  if(Out_density)
  {
    vtk_Out_density(Vtk_file, MPM_Mesh.Phi.rho, NumParticles);
  }

  /* Print particle damage */
  if(Out_damage)
  {
    vtk_Out_damage(Vtk_file, MPM_Mesh.Phi.chi, NumParticles);
  }

  /* Print particle nodal index */
  if(Out_nodal_idx)
  {
    vtk_Out_nodal_idx(Vtk_file, MPM_Mesh.I0, NumParticles);
  }

  /* Print particle material index */
  if(Out_material_idx)
  {
    vtk_Out_material_idx(Vtk_file, MPM_Mesh.MatIdx, NumParticles);
  }

  /* Print particle velocity */
  if(Out_velocity)
  {
    vtk_Out_vel(Vtk_file, MPM_Mesh.Phi.vel, NumParticles);
  }

  /* Print particle acceleration */
  if(Out_acceleration)
  {
    vtk_Out_acc(Vtk_file, MPM_Mesh.Phi.acc, NumParticles);
  }

  /* Print particle displacement */
  if(Out_displacement)
  {
    vtk_Out_dis(Vtk_file, MPM_Mesh.Phi.dis, NumParticles);
  }

  /* Print particle stress tensor */
  if(Out_stress)
  {
    vtk_Out_Stress(Vtk_file, MPM_Mesh.Phi.Stress, NumParticles);
  }

  /* Print particle eigenvalues of the stress tensor */
  if(Out_eigenvalues_stress)
  {
    vtk_Out_Stress_EV(Vtk_file, MPM_Mesh.Phi.Stress, NumParticles);
  }

  /* Print particle volumetric stress tensor */
  if(Out_volumetric_stress)
  {
    vtk_Out_Stress_P(Vtk_file, MPM_Mesh.Phi.Stress, NumParticles);
  }

  /* Print particle strain tensor */
  if(Out_strain)
  {
    vtk_Out_Stain(Vtk_file, MPM_Mesh.Phi.Strain, NumParticles);
  }

  /* Print particle eigenvalues of the strain tensor */
  if(Out_eigenvalues_strain)
  {
    vtk_Out_Strain_EV(Vtk_file, MPM_Mesh.Phi.Strain, NumParticles);
  }

  /* Print particle deformation gradient */
  if(Out_deformation_gradient)
  {
    vtk_Out_Deformation_Gradient(Vtk_file, MPM_Mesh.Phi.F_n, NumParticles);    
  }

  /* Print particle energy */
  if(Out_energy)
  {
    vtk_Out_Deformation_Energy(Vtk_file, MPM_Mesh.Phi.W, NumParticles);
    vtk_Out_Kinetic_Energy(Vtk_file, MPM_Mesh.Phi.vel, MPM_Mesh.Phi.mass, NumParticles);
  }  


  /* Close the file */
  fclose(Vtk_file);
  
}

/*********************************************************************/

void nodal_results_vtk__InOutFun__(char * Name_File, Mesh ElementMesh, Mask ActiveNodes, Matrix REACTIONS, int TimeStep_i, int ResultsTimeStep)
{

  /* Number of dimensions */
  int Ndim = NumberDimensions;
  FILE * Vtk_file;
  char Name_file_t[80];
  int NumberFields;
  char * FieldsList[MAXW] = {NULL};
  int i_mask;
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
  fprintf(Vtk_file,"POINTS %i double \n",ActiveNodes.Nactivenodes);
  for(int i = 0 ; i<ElementMesh.NumNodesMesh ; i++)
  {
    i_mask = ActiveNodes.Nodes2Mask[i];
    if(i_mask != -1)
      {
        for(int j = 0 ; j<3 ; j++)
          {
            if(j<Ndim)
              {
	               fprintf(Vtk_file,"%lf ",ElementMesh.Coordinates.nM[i][j]);
              }
            else
              {
	               fprintf(Vtk_file,"%lf ",0.0);
              }
          }
        fprintf(Vtk_file,"\n");
      }
  }

  /* Connectivity */
  fprintf(Vtk_file,"CELLS %i %i \n",ActiveNodes.Nactivenodes,Ndim*ActiveNodes.Nactivenodes);
  for(int i = 0 ; i<ActiveNodes.Nactivenodes ; i++)
  {
    fprintf(Vtk_file,"%i %i \n",1,i);
  }

  /* Type of element */
  fprintf(Vtk_file,"CELL_TYPES %i \n",ActiveNodes.Nactivenodes);
  for(int i = 0 ; i<ActiveNodes.Nactivenodes ; i++){
    fprintf(Vtk_file,"%i \n",1);
  }

  /* Point data */
  fprintf(Vtk_file,"POINT_DATA %i \n",ActiveNodes.Nactivenodes);
  fprintf(Vtk_file,"VECTORS %s double \n","X_GC");
  for(int i =  0 ; i<ElementMesh.NumNodesMesh ; i++)
  {
    i_mask = ActiveNodes.Nodes2Mask[i];
    if(i_mask != -1)
      {
        for(int j = 0 ; j<3 ; j++)
          {
            if(j<Ndim)
              {
	               fprintf(Vtk_file,"%lf ",ElementMesh.Coordinates.nM[i][j]);
              }
            else
              {
	               fprintf(Vtk_file,"%lf ",0.0);
              }
          }
        fprintf(Vtk_file,"\n");
      }
  }

//  /* Active nodes */
//  fprintf(Vtk_file,"SCALARS Active_Nod_Active int \n");
//  fprintf(Vtk_file,"LOOKUP_TABLE default \n");
//  for(int i = 0 ; i<ElementMesh.NumNodesMesh ; i++){
//    fprintf(Vtk_file,"%i\n",
//	    ElementMesh.NumParticles[i]);
//  } 
      
	fprintf(Vtk_file,"VECTORS %s double \n","REACTIONS");
	for(int i =  0 ; i<ElementMesh.NumNodesMesh ; i++)
  {
    i_mask = ActiveNodes.Nodes2Mask[i];
    if(i_mask != -1)
      {
	       /* Print the dimensions of the array */
    	   for(int j = 0 ; j<3 ; j++)
           {
	           if(j<Ndim)
                {
	                 fprintf(Vtk_file,"%lf ",REACTIONS.nM[i_mask][j]);
	              }
	           else
                {
            	      fprintf(Vtk_file,"%lf ",0.0);
	              }    
	         }
	         fprintf(Vtk_file,"\n");	
      }
	}
  
  /* Cell data */  
  fprintf(Vtk_file,"CELL_DATA %i \n",ActiveNodes.Nactivenodes);

  /* Close the file */
  fclose(Vtk_file);

}

/*********************************************************************/

static void vtk_Out_X_GC(FILE * Vtk_file, Matrix x_GC, int NumParticles)
{
  int Ndim = NumberDimensions;

  fprintf(Vtk_file,"VECTORS %s double \n","X_GC");
  for(int i =  0 ; i<NumParticles ; i++){
    for(int j = 0 ; j<3 ; j++){
      if(j<Ndim){
        fprintf(Vtk_file,"%lf ",x_GC.nM[i][j]);
      }
      else{
        fprintf(Vtk_file,"%lf ",0.0);
      }
    }
    fprintf(Vtk_file,"\n");
  }
}

/*********************************************************************/

static void vtk_Out_X_EC(FILE * Vtk_file, Matrix x_EC, int NumParticles)
{
  int Ndim = NumberDimensions;

  fprintf(Vtk_file,"VECTORS %s double \n","X_EC");
  for(int i =  0 ; i<NumParticles ; i++){
    for(int j = 0 ; j<3 ; j++){
      if(j<Ndim){
  fprintf(Vtk_file,"%lf ",x_EC.nM[i][j]);
      }
      else{
  fprintf(Vtk_file,"%lf ",0.0);
      }
    }
    fprintf(Vtk_file,"\n");
  }
}

/*********************************************************************/

static void vtk_Out_mass(FILE * Vtk_file, Matrix mass, int NumParticles)
{
  fprintf(Vtk_file,"SCALARS MASS double \n");
  fprintf(Vtk_file,"LOOKUP_TABLE default \n");
  for(int i =  0 ; i<NumParticles ; i++){
    fprintf(Vtk_file,"%lf \n",mass.nV[i]);
  }
}

/*********************************************************************/

static void vtk_Out_density(FILE * Vtk_file, Matrix rho, int NumParticles)
{
  fprintf(Vtk_file,"SCALARS DENSITY double \n");
  fprintf(Vtk_file,"LOOKUP_TABLE default \n");
  for(int i =  0 ; i<NumParticles ; i++){
    fprintf(Vtk_file,"%lf \n",rho.nV[i]); 
  }
}

/*********************************************************************/

static void vtk_Out_damage(FILE * Vtk_file, Matrix chi, int NumParticles)
{
  fprintf(Vtk_file,"SCALARS Damage double \n");
  fprintf(Vtk_file,"LOOKUP_TABLE default \n");
  for(int i =  0 ; i<NumParticles ; i++){
    fprintf(Vtk_file,"%lf \n",chi.nV[i]); 
  }
}

/*********************************************************************/

static void vtk_Out_nodal_idx(FILE * Vtk_file, int * I0, int NumParticles)
{
  fprintf(Vtk_file,"SCALARS ELEM_i integer \n");
  fprintf(Vtk_file,"LOOKUP_TABLE default \n");
  for(int i =  0 ; i<NumParticles ; i++){
    fprintf(Vtk_file,"%i \n",I0[i]);
  }
}

/*********************************************************************/
static void vtk_Out_material_idx(FILE * Vtk_file, int * MatIdx, int NumParticles)
{
  fprintf(Vtk_file,"SCALARS MatIdx integer \n");
  fprintf(Vtk_file,"LOOKUP_TABLE default \n");
  for(int i =  0 ; i<NumParticles ; i++){
    fprintf(Vtk_file,"%i \n",MatIdx[i]);
  }
}

/*********************************************************************/

static void vtk_Out_vel(FILE * Vtk_file, Matrix vel, int NumParticles)
{
  int Ndim = NumberDimensions;

    fprintf(Vtk_file,"VECTORS VELOCITY double \n");
  for(int i =  0 ; i<NumParticles ; i++){
    for(int j = 0 ; j<3 ; j++){
      if(j<Ndim){
  fprintf(Vtk_file,"%lf ",vel.nM[i][j]);
      }
      else{
  fprintf(Vtk_file,"%lf ",0.0);
      }
    }
    fprintf(Vtk_file,"\n");
  }
}

/*********************************************************************/

static void vtk_Out_acc(FILE * Vtk_file, Matrix acc, int NumParticles)
{
  int Ndim = NumberDimensions;

    fprintf(Vtk_file,"VECTORS ACCELERATION double \n");
  for(int i =  0 ; i<NumParticles ; i++){
    for(int j = 0 ; j<3 ; j++){
      if(j<Ndim){
  fprintf(Vtk_file,"%lf ",acc.nM[i][j]);
      }
      else{
  fprintf(Vtk_file,"%lf ",0.0);
      }
    }
    fprintf(Vtk_file,"\n");
  }
}

/*********************************************************************/

static void vtk_Out_dis(FILE * Vtk_file, Matrix dis, int NumParticles)
{
  int Ndim = NumberDimensions;

    fprintf(Vtk_file,"VECTORS DISPLACEMENT double \n");
  for(int i =  0 ; i<NumParticles ; i++){
    for(int j = 0 ; j<3 ; j++){
      if(j<Ndim){
  fprintf(Vtk_file,"%lf ",dis.nM[i][j]);
      }
      else{
  fprintf(Vtk_file,"%lf ",0.0);
      }
    }
    fprintf(Vtk_file,"\n");
  }
}

/*********************************************************************/

static void vtk_Out_Stress(FILE * Vtk_file, Matrix Stress, int NumParticles)
{
  int Ndim = NumberDimensions;

  fprintf(Vtk_file,"TENSORS STRESS double \n");
  for(int i =  0 ; i<NumParticles ; i++){
    for(int j = 0 ; j<3 ; j++){
      for(int k = 0 ; k<3 ; k++){
  if((j<Ndim) && (k<Ndim)){
    fprintf(Vtk_file,"%lf ",Stress.nM[i][j*Ndim+k]);
  }
  else{
    fprintf(Vtk_file,"%lf ",0.0);
  }
      }
      fprintf(Vtk_file,"\n");
    }
    fprintf(Vtk_file,"\n");
  }
}

/*********************************************************************/

static void vtk_Out_Stress_EV(FILE * Vtk_file, Matrix Stress, int NumParticles)
{
  int Ndim = NumberDimensions;
  Tensor Stress_p;
  Tensor EV_Stress_p; 

  fprintf(Vtk_file,"VECTORS STRESS_EV double \n");
  for(int i =  0 ; i<NumParticles ; i++){
    Stress_p = memory_to_tensor__TensorLib__(Stress.nM[i], 2);
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
}

/*********************************************************************/

static void vtk_Out_Stress_P(FILE * Vtk_file, Matrix Stress, int NumParticles)
{
  int Ndim = NumberDimensions;
  double P_GP; /* Trace of the stress tensor (volumetric) */

  fprintf(Vtk_file,"SCALARS P double \n");
  fprintf(Vtk_file,"LOOKUP_TABLE default \n");
  for(int i =  0 ; i<NumParticles ; i++){
    P_GP = 0;
    for(int j = 0 ; j<3 ; j++ ){
      if(j < Ndim){
        P_GP += Stress.nM[i][j];
      }
      else{
        P_GP = 0;
      }
    }    
    fprintf(Vtk_file,"%lf \n",P_GP/Ndim);
  }

}

/*********************************************************************/

static void vtk_Out_Stain(FILE * Vtk_file, Matrix Strain, int NumParticles)
{
    int Ndim = NumberDimensions;

    fprintf(Vtk_file,"TENSORS STRAIN double \n");
  for(int i =  0 ; i<NumParticles ; i++){
    for(int j = 0 ; j<3 ; j++){
      for(int k = 0 ; k<3 ; k++){
  if((j<Ndim) && (k<Ndim)){
    fprintf(Vtk_file,"%lf ",Strain.nM[i][j*Ndim+k]);
  }
  else{
    fprintf(Vtk_file,"%lf ",0.0);
  }
      }
      fprintf(Vtk_file,"\n");
    }
    fprintf(Vtk_file,"\n");
  }
}

/*********************************************************************/

static void vtk_Out_Strain_EV(FILE * Vtk_file, Matrix Strain, int NumParticles)
{
  int Ndim = NumberDimensions;
  Tensor Strain_p;
  Tensor EV_Strain_p; 

    fprintf(Vtk_file,"VECTORS STRAIN_EV double \n");
  for(int i =  0 ; i<NumParticles ; i++){
    Strain_p = memory_to_tensor__TensorLib__(Strain.nM[i], 2);
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
}

/*********************************************************************/

static void vtk_Out_Deformation_Gradient(FILE * Vtk_file, Matrix F_n, int NumParticles)
{
  int Ndim = NumberDimensions;

    fprintf(Vtk_file,"TENSORS DEFORMATION-GRADIENT double \n");
  for(int i =  0 ; i<NumParticles ; i++){
    for(int j = 0 ; j<3 ; j++){
      for(int k = 0 ; k<3 ; k++){
  if((j<Ndim) && (k<Ndim)){
    fprintf(Vtk_file,"%lf ",F_n.nM[i][j*Ndim+k]);
  }
  else{
    fprintf(Vtk_file,"%lf ",0.0);
  }
      }
      fprintf(Vtk_file,"\n");
    }
    fprintf(Vtk_file,"\n");
  }
}

/*********************************************************************/

static void vtk_Out_Deformation_Energy(FILE * Vtk_file, Matrix W, int NumParticles)
{
  fprintf(Vtk_file,"SCALARS W double \n");
  fprintf(Vtk_file,"LOOKUP_TABLE default \n");
  for(int i =  0 ; i<NumParticles ; i++)
  {
    fprintf(Vtk_file,"%lf \n",W.nV[i]); 
  }
}

/*********************************************************************/

static void vtk_Out_Kinetic_Energy(FILE * Vtk_file, Matrix vel, Matrix mass, int NumParticles)
{
  int Ndim = NumberDimensions;
  double K_p;

  fprintf(Vtk_file,"SCALARS K double \n");
  fprintf(Vtk_file,"LOOKUP_TABLE default \n");
  for(int i =  0 ; i<NumParticles ; i++)
  {
    K_p = 0 ;
    for(int j = 0 ; j<Ndim ; j++)
    {
      K_p += DSQR(vel.nM[i][j]);
    }
    fprintf(Vtk_file,"%lf \n",0.5*K_p*mass.nV[i]); 
  }
}

/*********************************************************************/
