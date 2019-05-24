#define USE_CONST
#include <stdlib.h>
#include "gidpost.h"

typedef struct {
  int id;
  int uctrlp;
  int vctrlp;
  double* ctrlpoints;
} SNurbSurface;

void GenNurb(SNurbSurface * nsurf,int id,int uctrlp,int vctrlp){
  nsurf->id = id;
  nsurf->uctrlp = uctrlp;
  nsurf->vctrlp = vctrlp;
  nsurf->ctrlpoints=NULL;//not used in this example
}

float Random() {
  return rand()/(float)(RAND_MAX);
}

int main() {
  int output_format=2;//0 ASCII, 2 binary, 3 HDF5
  const char* analysis="Analysis_example_IGA";
  const char* velocity_component_names[]={"Velocity-x","Velocity-y","Velocity-z"};
  double* values;
  int i,nval;
  int i_surf,n_surf=2;
  SNurbSurface nsurfs[2];
  double time_step=1.0;//for results along the time
  GenNurb(&nsurfs[0],1,3,2); /* id=1, num control points=3x2 */
  GenNurb(&nsurfs[1],2,2,2); /* id=2, num control points=2x2 */
  /* ids and dimensions of these surfaces must match the ones defined in the file IGA_test_1.geo , 
     associated with the file IGA_test_1.post.res with random results that will be created by this example */

  GiD_PostInit();
  
  if(output_format==0){
    GiD_OpenPostResultFile( "IGA_test_1.post.res", GiD_PostAscii);
  } else if(output_format==2){
    GiD_OpenPostResultFile( "IGA_test_1.post.res", GiD_PostBinary);
  } else if(output_format==3){
    GiD_OpenPostResultFile( "IGA_test_1.post.res", GiD_PostHDF5);
  } else {
    return 1;
  }

  GiD_BeginResultHeader("Pressure",analysis,time_step,GiD_Scalar,GiD_OnNurbsSurface,NULL);
  GiD_ResultUnit("Pa");
  for(i_surf=0;i_surf<n_surf;i_surf++){
    nval=nsurfs[i_surf].uctrlp*nsurfs[i_surf].vctrlp;
    values = (double*)malloc(nval*sizeof(double));
    for(i=0;i<nval;i++){
      values[i]=Random();
    }
    GiD_WriteNurbsSurface(nsurfs[i_surf].id,nval,values);
    free(values);
  }
  GiD_EndResult();

  //vector 3 components
  GiD_BeginResultHeader("Velocity",analysis,time_step,GiD_Vector,GiD_OnNurbsSurface,NULL);
  GiD_ResultComponents(3,velocity_component_names);
  GiD_ResultUnit("m/s");
  for(i_surf=0;i_surf<n_surf;i_surf++){
    nval=nsurfs[i_surf].uctrlp*nsurfs[i_surf].vctrlp;
    values=(double*)malloc(nval*3*sizeof(double));
    for(i=0;i<nval*3;i++){
      values[i]=Random();
    }
    GiD_WriteNurbsSurfaceVector(nsurfs[i_surf].id,nval,3,values);
    free(values);
  }
  GiD_EndResult();

  GiD_ClosePostResultFile();
  GiD_PostDone();
  return 0;
}
