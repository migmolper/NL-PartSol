#include "nl-partsol.h"


/*********************************************************************/

Matrix body_loads__Particles__(Load * B, int NumLoads, int NumGP, int TimeStep)
/*
  Evaluate body forces in a time step.
*/
{
  Matrix Body_Forces_t = allocZ__MatrixLib__(NumberDimensions,NumGP);
  int GP_Force;
  
  if(NumLoads>0){
    for(int i = 0 ; i<NumLoads; i++){
      for(int j = 0 ; j<B[i].NumNodes ; j++){
	GP_Force = B[i].Nodes[j];
	for(int k = 0 ; k<B[i].Dim ; k++){
	  if( (B[i].Dir[k] == 1) ||
	      (B[i].Dir[k] == -1)){
	    if( (TimeStep < 0) || (TimeStep > B[i].Value[k].Num)){
	      printf("%s : %s\n",
		     "Error in body_loads__Particles__()",
		     "The time step is out of the curve !!");
	      exit(EXIT_FAILURE);
	    }
	    Body_Forces_t.nM[k][GP_Force] +=
	      B[i].Value[k].Fx[TimeStep]*
	      (double)B[i].Dir[k];
	  }
	}
      }
    }
  }

  return Body_Forces_t;

}

/*********************************************************************/

Matrix contact_loads__Particles__(Load * F, int NumLoads, int NumGP, int TimeStep)
/*
  Evaluate contact forces in a time step
 */
{
  int GP_Force;
  Matrix Contact_Forces_t = allocZ__MatrixLib__(NumberDimensions,NumGP);
  
  if(NumLoads>0){
    for(int i = 0 ; i<NumLoads; i++){
      for(int j = 0 ; j<F[i].NumNodes ; j++){
	GP_Force = F[i].Nodes[j];
	for(int k = 0 ; k<F[i].Dim ; k++){
	  if( (F[i].Dir[k] == 1) ||
	      (F[i].Dir[k] == -1)){
	    if( (TimeStep < 0) || (TimeStep > F[i].Value[k].Num)){
	      printf("%s : %s \n",
		     "Error in contact_loads__Particles__()",
		     "The time step is out of the curve !!");
	      exit(EXIT_FAILURE);
	    }
	    Contact_Forces_t.nM[k][GP_Force] +=
	      F[i].Value[k].Fx[TimeStep]*
	      (double)F[i].Dir[k];
	  }
	}
      }
    }
  }

  return Contact_Forces_t;
}

/**********************************************************************/
