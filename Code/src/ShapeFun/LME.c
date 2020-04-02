#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "grams.h"

/***********************************************/
/************* Local Maximum-Entropy **************/
/**************************************************/

/*!
  Shape functions based in :
  "" Local maximum-entropy approximation schemes : a seamless 
  bridge between finite elements and meshfree methods ""
  by M.Arroyo and M.Ortiz, 2006.

  Here we employ the same nomenclature as in the paper. With the single
  different of the "l" variable wich represents the distances between the
  evaluation point and the neighborhood nodes.

  List of functions :
  - LME_Init_lambda
  - LME_lambda_NR
  - LME_fa
  - LME_p
  - LME_r
  - LME_J
  - LME_dp
  - LME_Tributary_Nodes
*/

/****************************************************************************/

void LME_Initialize(GaussPoint MPM_Mesh, Mesh FEM_Mesh)
{

  int Ndim = 3;

  /* Variables for the GP coordinates */  
  Matrix X_GC_GP = MatAssign(Ndim,1,NAN,NULL,NULL);
  Matrix X_EC_GP = MatAssign(Ndim,1,NAN,NULL,NULL);

  /* Variables for the poligon */
  int NumVertex;
  int * Poligon_Connectivity;
  Matrix Poligon_Coordinates;
  ChainPtr ListNodes_I;

  /* Auxiliar variables for LME */
  Matrix lambda_GP = /* Lagrange multipliers */
    MatAssign(Ndim,1,NAN,NULL,NULL);
  Matrix Delta_Xip; /* Distance from GP to the nodes */
  Matrix Dist;
  Matrix Beta_GP = /* Tunning parameter */
      MatAssign(Ndim,1,NAN,NULL,NULL);
  int NumNodes_GP; /* Number of neibourghs */
  int * ListNodes; /* List of nodes */
  int I_iGP; /* Iterator for the neibourghs */

  /* */
  ChainPtr List_Ord, List_Dis;
  int * List;

  /* Auxiliar variables to initialize lambda */
  Matrix A = MatAlloc(Ndim,Ndim);
  Matrix B = MatAlloc(Ndim,1);

  /* 1º Set to zero the active/non-active node, and the GPs in each 
   element */
  for(int i = 0 ; i<FEM_Mesh.NumNodesMesh ; i++){
    FEM_Mesh.ActiveNode[i] = 0;
  }
  
  for(int i = 0 ; i<FEM_Mesh.NumElemMesh ; i++){
    FreeChain(FEM_Mesh.GPsElements[i]);
    FEM_Mesh.GPsElements[i] = NULL;
  }

  for(int i = 0 ; i<MPM_Mesh.NumGP ; i++){

    /* 2º Assign the value to this auxiliar pointer */ 
    X_GC_GP.nV = MPM_Mesh.Phi.x_GC.nM[i];

    /* 3º Initialize Beta and assign to an auxiliar pointer */
    for(int j = 0 ; j<Ndim ; j++){
      MPM_Mesh.Beta.nM[i][j] =
	gamma_LME/(FEM_Mesh.DeltaX*FEM_Mesh.DeltaX);
    }
    Beta_GP.nV = MPM_Mesh.Beta.nM[i];


    for(int j = 0 ; j<FEM_Mesh.NumElemMesh ; j++){

      /* 4º Connectivity of the Poligon */
      NumVertex = FEM_Mesh.NumNodesElem[j];
      Poligon_Connectivity =
	ChainToArray(FEM_Mesh.Connectivity[j],NumVertex);

      /* 5º Get the coordinates of the element */
      Poligon_Coordinates =
	ElemCoordinates(FEM_Mesh,Poligon_Connectivity,NumVertex);
      
      /* 6º Check out if the GP is in the Element */
      if(InOut_Poligon(X_GC_GP,Poligon_Coordinates) == 1){

	/* 7º Asign to the GP a element in the background mesh, just for 
	   searching porpuses */
	MPM_Mesh.Element_id[i] = j;
	PushNodeTop(&FEM_Mesh.GPsElements[j],i);

	/* 8º If the GP is in the element, get its natural coordinates */
	X_EC_GP.nV = MPM_Mesh.Phi.x_EC.nM[i];
	Q4_X_to_Xi(X_EC_GP,X_GC_GP,Poligon_Coordinates);

	/* 9º Get list of nodes near to the GP */
	FreeChain(MPM_Mesh.ListNodes[i]);
	MPM_Mesh.ListNodes[i] = NULL;
 
	/* 10º Calculate connectivity */
	MPM_Mesh.ListNodes[i] =
	  LME_Tributary_Nodes(X_GC_GP,Beta_GP,
			      MPM_Mesh.Element_id[i],FEM_Mesh);
    
	/* 11º Calculate number of nodes */
	MPM_Mesh.NumberNodes[i] = LenghtChain(MPM_Mesh.ListNodes[i]);

	/* Generate nodal distance list */
	NumNodes_GP = MPM_Mesh.NumberNodes[i];
	ListNodes = ChainToArray(MPM_Mesh.ListNodes[i],NumNodes_GP);
	Delta_Xip = MatAlloc(NumNodes_GP,Ndim);
	/* Allocate the distance vector of each node */
	Dist = MatAllocZ(NumNodes_GP,1);
	/* Fill both vectors */
	for(int k = 0 ; k<NumNodes_GP ; k++){
	  I_iGP = ListNodes[k];
	  for(int l = 0 ; l<Ndim ; l++){
	    Delta_Xip.nM[k][l] =
	      X_GC_GP.nV[l]-
	      FEM_Mesh.Coordinates.nM[I_iGP][l];
	    /* Distance modulus */
	    Dist.nV[k] += Delta_Xip.nM[k][l]*Delta_Xip.nM[k][l];
	  }
	  Dist.nV[k] = pow(Dist.nV[k],0.5);
	}
	free(ListNodes);

	/* Update the value of beta */
	Beta_GP = LME_Beta(Beta_GP, Delta_Xip, gamma_LME);

	/* Ordenate distances  */
	List_Ord = NULL, List_Dis = NULL;
	List_Dis = RangeChain(0,NumNodes_GP-1);
	OrderList(&List_Ord,&List_Dis,Dist);

	/* Transform the list in to an array */
	List = ChainToArray(List_Ord,NumNodes_GP);

	/* Fill A */
	A.nM[0][0] = Delta_Xip.nM[List[1]][0] - Delta_Xip.nM[List[0]][0];
	A.nM[0][1] = Delta_Xip.nM[List[1]][1] - Delta_Xip.nM[List[0]][1];
	A.nM[1][0] = Delta_Xip.nM[List[2]][0] - Delta_Xip.nM[List[0]][0];
	A.nM[1][1] = Delta_Xip.nM[List[2]][1] - Delta_Xip.nM[List[0]][1];
  
	/* Fill B */
	B.nV[0] = -Beta_GP.nV[0]*(pow(Dist.nV[List[0]],2) -
				  pow(Dist.nV[List[1]],2));
	B.nV[1] = -Beta_GP.nV[1]*(pow(Dist.nV[List[0]],2) -
				  pow(Dist.nV[List[2]],2));

	/* /\* Check the conditioning number of A *\/ */
	/* if (fabs(Cond_Mat(A,TOL_NR)) > 10){ */
	/*   printf(" %s (%s %i) : %s \n", */
	/* 	 "Error in LME_Initialize","GP",i, */
	/* 	 "A is near to singular matrix"); */
	/*   exit(0); */
	/* } */
	/* Initialize lambda */
	lambda_GP.nV = MPM_Mesh.lambda.nM[i];
	/* lambda_GP = Solve_Linear_Sistem(A,B,lambda_GP); */
	
	/* Calculate lagrange multipliers with Newton-Rapson */
	lambda_GP = LME_lambda_NR(Delta_Xip, lambda_GP, Beta_GP);

	/* Free memory */
	FreeMat(Delta_Xip), FreeChain(List_Ord), free(List);
	
	/* 9º Active those nodes that interact with the GP */
	ListNodes_I = MPM_Mesh.ListNodes[i];
	while(ListNodes_I != NULL){
	  FEM_Mesh.ActiveNode[ListNodes_I->I] += 1;
	  ListNodes_I = ListNodes_I->next; 
	}

	/* 10º Free memory and go for the next GP */
	free(Poligon_Connectivity);
	FreeMat(Poligon_Coordinates);
	break;
	
      }
      
      /* 11º Free memory */
      free(Poligon_Connectivity);
      FreeMat(Poligon_Coordinates);
      
    }

  }

  /* 12º Free memory */
  FreeMat(A), FreeMat(B);
  
}

/****************************************************************************/

Matrix LME_Beta(Matrix Beta, Matrix l, double Gamma)
/*!
  Function to update the value of beta
*/
{

  int Ndim = 3;
  int NumNodes_GP = l.N_rows;
  double h = 0;
  Matrix l_GP_I = MatAssign(Ndim,1,NAN,NULL,NULL);
  
  /* Get the mean distande */
  for(int i = 0 ; i<NumNodes_GP ; i++){
    l_GP_I.nV = l.nM[i];
    h += Norm_Mat(l_GP_I,2);
  }
  h = h/NumNodes_GP;

  /* Fill Beta */
  for(int j = 0 ; j<Ndim ; j++){
    Beta.nV[j] = Gamma/(h*h);
  }

  return Beta;
}

/****************************************************************************/

Matrix LME_lambda_NR(Matrix l, Matrix lambda, Matrix Beta)
/*!
  Get the lagrange multipliers "lambda" (1 x dim) for the LME 
  shape function. The numerical method for that is the Newton-Rapson.

  Input parameters :
  -> l : Matrix with the distances to the
  neighborhood nodes (neighborhood x dim).
  -> lambda : Initial value of the
  lagrange multipliers (1 x dim).
  -> Beta : Tunning parameter (scalar).
  -> h : Grid spacing (scalar).
  -> TOL_zero : Tolerance for Newton-Rapson.
*/
{
  /* Definition of some parameters */
  int MaxIter = 100;
  int Ndim = 3;
  int NumIter = 0; /* Iterator counter */
  double norm_r = 10; /* Value of the norm */
  Matrix p; /* Shape function vector */
  Matrix r; /* Gradient of log(Z) */
  Matrix J; /* Hessian of log(Z) */
  Matrix D_lambda; /* Increment of lambda */
 
  /* Start with the Newton-Rapson method */
  while(norm_r > TOL_NR){
	
    /* Get vector with the shape functions evaluated in the nodes */
    p = LME_p(l,lambda,Beta);

    /* Get the gradient of log(Z) and its norm */
    r = LME_r(l,p), norm_r = Norm_Mat(r,2);

    /* Get the Hessian of log(Z) and update it with +||r||*I 
       according with Dennis M.Kochmann et al. 2019 (CMAME) */    
    J = LME_J(l,p,r);
    for(int i = 0 ; i<Ndim ; i++){
      J.nM[i][i] += norm_r;
    }

    /* Check the conditioning number of the Hessian */
    if (fabs(Cond_Mat(J,TOL_NR)) > 10){
      printf(" %s (%s %i) : %s \n",
    	     "Error in LME_lambda_NR","Iter",NumIter,
    	     "The Hessian is near to singular matrix");
      exit(0);
    }
    
    /* Get the increment of lambda */
    D_lambda = Solve_Linear_Sistem(J,r,MatAllocZ(Ndim,1));

    /* Update the value of lambda */
    for(int i = 0 ; i<Ndim ; i++){
      lambda.nV[i] -= D_lambda.nV[i];
    }

    /* Free memory */
    FreeMat(p), FreeMat(r), FreeMat(J), FreeMat(D_lambda);
    
    /* Update the number of iterations */
    NumIter ++;
    if(NumIter >= MaxIter){
      printf("%s : %s %i/%i %s \n",
      	     "Warning in LME_lambda",
      	     "No convergence in",NumIter,MaxIter,"iterations");
      printf("%s : %f \n",
	     "Error",norm_r);
      exit(0);
    }
  
  }
  
  /* Once the stopping criteria is reached, 
     return the lagrange multipliers value */
  return lambda;
}

/****************************************************************************/

double LME_fa(Matrix la, Matrix lambda, Matrix Beta)
/*!
  Output :
  -> fa : the function fa that appear in [1] (scalar).
  Input parameters :
  -> la : Matrix with the distance to the neighborhood node ''a'' (1 x dim).
  -> lambda : Initial value of the lagrange multipliers (dim x 1).
  -> Gamma : Tunning parameter (scalar).
*/
{  
  int Ndim = 3;
  double fa = 0;

  for(int i = 0 ; i<Ndim ; i++){
    fa += - Beta.nV[i]*la.nV[i]*la.nV[i] + la.nV[i]*lambda.nV[i];
  }
    
  /* Return the value of fa */
  return fa;
}

/****************************************************************************/

Matrix LME_p(Matrix l, Matrix lambda, Matrix Beta)
/*!
  Get the value of the shape function "pa" (1 x neighborhood) in the
  neighborhood nodes.

  Input parameters :
  -> l : Matrix with the distances to the
  neighborhood nodes (neiborghood x dim).
  -> lambda : Initial value of the lagrange multipliers (1 x dim).
  -> Beta : Tunning parameter (scalar).
*/
{
  
  /* Definition of some parameters */
  int N_a = l.N_rows;
  int Ndim = 3;
  double Z = 0, Z_m1 = 0;
  Matrix p = /* Vector with the values of the shape-function in the nodes */
    MatAlloc(1,N_a); 
  Matrix la = /* Distance to the neighbour (x-x_a) */
    MatAssign(1,Ndim,NAN,NULL,NULL);

  /* Get Z and the numerator */
  for(int a = 0 ; a<N_a ; a++){
    la.nV = l.nM[a];
    p.nV[a] = exp(LME_fa(la,lambda,Beta));
    Z += p.nV[a];
  }

  /* Get the inverse of Z */
  Z_m1 = (double)1/Z;

  /* Divide by Z and get the final value */
  for(int a = 0 ; a<N_a ; a++){
    p.nV[a] *= Z_m1;
  }
  
  /* Return the value of the shape function */  
  return p;
}

/****************************************************************************/

Matrix LME_r(Matrix l, Matrix p)
/*!
  Get the gradient "r" (dim x 1) of the function log(Z) = 0.
  Input parameters :
  -> l : Matrix with the distances to the 
  neighborhood nodes (neighborhood x dim).
  -> p : Shape function value in the
  neighborhood nodes (1 x neighborhood).
*/
{  
  /* Definition of some parameters */
  int N_a = l.N_rows;
  int Ndim = 3;
  Matrix r /* Value of the gradient */
    = MatAllocZ(Ndim,1);

  /* Fill ''r'' */
  for(int i = 0 ; i<Ndim ; i++){
    for(int a = 0 ; a<N_a ; a++){
      r.nV[i] += p.nV[a]*l.nM[a][i];
    }
  }

  /* Return the value of the gradient */
  return r;
}

/****************************************************************************/

Matrix LME_J(Matrix l, Matrix p, Matrix r)
/*!
  Get the Hessian "J" (dim x dim) of the function log(Z) = 0.
  Input parameters :
  -> l : Matrix with the distances to the
  neighborhood nodes (neighborhood x dim).
  -> p : Shape function value in the
  neighborhood nodes (neighborhood x 1).
  -> r : Gradient of log(Z) (dim x 1).
*/
{  
  /* Definition of some parameters */
  int N_a = l.N_rows;
  int Ndim = 3;
  Matrix J; /* Hessian definition */
  
  /* Allocate Hessian */
  J = MatAllocZ(Ndim,Ndim);

  /* Fill the Hessian */
  for(int i = 0 ; i<Ndim ; i++){
    for(int j = 0 ; j<Ndim ; j++){
      for(int a = 0 ; a<N_a ; a++){
	/* Get the first component of the Hessian looping
	   over the neighborhood nodes. */
	J.nM[i][j] += p.nV[a]*l.nM[a][i]*l.nM[a][j];
      }
      /* Get the second value of the Hessian */
      J.nM[i][j] -= r.nV[i]*r.nV[j];
    }
  }

  /* Return the value of the Hessian */
  return J;
}

/****************************************************************************/

Matrix LME_dp(Matrix l, Matrix p)
/*!
  Value of the shape function gradient "dp" (dim x neighborhood) in 
  the neighborhood nodes.
  Input parameters :
  -> l : Matrix with the distances to the
  neighborhood nodes (neighborhood x dim).
  -> p : Shape function value in the
  neighborhood nodes (neighborhood x 1).
*/
{  
  /* Definition of some parameters */
  int N_a = l.N_rows;
  int Ndim = 3;
  Matrix dp = MatAllocZ(N_a,Ndim);
  Matrix r; /* Gradient of log(Z) */
  Matrix J; /* Hessian of log(Z) */
  Matrix Jm1; /* Inverse of J */
  Matrix Jm1_la; /* Auxiliar vector */
  Matrix la = /* Distance to the neighbour (x-x_a) */
    MatAssign(Ndim,1,NAN,NULL,NULL); 
  
  /* Get the Gradient and the Hessian of log(Z) */
  r = LME_r(l,p);
  J = LME_J(l,p,r);
  
  /* Inverse of the Hessian */
  Jm1 = Get_Inverse(J);
  
  /* Free memory */
  FreeMat(r), FreeMat(J);

  /* Fill the gradient for each node */
  for(int a = 0 ; a<N_a ; a++){
    la.nV = l.nM[a]; 
    Jm1_la = Scalar_prod(Jm1,la);    
    for(int i = 0 ; i<Ndim ; i++){
      dp.nM[a][i] = - p.nV[a]*Jm1_la.nV[i];
    }
    FreeMat(Jm1_la);
  }

  /* Free memory */
  FreeMat(Jm1);
  
  /* Return the value of the shape function gradient */  
  return dp;  
}

/****************************************************************************/

ChainPtr LME_Tributary_Nodes(Matrix X_GP, Matrix Beta,
			     int Elem_GP, Mesh FEM_Mesh){

  int Ndim = 3;
  Matrix Distance; /* Distance between node and GP */
  Matrix X_I = MatAssign(Ndim,1,NAN,NULL,NULL);
  ChainPtr * Table_Elem = NULL;
  ChainPtr Triburary_Nodes = NULL;
  ChainPtr List_Nodes = NULL;
  ChainPtr * Table_ElemNodes = NULL;
  ChainPtr Triburary_Elements = NULL;
  ChainPtr iPtr = NULL;
  int NumNodesElem = /* Number of nodes of the element */
    FEM_Mesh.NumNodesElem[Elem_GP];
  int Num_Elem;
  int * List_Elements;
  int * NodesElem = /* List of nodes of the element */
    ChainToArray(FEM_Mesh.Connectivity[Elem_GP],NumNodesElem);
  double Ra = /* Get the search radius */
    sqrt(-log(TOL_lambda)/Beta.nV[0]);


  /* Chain with the tributary elements, this is the list of element near the
     gauss point, including where it is */

  /* Iterate in the list and select the union of the sets of nodes */
  Table_Elem = malloc(NumNodesElem*sizeof(ChainPtr));
  for(int i = 0 ; i<NumNodesElem ; i++){
    Table_Elem[i] = FEM_Mesh.NodeNeighbour[NodesElem[i]];
  }
  Triburary_Elements = ChainUnion(Table_Elem,NumNodesElem);
  /* Free memory */
  free(NodesElem);
  free(Table_Elem);
  Table_Elem = NULL;
  
  /* List with the tributary nodes */
  Num_Elem = LenghtChain(Triburary_Elements);
  List_Elements = ChainToArray(Triburary_Elements,Num_Elem);

  /* Free the chain with the tributary elements */
  FreeChain(Triburary_Elements);
  
  /* Fill the chain with the preliminary tributary nodes */
  Table_ElemNodes = malloc(Num_Elem*sizeof(ChainPtr));
  for(int i = 0 ; i<Num_Elem ; i++){
    Table_ElemNodes[i] = FEM_Mesh.Connectivity[List_Elements[i]];
  }

  List_Nodes = ChainUnion(Table_ElemNodes,Num_Elem);
  
  /* Free the array wit the list of tributary elements */
  free(List_Elements);
  free(Table_ElemNodes);
  Table_ElemNodes = NULL;
  
  /* Initialize the iterator to iterate over the list of tributary nodes */
  iPtr = List_Nodes;

  /* Loop over the chain with the tributary nodes */
  while(iPtr != NULL){

    /* Assign to a pointer the coordinates of the nodes */
    X_I.nV = FEM_Mesh.Coordinates.nM[iPtr->I];

    /* Get a vector from the GP to the node */
    Distance = Sub_Mat(X_GP,X_I);

    /* If the node is near the GP push in the chain */
    if(Norm_Mat(Distance,2) <= Ra){
      PushNodeTop(&Triburary_Nodes,iPtr->I);
    }

    /* Free memory of the distrance vector */
    FreeMat(Distance);

    /* Update pointer index */
    iPtr = iPtr->next;
  }
  /* Free memory */
  FreeChain(List_Nodes);

  return Triburary_Nodes;
}

