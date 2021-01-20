#include "nl-partsol.h"


/**************************************************************/

Mask generate_NodalMask__MeshTools__(Mesh FEM_Mesh)
{
  int Nnodes = FEM_Mesh.NumNodesMesh;
  int Nactivenodes = 0;
  int * Nodes2Mask = (int *)Allocate_ArrayZ(Nnodes,sizeof(int));
  Mask M;

  for(int A = 0 ; A<Nnodes ; A++)
    {
      
      if(FEM_Mesh.NumParticles[A] > 0)
      {
        Nodes2Mask[A] = Nactivenodes;
        Nactivenodes++;
      }
      else
      {
        Nodes2Mask[A] = - 1;
      }
    }

  M.Nactivenodes = Nactivenodes;
  M.Nodes2Mask = Nodes2Mask;  
 
  return M;
}

/**************************************************************/

Mask generate_Mask_for_static_condensation__MeshTools__(Mask ActiveNodes, Mesh FEM_Mesh)
{

  /* 
    Define auxilar variables 
  */
  int Ndof = NumberDOF;
  int Nnodes_mask = ActiveNodes.Nactivenodes;
  int Order = Nnodes_mask*Ndof;
  int Number_of_BCC = FEM_Mesh.Bounds.NumBounds;
  int NumNodesBound; /* Number of nodes of the bound */
  int NumDimBound; /* Number of dimensions */
  int Id_BCC; /* Index of the node where we apply the BCC */
  int Id_BCC_mask;
  int Id_BCC_mask_k;

  /*
    Generate mask for the static condensation.
  */
  int Nactivenodes = 0;
  int * Nodes2Mask = (int *)Allocate_ArrayZ(Order,sizeof(int));
  Mask Free_and_Restricted_Dofs;

  /* 
    Loop over the the boundaries to find the constrained dofs
  */
  for(int i = 0 ; i<Number_of_BCC ; i++)
    {

      /* 
        Get the number of nodes of this boundary
      */
      NumNodesBound = FEM_Mesh.Bounds.BCC_i[i].NumNodes;

      /* 
        Get the number of dimensions where the BCC it is applied 
      */
      NumDimBound = FEM_Mesh.Bounds.BCC_i[i].Dim;

      for(int j = 0 ; j<NumNodesBound ; j++)
        {
          /* 
            Get the index of the node 
          */
          Id_BCC = FEM_Mesh.Bounds.BCC_i[i].Nodes[j];
          Id_BCC_mask = ActiveNodes.Nodes2Mask[Id_BCC];

          /*
            If the boundary condition is under an active node 
          */
          if(Id_BCC_mask != -1)
          {
            /* 
              Loop over the dimensions of the boundary condition 
            */
            for(int k = 0 ; k<NumDimBound ; k++)
              {

                /* 
                  Apply only if the direction is active
                */
                if(FEM_Mesh.Bounds.BCC_i[i].Dir[k] == 1)
                  {
                    Id_BCC_mask_k = Id_BCC_mask*Ndof + k;
                    Nodes2Mask[Id_BCC_mask_k] = -1;
                  }
              }
          }

        }    
    }

  /* 
    Generate mask using the location of the constrained dofs
  */
  for(int A_i = 0 ; A_i<Order ; A_i++)
    {
     if(Nodes2Mask[A_i] != -1)
      {
        Nodes2Mask[A_i] = Nactivenodes;
        Nactivenodes++;
      } 
    }


  /*
    Assign information to the output
  */
  Free_and_Restricted_Dofs.Nactivenodes = Nactivenodes;
  Free_and_Restricted_Dofs.Nodes2Mask = Nodes2Mask;  

  return Free_and_Restricted_Dofs;
}


/**************************************************************/


Matrix get_set_field__MeshTools__(Matrix Field, Element Nodes_p, Mask ActiveNodes)
/*
  This function performs two operations. First takes the nodal connectivity of the particle, 
  and translate it to the mask numeration. Second, generate a Matrix with the nodal values.
  To help in the future computations. Nodal data is substracted in the shape (nodesxndofs).
 */
{
  int Nnodes = Nodes_p.NumberNodes;
  int Ndim = Field.N_cols;
  Matrix Field_Ap = allocZ__MatrixLib__(Nnodes,Ndim);
  int Ap;
  int A_mask;

  if(Ndim > 1)
  {
   for(int A = 0 ; A<Nnodes ; A++)
   {
     
	    /* 
	       Get the node in the mass matrix with the mask
	    */
     Ap = Nodes_p.Connectivity[A];
     A_mask = ActiveNodes.Nodes2Mask[Ap];
     
     for(int i = 0 ; i<Ndim ; i++)
     {
      Field_Ap.nM[A][i] = Field.nM[A_mask][i];
    }
  }
}
else
{
	for(int A = 0 ; A<Nnodes ; A++)
 {
	 /* 
	   Get the node in the mass matrix with the mask
	 */
   Ap = Nodes_p.Connectivity[A];
   A_mask = ActiveNodes.Nodes2Mask[Ap];
   
   Field_Ap.nV[A] = Field.nV[A_mask];
 }
}

return Field_Ap;
}

/*********************************************************************/

ChainPtr get_nodal_locality__MeshTools__(int I, Mesh FEM_Mesh){

  /* Define output */
  ChainPtr Nodes = NULL;

  /* Number of elements sourronding the node */
  int NumNeighbour = FEM_Mesh.NumNeighbour[I];
  
  /* Index of the elements sourronding the node */
  int * NodeNeighbour =
    set_to_memory__SetLib__(FEM_Mesh.NodeNeighbour[I],NumNeighbour);
  
  /* Table with the nodes of each element */
  ChainPtr * Table_ElemNodes = malloc(NumNeighbour*sizeof(ChainPtr));

  /* Fill each position of the table with a list of nodes in the element */
  for(int i = 0 ; i<NumNeighbour ; i++){
    Table_ElemNodes[i] = FEM_Mesh.Connectivity[NodeNeighbour[i]];
  }
  
  /* Free table with elements */
  free(NodeNeighbour);
  
  /* Get the union of this nodes */
  Nodes = union__SetLib__(Table_ElemNodes,NumNeighbour);
  
  /* Free table with the nodes of each elements */
  free(Table_ElemNodes);

  /* Return nodes close to the node I */
  return Nodes;
}

/*********************************************************************/

void get_nodal_connectivity__MeshTools__(Mesh FEM_Mesh){

  /* Variable declaration */
  int * Element_Connectivity;
  int NumNodesElem;
  
  /* 1º Start the search of neighbour for each node */
  for(int i = 0 ; i<FEM_Mesh.NumNodesMesh ; i++){
    /* 2º Loop over all the elements in the mesh */
    for(int j = 0 ; j<FEM_Mesh.NumElemMesh ; j++){
      NumNodesElem = FEM_Mesh.NumNodesElem[j];
      Element_Connectivity = set_to_memory__SetLib__(FEM_Mesh.Connectivity[j],NumNodesElem);
      /* 3º Loop over the all the node in an element */
      for(int k = 0 ; k<NumNodesElem ; k++){
	/* 4º If my node belong to the element */
	if(Element_Connectivity[k] == i){
	  /* 5º Introduce the element in the chain */
	  push__SetLib__(&FEM_Mesh.NodeNeighbour[i], j);
	  /* 6º Update the counter */
	  FEM_Mesh.NumNeighbour[i] += 1;	  
	}
      }
      /* Free memory */
      free(Element_Connectivity);
    }
  }
  
}

/*********************************************************************/

Matrix compute_N__MeshTools__(Element GP_Element,GaussPoint MPM_Mesh,Mesh FEM_Mesh) 
{ 
  int i_GP = GP_Element.i_GP;
  int GP_NumNodes = GP_Element.NumberNodes;
  int * GP_Connect = GP_Element.Connectivity;
  
  Matrix GP_ElemCoord; /* Coordinates of the nodes */
  int GP_I; /* Get the node for the GP */

  int Ndim = NumberDimensions;
  
  /* Gauss-Point properties */
  Matrix X_GP = /* Element coordinates of the Gauss-Point */
    memory_to_matrix__MatrixLib__(Ndim,1,NULL); 
  Matrix lp; /* Just for GIMP -> Particle voxel */
  Matrix Beta_GP =  /* Tunning parameter for LME */
    memory_to_matrix__MatrixLib__(Ndim,1,NULL);
  Matrix Delta_Xip; /* Just for GIMP -> Distance from GP to the nodes */
  Matrix lambda_GP = /* Just for LME/LME -> Lagrange multipliers */
    memory_to_matrix__MatrixLib__(Ndim,1,NULL);
  
  Matrix ShapeFunction_p; /* Matrix with the nodal shape functions */
  
  if(strcmp(ShapeFunctionGP,"MPMQ4") == 0){

    /* Fill the poligon */
    GP_ElemCoord = allocZ__MatrixLib__(GP_NumNodes,Ndim);
    for(int k = 0; k<GP_NumNodes ; k++){
      /* Get the node for the GP */
      GP_I = GP_Connect[k];
      for(int l = 0 ; l<Ndim ; l++){
	GP_ElemCoord.nM[k][l] =
	  FEM_Mesh.Coordinates.nM[GP_I][l];
      }
    }
    
    /* Get the element coordinates of the GP */
    X_GP.nV = MPM_Mesh.Phi.x_EC.nM[i_GP];

    /* Evaluate the shape function */
    ShapeFunction_p = N__Q4__(X_GP);
    
    /* Free memory */
    free__MatrixLib__(GP_ElemCoord);
  }
  else if(strcmp(ShapeFunctionGP,"uGIMP") == 0){
    /* Generate a matrix with the distances to the nodes */
    Delta_Xip = alloc__MatrixLib__(GP_NumNodes,Ndim);
    for(int k = 0 ; k<GP_NumNodes ; k++){
      /* Get the node for the GP */
      GP_I = GP_Connect[k];
      for(int l = 0 ; l<Ndim ; l++){
	Delta_Xip.nM[k][l] =
	  MPM_Mesh.Phi.x_GC.nM[i_GP][l]-
	  FEM_Mesh.Coordinates.nM[GP_I][l];
      }
    }
    
    /* Get the GP voxel */
    lp.nV = MPM_Mesh.lp.nM[i_GP];

    /* Evaluate the shape function */
    ShapeFunction_p = N__GIMP__(Delta_Xip,lp,FEM_Mesh.DeltaX);

    /* Free memory */
    free__MatrixLib__(Delta_Xip);
  }
  else if(strcmp(ShapeFunctionGP,"LME") == 0){
    /* Get the distance of the GP to the nodes */
    Delta_Xip = alloc__MatrixLib__(GP_NumNodes,Ndim);
    for(int k = 0 ; k<GP_NumNodes ; k++){
      GP_I = GP_Connect[k];
      for(int l = 0 ; l<Ndim ; l++){
	Delta_Xip.nM[k][l] =
	  MPM_Mesh.Phi.x_GC.nM[i_GP][l]-
	  FEM_Mesh.Coordinates.nM[GP_I][l];
      }
    }      
    /* Asign lambda and beta */
    lambda_GP.nV = MPM_Mesh.lambda.nM[i_GP];
    Beta_GP.nV = MPM_Mesh.Beta.nM[i_GP];
   
    /* Evaluate the shape function */
    ShapeFunction_p = p__LME__(Delta_Xip, lambda_GP,Beta_GP);
    
    /* Free memory */
    free__MatrixLib__(Delta_Xip);
  }

  else if(strcmp(ShapeFunctionGP,"aLME") == 0)
  {
    /* Get the distance of the GP to the nodes */
    Delta_Xip = alloc__MatrixLib__(GP_NumNodes,Ndim);
    for(int k = 0 ; k<GP_NumNodes ; k++)
    {
      GP_I = GP_Connect[k];
      for(int l = 0 ; l<Ndim ; l++)
      {
        Delta_Xip.nM[k][l] = MPM_Mesh.Phi.x_GC.nM[i_GP][l] - FEM_Mesh.Coordinates.nM[GP_I][l];
      }
    }  

    /* Asign lambda and beta */
    Tensor lambda_p = memory_to_tensor__TensorLib__(MPM_Mesh.lambda.nM[i_GP],1);
    Tensor Beta_p   = memory_to_tensor__TensorLib__(MPM_Mesh.Beta.nM[i_GP],2);
    
    /* Evaluate the shape function */
    ShapeFunction_p = p__aLME__(Delta_Xip, lambda_p, Beta_p);

    /* Free memory */
    free__MatrixLib__(Delta_Xip);
  }

  else{
    printf("%s : %s %s %s \n",
	   "Error in Get_Operator()",
	   "The shape-function ",
	   ShapeFunctionGP,
	   "is not implemented");      
    exit(EXIT_FAILURE);
  }

  return ShapeFunction_p;
}


/*********************************************************************/

Matrix compute_dN__MeshTools__(Element GP_Element,GaussPoint MPM_Mesh,
			      Mesh FEM_Mesh) 
{ 
  int i_GP = GP_Element.i_GP;
  int GP_NumNodes = GP_Element.NumberNodes;
  int * GP_Connect = GP_Element.Connectivity;

  int Ndim = NumberDimensions;
  
  Matrix GP_ElemCoord; /* Coordinates of the nodes */
  int GP_I; /* Get the node for the GP */
  
  /* Gauss-Point properties */
  Matrix X_GP = /* Element coordinates of the Gauss-Point */
    memory_to_matrix__MatrixLib__(Ndim,1,NULL); 
  Matrix lp; /* Just for GIMP -> Particle voxel */
  Matrix Beta_GP =  /* Tunning parameter for LME */
    memory_to_matrix__MatrixLib__(Ndim,1,NULL);
  Matrix Delta_Xip; /* Just for GIMP -> Distance from GP to the nodes */
  
  Matrix lambda_GP = /* Just for LME/LME -> Lagrange multipliers */
    memory_to_matrix__MatrixLib__(Ndim,1,NULL);
  
  Matrix ShapeFunction_p; /* Matrix with the nodal shape functions */
  Matrix Gradient_p; /* Matrix with the nodal derivatives */
  
  if(strcmp(ShapeFunctionGP,"MPMQ4") == 0){
    /* Fill the poligon */
    GP_ElemCoord = allocZ__MatrixLib__(GP_NumNodes,Ndim);
    
    for(int k = 0; k<GP_NumNodes ; k++){
      /* Get the node for the GP */
      GP_I = GP_Connect[k];
      for(int l = 0 ; l<Ndim ; l++){
	GP_ElemCoord.nM[k][l] =
	  FEM_Mesh.Coordinates.nM[GP_I][l];
      }
    }
    
    /* Get the element coordinates of the GP */
    X_GP.nV = MPM_Mesh.Phi.x_EC.nM[i_GP];

    /* Evaluate the shape function gradient */
    Gradient_p = dN__Q4__(X_GP,GP_ElemCoord);
    
    /* Free memory */
    free__MatrixLib__(GP_ElemCoord);
  }
  else if(strcmp(ShapeFunctionGP,"uGIMP") == 0){
    /* Generate a matrix with the distances to the nodes */
    Delta_Xip = alloc__MatrixLib__(GP_NumNodes,Ndim);
    for(int k = 0 ; k<GP_NumNodes ; k++){
      /* Get the node for the GP */
      GP_I = GP_Connect[k];
      for(int l = 0 ; l<Ndim ; l++){
	Delta_Xip.nM[k][l] =
	  MPM_Mesh.Phi.x_GC.nM[i_GP][l]-
	  FEM_Mesh.Coordinates.nM[GP_I][l];
      }
    }
    
    /* Get the GP voxel */
    lp.nV = MPM_Mesh.lp.nM[i_GP];
    
    /* Evaluate the shape function gradient */
    Gradient_p = dN__GIMP__(Delta_Xip,lp,FEM_Mesh.DeltaX);

    /* Free memory */
    free__MatrixLib__(Delta_Xip);
  }
  else if(strcmp(ShapeFunctionGP,"LME") == 0){
    /* Get the distance of the GP to the nodes */
    Delta_Xip = alloc__MatrixLib__(GP_NumNodes,Ndim);
    for(int k = 0 ; k<GP_NumNodes ; k++){
      GP_I = GP_Connect[k];
      for(int l = 0 ; l<Ndim ; l++){
	Delta_Xip.nM[k][l] =
	  MPM_Mesh.Phi.x_GC.nM[i_GP][l]-
	  FEM_Mesh.Coordinates.nM[GP_I][l];
      }
    }      
    /* Asign lambda and beta */
    lambda_GP.nV = MPM_Mesh.lambda.nM[i_GP];
    Beta_GP.nV = MPM_Mesh.Beta.nM[i_GP];
    
    /* Evaluate the shape function gradient */
    ShapeFunction_p = p__LME__(Delta_Xip, lambda_GP,Beta_GP);
    Gradient_p = dp__LME__(Delta_Xip, ShapeFunction_p);
   
    /* Free memory */
    free__MatrixLib__(ShapeFunction_p);
    free__MatrixLib__(Delta_Xip);
  }

  else if(strcmp(ShapeFunctionGP,"aLME") == 0)
  {
    /* Get the distance of the GP to the nodes */
    Delta_Xip = alloc__MatrixLib__(GP_NumNodes,Ndim);
    for(int k = 0 ; k<GP_NumNodes ; k++)
    {
      GP_I = GP_Connect[k];
      for(int l = 0 ; l<Ndim ; l++)
      {
        Delta_Xip.nM[k][l] = MPM_Mesh.Phi.x_GC.nM[i_GP][l] - FEM_Mesh.Coordinates.nM[GP_I][l];
      }
    }  

    /* Asign lambda and beta */
    Tensor lambda_p = memory_to_tensor__TensorLib__(MPM_Mesh.lambda.nM[i_GP],1);
    Tensor Beta_p   = memory_to_tensor__TensorLib__(MPM_Mesh.Beta.nM[i_GP],2);
    
    /* Evaluate the shape function gradient */
    ShapeFunction_p = p__aLME__(Delta_Xip, lambda_p, Beta_p);
    Gradient_p = dp__LME__(Delta_Xip, ShapeFunction_p);

    /* Free memory */
    free__MatrixLib__(ShapeFunction_p);
    free__MatrixLib__(Delta_Xip);
  }

  else{
    printf("%s : %s %s %s \n",
	   "Error in Get_Operator()",
	   "The shape-function ",
	   ShapeFunctionGP,
	   "is not implemented");      
    exit(EXIT_FAILURE);
  }

  return Gradient_p;
}

/*********************************************************************/

double mesh_size__MeshTools__(Mesh FEM_Mesh)
/*
  Function to get the minimum mesh size.
*/
{

  /* Auxiliar variables of the function */
  int NumElemMesh = FEM_Mesh.NumElemMesh;
  int NumNodesElem; /* Number of nodes of each element */
  int * Connectivity; /* Connectivity of the element */
  Matrix Poligon; /* Element Poligon */
  Matrix X_eval = allocZ__MatrixLib__(1,2); /* Where to evaluate the shape function */
  X_eval.nV[0] = 0;
  X_eval.nV[1] = 0;
  Matrix dNdx; /* Gradient of the shapefunction for each node */
  double MinElementSize_aux;
  double MinElementSize = 10e16;

  /* 1º Loop over the elements in the mesh */
  for(int i = 0 ; i<NumElemMesh ; i++){

    /* 2º Connectivity of the Poligon */
    NumNodesElem = FEM_Mesh.NumNodesElem[i];
    Connectivity = set_to_memory__SetLib__(FEM_Mesh.Connectivity[i],NumNodesElem);
    
    /* 4º Get the gradient of the element for each node */
    if((NumNodesElem == 3) &&
       (NumberDimensions == 2)){ /* Triangular element */
      /* The poligon is a triangle */
      Poligon = allocZ__MatrixLib__(3,2);
      /* Fill the triangle */
      for(int k = 0; k<3; k++){
	for(int l = 0 ; l<2 ; l++){
	  Poligon.nM[k][l] = FEM_Mesh.Coordinates.nM[Connectivity[k]][l];
	}
      }
      /* Get the gradient of the triangle */
      dNdx = dN__T3__(X_eval,Poligon);
      free__MatrixLib__(Poligon);
      
      /* Get the minimum minimum height of the triangle */
      for(int j = 0 ; j<3 ; j++){
	MinElementSize_aux =
	  1/pow(dNdx.nM[0][j]*dNdx.nM[0][j] +
		dNdx.nM[1][j]*dNdx.nM[1][j],0.5);
	MinElementSize = DMIN(MinElementSize,MinElementSize_aux);
      }
      /* Free memory */
      free__MatrixLib__(dNdx);
      
    }
    else if((NumNodesElem == 4) &&
	    (NumberDimensions == 2)){ /* Quadrilateral element */
      /* The poligon is a quadrilateral */
      Poligon = allocZ__MatrixLib__(4,2);

      /* Fill the poligon with vectors */
      for(int k = 0; k<3; k++){
	for(int l = 0 ; l<2 ; l++){
	  Poligon.nM[k][l] =
	    FEM_Mesh.Coordinates.nM[Connectivity[k+1]][l] -
	    FEM_Mesh.Coordinates.nM[Connectivity[k]][l];
	}
      }
      for(int l = 0 ; l<2 ; l++){
	Poligon.nM[3][l] = FEM_Mesh.Coordinates.nM[Connectivity[0]][l] -
	  FEM_Mesh.Coordinates.nM[Connectivity[3]][l];
      }
      
      /* Get the minimum minimum height of the triangle */
      for(int k = 0 ; k<4 ; k++){
	MinElementSize_aux = pow(Poligon.nM[k][0]*Poligon.nM[k][0] +
				 Poligon.nM[k][1]*Poligon.nM[k][1] , 0.5);
	MinElementSize = DMIN(MinElementSize,MinElementSize_aux);
      }

      /* Free memory */
      free__MatrixLib__(Poligon);

    }
    else{
      printf("%s : %s %i %s \n",
	     "Error in mesh_size__MeshTools__",
	     "Element with ",
	     NumNodesElem,
	     "nodes is not implemented !!!" );
      exit(EXIT_FAILURE);
    }

    /* Free memory */
    free(Connectivity);
    
  }

  /* Free memory */
  free__MatrixLib__(X_eval);

  return MinElementSize;

}

/*********************************************************************/

Matrix get_nodes_coordinates__MeshTools__(ChainPtr Element_p, Matrix Coordinates)
/*
  Get the matrix with the coordinates of an element
*/
{

  int Ndim = NumberDimensions;
  int NumVertex = lenght__SetLib__(Element_p);
  Matrix Element_Coordinates = allocZ__MatrixLib__(NumVertex,Ndim);
  ChainPtr Idx = NULL;
  int I_Idx = 0;

  Idx = Element_p;

  while(Idx != NULL){

    /* Fill the elelemtn coordinates */
    for(int l = 0 ; l<Ndim ; l++){
      Element_Coordinates.nM[I_Idx][l] = Coordinates.nM[Idx->I][l];
    }
    
    /* Cycle */
    Idx = Idx->next;
    I_Idx++;
  }
  
  return Element_Coordinates;
}

/*********************************************************************/


Matrix compute_distance__MeshTools__(ChainPtr Set, Matrix X0, Matrix Coordinates)
{
  int Ndim = NumberDimensions;
  int SizeSet = lenght__SetLib__(Set);
  int I = 0;

  /* Allocate output */
  Matrix Set_Coordinates = allocZ__MatrixLib__(SizeSet,Ndim);

  /* Loop in the set */
  ChainPtr Aux_Set = Set;
  while (Aux_Set != NULL){ 
    /* Get coordinates local coodinates of each node in the set */
    for(int i = 0 ; i<Ndim ; i++){
      Set_Coordinates.nM[I][i] =  X0.nV[i] - Coordinates.nM[Aux_Set->I][i];
    }
    /* Update index */
    I++;	
    Aux_Set = Aux_Set->next; 
  }
  
  return Set_Coordinates;
}

/*********************************************************************/

Matrix get_set_field_old__MeshTools__(Matrix Nodal_Field, Element GP_Element){

  int * Element_Connectivity = GP_Element.Connectivity;
  int NumNodes = GP_Element.NumberNodes;
  Matrix Element_Field;
  int Ndim = NumberDimensions;
  int Ie;

  /* Allocate a matrix to store the nodal quatities in the element */
  Element_Field = alloc__MatrixLib__(NumNodes,Ndim);

  /* Loop over the nodes of the element */
  for(int I = 0; I<NumNodes; I++){
    /* Get the node of the element */
    Ie = Element_Connectivity[I];
    /* Fill each dimension of the nodal quantitie */
    for(int i = 0 ; i<Ndim ; i++){
      Element_Field.nM[I][i] = Nodal_Field.nM[Ie][i];
    }
  }
  
  return Element_Field;
}
 
/*********************************************************************/

int get_closest_node__MeshTools__(Matrix X_p, ChainPtr Nodes, Matrix Coordinates)
/*
  Ordenate recursively and array with distances and get a chain 
  with the positions in orden
*/
{
  /* Number of dimensions */
  int Ndim = NumberDimensions;
  /* Index of the current node */
  int I;
  /* Corrdinates of the current node */
  Matrix X_I;
  /* Distance from the particle to the node */
  double Distance_I;
  /* Interator pointer */
  ChainPtr Node_I = NULL;
  /* Get the value of the maximum distance from the particle to the pointer */
  double DistMin;
  int I_DistMin;

  /* Initialize interator with the first node */
  Node_I = Nodes;

  /* Get the index of the first node */
  I = Node_I->I;
  
  /* Get the coordinates of the first node */
  X_I = memory_to_matrix__MatrixLib__(Ndim,1,Coordinates.nM[I]);

  /* Get the distance from the node to the particle */
  Distance_I = point_distance__MatrixLib__(X_p, X_I);
      
  /* Get the distance from the node to the particle */
  DistMin = Distance_I;
  I_DistMin = Node_I->I;
      
  /* Search in the reamaining nodes */      
  Node_I = Node_I->next;
  
  while(Node_I != NULL){

    /* Get the index of the node */
    I = Node_I->I;

    /* Get the coordinates of the node */
    X_I = memory_to_matrix__MatrixLib__(Ndim,1,Coordinates.nM[I]);

    /* Get the distance from the node to the particle */
    Distance_I = point_distance__MatrixLib__(X_p, X_I);
      
    /* Get the max distance of the matrix */
    if(Distance_I < DistMin){
      DistMin = Distance_I;
      I_DistMin = Node_I->I;
    }
      
    /* Continue iterating */      
    Node_I = Node_I->next;
      
  }
   
  return I_DistMin;
}

/*********************************************************************/

bool inout_convex_set__MeshTools__(Matrix X_p, ChainPtr Elem_p, Matrix Coordinates)
{

  bool Is_In_Element = false;
  Matrix Element_Coordinates;
  
  Element_Coordinates = get_nodes_coordinates__MeshTools__(Elem_p, Coordinates);
  
  if (inout__MatrixLib__(X_p, Element_Coordinates) == 1){
    Is_In_Element = true;
  }

  free__MatrixLib__(Element_Coordinates);

  return Is_In_Element;
}

/*********************************************************************/

double interpolate_scalar__MeshTools__(Matrix A, Matrix N_p)
{
  int Nnodes_p = A.N_rows;
  double A_p = 0;

  for(int I = 0 ; I<Nnodes_p ; I++)
  {
    A_p += A.nV[I]*N_p.nV[I];
  }

  return A_p; 
}

/*********************************************************************/

