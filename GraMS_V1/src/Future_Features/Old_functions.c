



/* void Get_Deformation_Gradient(Element * Elem, /\* Gauss point *\/ */
/* 			      GaussPoint * GP) /\* Element *\/ */
/* { */
  
/*   /\* Deformation reference gradient F_ref *\/ */
/*   Get_RefDeformation_Gradient(GP,Elem); */

/*   /\* Get the inverse of the deformation gradient *\/ */
/*   Tensor F_ref_m1; */
/*   F_ref_m1 = Get_Inverse(GP[0].F_ref); */

/*   /\* Iterate over the nodes in the Element to get B *\/ */
/*   for(int i = 0 ; i<Element[0].NumberNodes ; i++){ */

    
    
/*   }/\* Nodes loop *\/ */

/*   /\* Once we have calculated the deformation gradient */
/*      matrix we dont need anymore the inverse of the */
/*      reference deformation gradient so we free memory *\/ */
/*   free(F_ref_m1.n); */
  
/* } */

/*********************************************************************/


/* Matrix Get_Stiffness_Matrix(Element * Elem) /\* Element to analyze *\/ */
/* /\*  */
/*    Calcule the element-stiffness matrix K : */
/*    Inputs : Constitutive matrix */
/*    Outputs : Element */
/* *\/ */
/* { */

/*   Matrix dNdX; */
/*   Matrix dNdX_T; */
/*   Matrix K; */
/*   Matrix K_GP; */
/*   Matrix D; */
/*   Matrix J; */

/*   /\* Auxiliar pointer *\/ */
/*   GaussPoint * GP_ei; */
  
/*   /\* Allocate and initialize K to zero *\/ */
/*   K = MatAllocZ(Elem->NumberNodes * Elem->NumberDOF, */
/* 		 Elem->NumberNodes * Elem->NumberDOF); */
    
/*   /\* Iterate over the Gauss Points *\/ */
/*   for(int i_GP = 0 ; i_GP<Elem->NumberGP ; i_GP++){ */

/*     /\* Get the derivative Matrix and it transposse evaluated in the gauss point *\/ */

/*     /\* Pointer to the GP in the element *\/ */
/*     GP_ei = &Elem->GP_e[i_GP]; */
    
/*     dNdX = Get_dNdx_matrix(Elem,GP_ei); */
/*     dNdX_T= Transpose_Mat(Get_dNdx_matrix(Elem,GP_ei)); */

/*     /\* Get the jacobian of the transformation and store it as a scalar *\/ */
/*     J.n = Get_Determinant(GP_ei->F_ref); */
/*     J.nM = NULL; */
/*     J.nV = NULL; */
/*     J.N_cols = 1; */
/*     J.N_rows = 1; */

/*     /\* Multiply the constitutive matrix by the derivative matrix *\/ */
/*     D = CopyMat(GP_ei->D);     */
/*     K_GP = Scalar_prod(D,dNdX); */

/*     /\* Multiply the result by the derivative matrix transpose *\/ */
/*     K_GP = Scalar_prod(dNdX_T,K_GP); */
    
/*     /\* Multiply by the jacobian*\/ */
/*     K_GP = Scalar_prod(K_GP,J); */

/*     /\* Acumulate for the integral rule *\/ */
/*     K = Add_Mat(K,K_GP); */
    
/*   } */
  
/*   return K; */
/* } */

/*********************************************************************/

/* Matrix Get_Geom_Mass_Matrix(GaussPoint GP_Mesh, */
/* 			    Mesh ElementMesh)  */
/* /\* */
/*    Calcule the mass matrix M : */
/*    Inputs : Element */
/*    Outputs : Mass matrix (M) */
/* *\/ */
/* { */
/*   /\* Mass matrix *\/ */
/*   Matrix M; */
/*   /\* Index of the element where it is the GP *\/ */
/*   int i_GP_Elem; */
  
/*   /\* Coordenates of the element nodes *\/ */
/*   Matrix X_ElemNod; */
/*   /\* Pointer to the Connectivity of the element *\/ */
/*   int * Id_ElemNod; */
/*   /\* Derivative Matrix *\/ */
/*   Matrix dNdX_Ref_GP; */

/*   /\* Jacobian *\/ */
/*   Matrix J_GP; */
  
/*   /\* Allocate and initialize M to zero *\/ */
/*   M = MatAllocZ(ElementMesh.NumNodesMesh, */
/* 		ElementMesh.NumNodesMesh); */

/*   /\* Allocate the mesh with the nodal coordinates *\/ */
/*   X_ElemNod = MatAlloc(ElementMesh.NumNodesElem, */
/* 		       ElementMesh.Dimension); */

/*   /\* Iterate over the Gauss Points *\/ */
/*   for(int i_GP = 0 ; i_GP<GP_Mesh.NumGP ; i_GP++){ */
    
/*     i_GP_Elem = GP_Mesh.Element_id[i_GP]; */

/*     /\* Calcule the Jacobian of the element evaluated in the GP : *\/ */
/*     /\* Get the element connectivity and take in to account that  */
/*        in the C programming languaje the index starts in 0 *\/ */
/*     Id_ElemNod = ElementMesh.Connectivity[i_GP_Elem]; */
/*     /\* Get the coordinates of the nodes *\/ */
/*     X_ElemNod.nV[0] = ElementMesh.Coordinates.nM[Id_ElemNod[0] - 1][0]; */
/*     X_ElemNod.nV[1] = ElementMesh.Coordinates.nM[Id_ElemNod[1] - 1][0]; */
/*     dNdX_Ref_GP = ElementMesh.dNdX_ref(X_ElemNod); */
/*     J_GP = Get_RefDeformation_Gradient(X_ElemNod,dNdX_Ref_GP); */

/*     M.nM[Id_ElemNod[0]-1][Id_ElemNod[0]-1] += (double)1/3 * J_GP.n; */
/*     M.nM[Id_ElemNod[0]-1][Id_ElemNod[1]-1] += (double)1/6 * J_GP.n; */
/*     M.nM[Id_ElemNod[1]-1][Id_ElemNod[0]-1] += (double)1/6 * J_GP.n; */
/*     M.nM[Id_ElemNod[1]-1][Id_ElemNod[1]-1] += (double)1/3 * J_GP.n; */
    
/*   } */

/*   return M; */


/* } */


/*********************************************************************/




/* /\*********************************************************************\/ */

/* Matrix GetMassMatrix_L(Mesh ElementMesh, */
/* 		       GaussPoint GP_Mesh){ */
  
/*   Matrix M_l = MatAllocZ(ElementMesh.NumNodesMesh,1); */
/*   int Elem_GP_i; /\* Element where the GP is placed *\/ */
/*   double M_GP; /\* Mass of the Gauss-Point *\/ */
/*   int * Nodes_Elem_GP_i; /\* Nodes of the element where the GP is placed *\/ */
/*   Matrix GP_i_XE; /\* Element coordinates of the Gauss Points *\/ */
/*   Matrix N_ref_XG; /\* Element function evaluated in the Gauss-Points *\/ */

/*   for(int i=0 ; i<GP_Mesh.NumGP ; i++){ */
     
/*     /\* 4º Index of the element where the G-P is located *\/ */
/*     Elem_GP_i = GP_Mesh.Element_id[i]; */
    
/*     /\* 5º List of nodes of the element *\/ */
/*     Nodes_Elem_GP_i = ElementMesh.Connectivity[Elem_GP_i]; */
    
/*     /\* 6º Element coordinates of the G-P *\/ */
/*     GP_i_XE.n = GP_Mesh.Phi.x_EC.nV[i]; */
    
/*     /\* 7º Evaluate the shape functions of the element in the  */
/*        coordinates of the G-P *\/ */
/*     N_ref_XG = ElementMesh.N_ref(GP_i_XE); */

/*     /\* 8º Get the mass of the material point *\/ */
/*     M_GP = GP_Mesh.Phi.mass.nV[i]; */
    
/*     /\* 8º Iterate over the nodes of the element *\/ */
/*     for(int j=0 ; j<ElementMesh.NumNodesElem ; j++){ */
/*       M_l.nV[Nodes_Elem_GP_i[j]-1] += M_GP*N_ref_XG.nV[j]; */
/*     } */
    
/*   } */

/*   return M_l; */
/* } */


/* /\*********************************************************************\/ */

/* void GaussPointsToMesh(Mesh ElementMesh, */
/* 		       GaussPoint GP_Mesh, */
/* 		       Matrix Phi_n_GP, */
/* 		       Matrix Phi_n_Nod, */
/* 		       Matrix M_l){ */

/*   printf(" * Transfer the Gauss-Points information to the mesh \n"); */
/*   /\************* Transfer the GP information to the mesh ****************\/ */
/*   /\**********************************************************************\/ */
/*   /\*  [N0]-------(e0)-------[N1]-------(e1)-------[N2]-- ;  t = n       *\/ */
/*   /\* Phi_n_GP0 <-- Phi_e0 --> Phi_n_GP1 <-- Phi_e1 --> Phi_n_GP         *\/ */
/*   /\**********************************************************************\/ */
  
/*   /\* 0º Define auxiliar variable to the transference from */
/*      the nodes to the Gauss Points *\/ */
/*   int NumDOF; /\* Degree of freedom of the GP *\/ */
/*   int Elem_GP_i; /\* Element where the GP is placed *\/ */
/*   double M_GP; /\* Mass of the Gauss-Point *\/ */
/*   int * Nodes_Elem_GP_i; /\* Nodes of the element where the GP is placed *\/ */
/*   Matrix GP_i_XE; /\* Element coordinates of the Gauss Points *\/ */
/*   Matrix N_ref_XG; /\* Element function evaluated in the Gauss-Points *\/ */
   
/*   /\* 1º Set the number of degree of freedom of the problem *\/ */
/*   NumDOF = Phi_n_GP.N_rows; */

/*   /\* 2º Check if it is necessary to do a resize of Phi_n_Nod  *\/ */
/*   if(ElementMesh.NumNodesMesh != Phi_n_Nod.N_cols){ */
/*     if(NumDOF >= 2){ /\* 2aº Vectorial *\/ */
/*       FreeMat(Phi_n_Nod); /\* Free previous data field *\/ */
/*       Phi_n_Nod = MatAllocZ(NumDOF,ElementMesh.NumNodesMesh); /\* Resize *\/ */
/*     } */
/*     else{/\* 2bº Scalar *\/ */
/*       FreeMat(Phi_n_Nod); /\* Free previous data field *\/ */
/*       Phi_n_Nod = MatAllocZ(NumDOF,ElementMesh.NumNodesMesh); /\* Resize *\/ */
/*     } */
/*   } */
/*   else{ /\* 3º If is not necessary to resize, initialize it to zero *\/ */
/*     for(int i = 0 ; i<ElementMesh.NumNodesMesh ; i++){ */
/*       if(NumDOF >= 2){ /\* 3aº Vectorial *\/ */
/* 	for(int j = 0 ; j<NumDOF ; j++){ */
/* 	  Phi_n_Nod.nM[j][i] = 0;  */
/* 	} */
/*       } */
/*       else{ /\* 3bº Scalar *\/ */
/* 	Phi_n_Nod.nV[i] = 0;  */
/*       }       */
/*     } */
/*   } */
  
/*   for(int i=0 ; i<GP_Mesh.NumGP ; i++){ */
     
/*     /\* 4º Index of the element where the G-P is located *\/ */
/*     Elem_GP_i = GP_Mesh.Element_id[i]; */
    
/*     /\* 5º List of nodes of the element *\/ */
/*     Nodes_Elem_GP_i = ElementMesh.Connectivity[Elem_GP_i]; */
    
/*     /\* 6º Element coordinates of the G-P *\/ */
/*     GP_i_XE.n = GP_Mesh.Phi.x_EC.nV[i]; */
    
/*     /\* 7º Evaluate the shape functions of the element in the  */
/*        coordinates of the G-P *\/ */
/*     N_ref_XG = ElementMesh.N_ref(GP_i_XE); */

/*     /\* 8º Get the mass of the material point *\/ */
/*     M_GP = GP_Mesh.Phi.mass.nV[i]; */
    
/*     /\* 8º Iterate over the nodes of the element *\/ */
/*     for(int j=0 ; j<ElementMesh.NumNodesElem ; j++){ */
/*       /\* 9º Iterate over the fields *\/ */
/*       for(int k=0 ; k<NumDOF ; k++){ */
/* 	/\* 10º Update the field value *\/; */
/* 	Phi_n_Nod.nM[k][Nodes_Elem_GP_i[j]-1] += Phi_n_GP.nM[k][i]*M_GP*N_ref_XG.nV[j]; */
/*       } */
/*     } */
    
/*   } */

/*   /\* 11º Fix the nodal value with the lumped mass matrix *\/ */
/*   for(int i = 0 ; i<ElementMesh.NumNodesMesh ; i++){ */
/*     for(int j = 0 ; j<NumDOF ; j++){ */
/*       Phi_n_Nod.nM[j][i] /= M_l.nV[i];  */
/*     } */
/*   } */


/* } */

/* /\*********************************************************************\/ */

/* void MeshToGaussPoints(Mesh ElementMesh, */
/* 		       GaussPoint GP_Mesh, */
/* 		       Matrix Phi_n_Nod, */
/* 		       Matrix Phi_n_GP, */
/* 		       Matrix M_l){ */

/*   printf(" * Transfer the nodal values to the Gauss-Points \n"); */
/*   /\************ Transfer the nodal values of U_n1 to the G-P *************\/ */
/*   /\***********************************************************************\/ */
/*   /\*   Phi_n ---> Phi_e <--- Phi_n ---> Phi_e <--- Phi_n                 *\/ */
/*   /\*   [N0]-------(e0)-------[N1]-------(e1)-------[N2]-- ;  t = n + 1   *\/ */
/*   /\***********************************************************************\/ */

/*   /\* 0º Define auxiliar variable to the transference from */
/*      the Gauss Points to the nodes *\/ */
/*   int NumDOF; /\* Degree of freedom of the GP *\/ */
/*   int Elem_GP_i; /\* Element where the GP is placed *\/ */
/*   int * Nodes_Elem_GP_i; /\* Nodes of the element where the GP is placed *\/ */
/*   Matrix GP_i_XE; /\* Element coordinates of the Gauss Points *\/ */
/*   Matrix N_ref_XG; /\* Element function evaluated in the Gauss-Points *\/ */

/*   NumDOF = Phi_n_Nod.N_rows; */
  
/*   /\* 1º Iterate over the G-P *\/ */
/*   for(int i=0 ; i<GP_Mesh.NumGP ; i++){ */
    
/*     /\* 2º Index of the element where the G-P is located *\/ */
/*     Elem_GP_i = GP_Mesh.Element_id[i]; */
    
/*     /\* 3º List of nodes of the element *\/ */
/*     Nodes_Elem_GP_i = ElementMesh.Connectivity[Elem_GP_i]; */
    
/*     /\* 4º Element coordinates of the G-P *\/ */
/*     GP_i_XE.n = GP_Mesh.Phi.x_EC.nV[i]; */

/*     /\* 5º Evaluate the shape functions of the element in the coordinates of the G-P *\/ */
/*     N_ref_XG = ElementMesh.N_ref(GP_i_XE); */

/*     /\* 6º Set to zero the field value *\/ */
/*     for(int k=0 ; k<NumDOF ; k++){ */
/*       /\* 7º Update the field value *\/ */
/*       Phi_n_GP.nM[k][i] = 0; */
/*     } */

/*     /\* 8º Iterate over the nodes of the element *\/ */
/*     for(int j=0 ; j<ElementMesh.NumNodesElem ; j++){ */
/*       /\* 9º Iterate over the fields *\/ */
/*       for(int k=0 ; k<NumDOF ; k++){ */
/*   	/\* 10º Update the field value *\/ */
/*   	Phi_n_GP.nM[k][i] += N_ref_XG.nV[j]*Phi_n_Nod.nM[k][Nodes_Elem_GP_i[j]-1]; */
/*       } */
/*     } */
/*   } */
/* } */

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


/*********************************************************************/

/* Mesh RectangularMesh(double X0, double Y0, */
/* 		     double X1, double Y1, */
/* 		     double Dx, double Dy, */
/* 		     char * ElementType){ */
/*   Mesh OutMesh; */
/*   int node_i = 0; */

/*   if( strcmp(ElementType,"Linear") == 0 ){ */
/*     strcpy(OutMesh.TypeElem,ElementType); */
/*     OutMesh.NumNodesElem = 4; */
/*     OutMesh.Dimension = 2; */
/*     OutMesh.N_ref = Q4; */
/*     OutMesh.dNdX_ref = dQ4; */
/*     OutMesh.NumNodesMesh = (int)(((X1-X0)/Dx + 1) * */
/* 				 ((Y1-Y0)/Dy + 1)); */
/*     OutMesh.NumElemMesh = (int)(((X1-X0)/Dx) * */
/* 				((Y1-Y0)/Dy)); */
/*     OutMesh.NumNodesBound = (int)(((X1-X0)/Dx - 1)*2 + */
/* 				  ((Y1-Y0)/Dy - 1)*2) + 4; */

/*     /\* Allocate matrix with the coordinates of the nodes *\/ */
/*     OutMesh.Coordinates = MatAlloc(OutMesh.NumNodesMesh,2); */
/*     OutMesh.Connectivity = (int **)Allocate_Matrix(OutMesh.NumElemMesh,4,sizeof(int)); */
/*     OutMesh.ActiveNode = (int *)Allocate_ArrayZ(OutMesh.NumNodesMesh,sizeof(int)); */

/*     /\* Fill matrix with coordinates *\/ */
/*     for(int i = 0 ; i<(int)((X1-X0)/Dx + 1) ; i++){ */
/*       for(int j = 0 ; j<(int)((Y1-Y0)/Dy + 1) ; j++){ */
/* 	OutMesh.Coordinates.nM[node_i][0] = Dx*i; */
/* 	OutMesh.Coordinates.nM[node_i][1] = Dy*j; */
/* 	node_i++; */
/*       } */
/*     } */
    
/*     /\* Fill conectivity mesh *\/ */
        
/*   } */
 
  
/*   return OutMesh; */
/* } */
