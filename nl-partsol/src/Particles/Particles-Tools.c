#include "nl-partsol.h"

#define MAXVAL(A, B) ((A) > (B) ? (A) : (B))
#define MINVAL(A, B) ((A) < (B) ? (A) : (B))

/*********************************************************************/

void initial_position__Particles__(Matrix X_p, Mesh FEM_Mesh, int GPxElement)
/*
 *
 */
{

  if (strcmp(FEM_Mesh.TypeElem, "Triangle") == 0) {
    element_to_particles__T3__(X_p, FEM_Mesh, GPxElement);
  } else if (strcmp(FEM_Mesh.TypeElem, "Quadrilateral") == 0) {
    element_to_particles__Q4__(X_p, FEM_Mesh, GPxElement);
  } else if (strcmp(FEM_Mesh.TypeElem, "Tetrahedra") == 0) {
    element_to_particles__T4__(X_p, FEM_Mesh, GPxElement);
  } else if (strcmp(FEM_Mesh.TypeElem, "Hexahedra") == 0) {
    element_to_particles__H8__(X_p, FEM_Mesh, GPxElement);
  } else {
    fprintf(stderr, "%s : %s \n", "Error in initial_position__Particles__()",
            "Wrong type of shape function");
    exit(EXIT_FAILURE);
  }
}

/*********************************************************************/

int search_particle_in_surrounding_elements__Particles__(int p, Matrix X_p,
                                                         ChainPtr ListElement,
                                                         Mesh FEM_Mesh)
/*

*/
{
  ChainPtr Ixd = NULL;
  int I_element = -999;
  int Nn; /* Number of nodes of the element */
  ChainPtr Nodes;
  Matrix Element_Coordinates;

  Ixd = ListElement;

  while (Ixd != NULL) {

    Nn = FEM_Mesh.NumNodesElem[Ixd->Idx];
    Nodes = FEM_Mesh.Connectivity[Ixd->Idx];

    Element_Coordinates =
        get_nodes_coordinates__MeshTools__(Nodes, FEM_Mesh.Coordinates);

    /* Check if the particle is in the element */
    if (FEM_Mesh.In_Out_Element(X_p, Element_Coordinates)) {
      free__MatrixLib__(Element_Coordinates);
      I_element = Ixd->Idx;
      break;
    }

    /* Cycle */
    free__MatrixLib__(Element_Coordinates);
    Ixd = Ixd->next;
  }

  return I_element;
}

/*********************************************************************/

Element nodal_set__Particles__(int i_GP, ChainPtr ListNodes, int NumNodes) {

  /* Define new element */
  Element GP_Element;

  /* Fill element */
  GP_Element.i_GP = i_GP;
  GP_Element.NumberNodes = NumNodes;
  GP_Element.Connectivity = set_to_memory__SetLib__(ListNodes, NumNodes);

  return GP_Element;
}

/*********************************************************************/

void asign_to_nodes__Particles__(int p, int E_p, int I0, ChainPtr ListNodes_p,
                                 Mesh FEM_Mesh) {

  /*!
   * Assign particle to an element of the background mesh
   * */
  //  push__SetLib__(&FEM_Mesh.List_Particles_Element[E_p],p);
  //  FEM_Mesh.Num_Particles_Element[E_p] += 1;

  /*!
   * Assign particle to a node of the background mesh
   * */
  push__SetLib__(&FEM_Mesh.List_Particles_Node[I0], p);
  FEM_Mesh.Num_Particles_Node[I0] += 1;
}

/*********************************************************************/
