#include <stdlib.h>
#include "gidpost.h"

/*
  #define ASCII_FORMAT
  #define DO_FLUSH
  #define EARLY_TERMINATE
*/

#define USE_GP_IN_GROUP

#if defined(USE_GP_IN_GROUP)
#define LPL 3
#define LOC elems
#else
#define LPL 1
#define LOC nodes
#endif

typedef struct {
  int id;
  float x, y, z;
} SCoord;

typedef struct {
  int id;
  int n1, n2, n3;
} SElem;

SCoord nodes[9];
SElem elems[8];

void GenMesh()
{
  int i, j, idx;
  SCoord * node;
  SElem * elm1, * elm2;

  idx = 1;
  for ( i = 0; i < 3; i++ ) {
    for ( j = 0; j < 3; j++, idx++ ) {
      node = nodes + idx -1;
      node->id = idx;
      node->x = (float)(i);
      node->y = (float)(j);
      node->z = 0;
    }
  }
  idx = 0;
  for ( i = 0; i < 2; i++ ) {
    for ( j = 0; j < 2; j++, idx+=2 ) {
      elm1 = elems + idx;
      elm2 = elm1 + 1;
      elm1->id = idx+1;
      elm1->n1 = i*3 + j + 1;
      elm1->n3 = i*3 + j + 1 + 3 ;
      elm1->n2 = elm1->n3 + 1;
      elm2->id = idx+2;
      elm2->n1 = elm1->n1;
      elm2->n2 = elm1->n1 + 1;
      elm2->n3 = elm1->n2;
    }
  }
}

float Random()
{
  return rand()/(float)(RAND_MAX);
}

int main()
{
  int i, j;
  SCoord * node;
  SElem * elm;
  int elmi[4];
  GiD_FILE fdm, fdr;

  GenMesh();

  GiD_PostInit();

#if defined(ASCII_FORMAT)
  fdm = GiD_fOpenPostMeshFile( "test_fd.post.msh", GiD_PostAscii);
  fdr = GiD_fOpenPostResultFile( "test_fd.post.res", GiD_PostAsciiZipped);
#else
  fdm = fdr = GiD_fOpenPostResultFile( "test_fd.post.res", GiD_PostBinary);
#endif

  /* write mesh info */
  GiD_fBeginMeshColor(fdm, "TestMsh", GiD_2D, GiD_Triangle, 3, 0, 0.99, 0);
  /* coordinates */
  GiD_fBeginCoordinates(fdm);
  for ( i = 0; i < 9; i++ ) {
    node = nodes + i;
    GiD_fWriteCoordinates(fdm, node->id, node->x, node->y, node->z ); 
  }
  GiD_fEndCoordinates(fdm);
  /* elements */
  GiD_fBeginElements(fdm);
  for ( i = 0; i < 8; i++ ) {
    elm = elems + i;
    elmi[0] = elm->n1;
    elmi[1] = elm->n2;
    elmi[2] = elm->n3;
    elmi[3] = 2;
    GiD_fWriteElementMat(fdm, elm->id, elmi);
  }
  GiD_fEndElements(fdm);
  GiD_fEndMesh(fdm);

  /* now results info */

  /*
    #define SIMPLE_RESULT
  */
#ifdef SIMPLE_RESULT
  GiD_fBeginResult(fdr,
                   "EscalarNodos", "Analysis", 1.0, GiD_Scalar, GiD_OnNodes,
		  NULL, NULL, 0, NULL);
  for ( i = 0; i < 9; i++ ) { 
    GiD_fWriteScalar(fdr, nodes[i].id, Random());
  }
  GiD_fEndResult(fdr);
#else
  /* write the gauss points */
  GiD_fBeginGaussPoint(fdr, "element_gp", GiD_Triangle, NULL, 3, 0, 1);  
  GiD_fEndGaussPoint(fdr);
  
  /* scalar result over nodes */
  /*
    #if defined(ASCII_FORMAT)  
  */
#if defined(USE_GP_IN_GROUP)
  GiD_fBeginResultGroup(fdr, "Analysis", 1.0, GiD_OnGaussPoints,"element_gp" );
#else
  GiD_fBeginResultGroup(fdr, "Analysis", 1.0, GiD_OnNodes, NULL);
#endif
  GiD_fResultDescription(fdr, "EscalarNodosG", GiD_Scalar);
  GiD_fResultDescriptionDim(fdr, "VectorNodosG", GiD_Vector, 4);
  GiD_fResultDescription(fdr, "MatrixG", GiD_Matrix);
  GiD_fResultDescription(fdr, "Local AxesG", GiD_LocalAxes);
  for ( i = 0; i < sizeof(LOC)/sizeof(LOC[0]); i++ ) {
    for (j = 0; j < LPL; j++) {
      GiD_fWriteScalar(fdr, LOC[i].id, Random());
      GiD_fWriteVectorModule(fdr, LOC[i].id, Random(), Random(), Random(),-1);
      GiD_fWrite3DMatrix(fdr, LOC[i].id, Random(), Random(), Random(),
                         Random(), Random(), Random());
      GiD_fWriteLocalAxes(fdr, LOC[i].id, Random(), Random(), Random());
    }
  }
  GiD_fEndResult(fdr);
  
  /* #else */

  /* scalar results over gauss points: Triangle + 3 GP */
  
  GiD_fBeginResult(fdr,
                   "SobreGPTriangle3", "Analysis", 1.0, GiD_Scalar,
                   GiD_OnGaussPoints, "element_gp", NULL, 0, NULL);
  for ( i = 0; i < 8; i++ ) {
    for ( j = 0; j < 3; j++ )
      GiD_fWriteScalar(fdr, i+1, Random());
  }
  GiD_fEndResult(fdr);
  GiD_fBeginResult(fdr,
                   "EscalarNodos", "Analysis", 1.0, GiD_Scalar, GiD_OnNodes,
		  NULL, NULL, 0, NULL);
  for ( i = 0; i < 9; i++ ) { 
    GiD_fWriteScalar(fdr, nodes[i].id, Random());
  }
  GiD_fEndResult(fdr);

  /* vector result over nodes */
  
  /*
  GiD_fBeginResult("VectorNodos", "Analysis", 1.0, GiD_Vector, GiD_OnNodes,
  NULL, NULL, 0, NULL);*/
  GiD_fBeginResultHeader(fdr,
                         "VectorNodos", "Analysis", 1.0, GiD_Vector, GiD_OnNodes, NULL);
  for ( i = 0; i < 9; i++ ) { 
    GiD_fWriteVector(fdr,
                     nodes[i].id, Random(), Random(), Random());
  }
  GiD_fEndResult(fdr);

#ifndef DISABLE_CODE

#ifdef DO_FLUSH
  GiD_fFlushPostFile(fdr);
#endif 
  /* matrix result */

  GiD_fBeginResult(fdr,
                   "Matrix", "Analysis", 1.0, GiD_Matrix, GiD_OnNodes,
                   NULL, NULL, 0, NULL);
  for ( i = 0; i < 9; i++ ) {
    GiD_fWrite3DMatrix(fdr, i+1, Random(), Random(), Random(),
                       Random(), Random(), Random());
  }
  GiD_fEndResult(fdr);  

#ifdef EARLY_TERMINATE
  abort();
#endif
  /* local axes result */

  GiD_fBeginResult(fdr,
                   "Local Axes", "Analysis", 1.0, GiD_LocalAxes, GiD_OnNodes,
                   NULL, NULL, 0, NULL);
  for ( i = 0; i < 9; i++ ) {
    GiD_fWriteLocalAxes(fdr, i+1, Random(), Random(), Random());
  }
  GiD_fEndResult(fdr);  
  /* #endif */ /* ASCII_FORMAT */
#endif /* DISABLE_CODE */
#endif /* SIMPLE_RESULT */
  
#ifdef ASCII_FORMAT
  GiD_fClosePostMeshFile(fdm);
#endif
  GiD_fClosePostResultFile(fdr);

  GiD_PostDone();

  return 0;
}
