#define USE_CONST
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include "gidpost.h"

#ifdef WIN32
#define strcasecmp  _stricmp
#define strdup  _strdup
#endif

/* #define DO_FLUSH */

#define USE_GP_IN_GROUP

#if defined(USE_GP_IN_GROUP)
#define LPL     3
#define LOC     elems
#define LAST_ID last_elem_id
#else
#define LPL     1
#define LOC     nodes
#define LAST_ID last_node_id
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
  SElem * elem1, * elem2;

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
      elem1 = elems + idx;
      elem2 = elem1 + 1;
      elem1->id = idx+1;
      elem1->n1 = i*3 + j + 1;
      elem1->n3 = i*3 + j + 1 + 3 ;
      elem1->n2 = elem1->n3 + 1;
      elem2->id = idx+2;
      elem2->n1 = elem1->n1;
      elem2->n2 = elem1->n1 + 1;
      elem2->n3 = elem1->n2;
    }
  }
}

float Random()
{
  return rand()/(float)(RAND_MAX);
}

void print_help_and_exit( const char *prog_name) {
  printf( "Usage: %s [ -h] [ -abort] [ -f ascii|bin|hdf5] filename_prefix\n", prog_name);
  printf( "    will create filename_prefix.post.{ msh,res | bin | h5} depending on the choosen format.\n");
  printf( "    Use: %s -- filename_prefix          if filename_prefix begins with a '-'\n", prog_name);
  printf( "            -abort                does not close the file properly, after last result it aborts\n" );
  exit( 0);
}

int main( int argc, char *argv[])
{
  int i, j;
  int last_node_id, last_elem_id;
  SCoord * node;
  SElem * elem;
  int elemi[4];
  const char *analysis = "Analysis_of whatever";
  const char *range_table_name = "simple range table";
  char buf[ 1024];
  char *prefix = NULL;
  char *format = NULL;
  int skip_options = 0;
  int early_terminate = 0;

  for ( int ia = 1; ia < argc; ia++) {
    if ( !skip_options && ( argv[ ia][ 0] == '-')) {
      char opt = tolower( argv[ ia][ 1]);
      if ( opt == 'h') {
        print_help_and_exit( argv[ 0]);
      } else if ( opt == 'f') {
        ia++;
        if ( ia < argc) {
          format = strdup( argv[ ia]);
        } else {
          printf( "Missing arguments.\n");
          print_help_and_exit( argv[ 0]);
        }
      } else if ( opt == 'a' ) {
        early_terminate = 1;
      } else if ( opt == '-' ) {
        skip_options = 1;
      } else {
        printf( "Unknwon option '%s'.\n", argv[ ia]);
        print_help_and_exit( argv[ 0]);
      }
    } else {
      prefix = strdup( argv[ ia]);
      break;
    }
  }

  GiD_PostInit();
  GenMesh();

  // default values
  if ( !prefix) prefix = strdup( "test");
  if ( !format) format = strdup( "ascii");

  if ( !strcasecmp( format, "ascii")) {
    GiD_PostSetFormatReal("%.8g");
    sprintf( buf, "%s.post.msh", prefix);
    GiD_OpenPostMeshFile( buf, GiD_PostAscii);
    sprintf( buf, "%s.post.res", prefix);
    GiD_OpenPostResultFile( buf, GiD_PostAscii);
    printf( "Creating ASCII gid post files '%s.post.msh' and '%s.post.res'.\n",
            prefix, prefix);
  } else if ( !strcasecmp( format, "bin")) {
    sprintf( buf, "%s.post.bin", prefix);
    GiD_OpenPostResultFile( buf, GiD_PostBinary);
    printf( "Creating BINary gid post file '%s'.\n", buf);
  } else if ( !strcasecmp( format, "hdf5")) {
    sprintf( buf, "%s.post.h5", prefix);
    printf( "Creating HDF5 gid post file '%s'.\n", buf );
    GiD_OpenPostResultFile( buf, GiD_PostHDF5);
  } else {
    printf( "Unkown format '%s'.\n", format);
    print_help_and_exit( argv[ 0]);
  }
  
  /* write mesh info */
  GiD_BeginMeshColor("TestMesh", GiD_2D, GiD_Triangle, 3, 0, 0.99, 0);
  /* coordinates */
  GiD_BeginCoordinates();
  last_node_id = 0;
  for ( i = 0; i < 9; i++ ) {
    node = nodes + i;
    last_node_id = last_node_id < node->id ? node->id : last_node_id;
    GiD_WriteCoordinates( node->id, node->x, node->y, node->z ); 
  }
  GiD_EndCoordinates();
  /* elements */
  GiD_BeginElements();
  last_elem_id = 0;
  for ( i = 0; i < 8; i++ ) {
    elem = elems + i;
    elemi[0] = elem->n1;
    elemi[1] = elem->n2;
    elemi[2] = elem->n3;
    elemi[3] = 2;
    last_elem_id = last_elem_id < elem->id ? elem->id : last_elem_id;
    GiD_WriteElementMat(elem->id, elemi);
  }
  GiD_EndElements();
  GiD_EndMesh();

  /* write mesh info */
  GiD_BeginMeshColor("TestMesh 2", GiD_2D, GiD_Triangle, 3, 0, 0, 0.99);
  /* coordinates */
  GiD_BeginCoordinates();
  for ( i = 0; i < 9; i++ ) {
    node = nodes + i;
    GiD_WriteCoordinates( node->id + last_node_id, node->x, node->y, node->z + 3.0 ); 
  }
  GiD_EndCoordinates();
  /* elements */
  GiD_BeginElements();
  for ( i = 0; i < 8; i++ ) {
    elem = elems + i;
    elemi[0] = elem->n1 + last_node_id;
    elemi[1] = elem->n2 + last_node_id;
    elemi[2] = elem->n3 + last_node_id;
    GiD_WriteElement( elem->id + last_elem_id, elemi);
  }
  last_elem_id += last_elem_id;
  GiD_EndElements();
  GiD_EndMesh();

  /* write mesh info */
  GiD_BeginMeshColor( "Spheres", GiD_3D, GiD_Sphere, 1, 1.0, 0.7, 0.3 );
  /* coordinates */
  GiD_BeginCoordinates();
  GiD_EndCoordinates();
  /* elements */
  GiD_BeginElements();
  for ( i = 0; i < 9; i++ ) {
    double radius = Random() * 0.3 + 0.3;
    GiD_WriteSphere( last_elem_id + i + 1, i + 1, radius );
  }
  last_elem_id += 9;
  GiD_EndElements();
  GiD_EndMesh();

  GiD_BeginMeshColor( "Circles", GiD_3D, GiD_Circle, 1, 0.8, 1.0, 0.2 );
  /* coordinates */
  GiD_BeginCoordinates();
  GiD_EndCoordinates();
  /* elements */
  GiD_BeginElements();
  for ( i = 0; i < 9; i++ ) {
    double radius = Random() * 0.3 + 0.3;
    GiD_WriteCircle( last_elem_id + i + 1, i + 9 + 1, radius, 0.0, 0.0, 1.0 );
  }
  last_elem_id += 9;
  GiD_EndElements();
  GiD_EndMesh();

  // last_node_id += last_node_id;
  // last_elem_id += last_elem_id;

  /* now results info */
  /* write the gauss points */
  GiD_BeginGaussPoint("element_3gp", GiD_Triangle, NULL, 3, 0, 1);  
  GiD_EndGaussPoint();

  GiD_BeginGaussPoint("element_1gp", GiD_Triangle, NULL, 1, 0, 1);  
  GiD_EndGaussPoint();

  /* results are between 0.0 and 1.0 */
  GiD_BeginRangeTable( range_table_name);
  GiD_WriteMinRange( 0.1, "abs min");
  GiD_WriteRange( 0.1, 0.4, "min");
  GiD_WriteRange( 0.4, 0.6, "middle");
  GiD_WriteRange( 0.6, 0.9, "max");
  GiD_WriteMaxRange( 0.9, "abs max");
  GiD_EndRangeTable();
  
  GiD_BeginResult("debug nodes", analysis,
                  1.0, GiD_Scalar, GiD_OnNodes,
		  NULL, NULL, 0, NULL);
  for ( i = 0; i < 9; i++ ) { 
    GiD_WriteScalar(nodes[i].id, ( double)( nodes[i].id));
    // second mesh has last_node_id and last_elem_id offset
    GiD_WriteScalar(nodes[i].id + last_node_id, ( double)( nodes[i].id + last_node_id));
  }
  GiD_EndResult();

  GiD_BeginResult("debug elements", analysis,
                  1.0, GiD_Scalar, GiD_OnGaussPoints,
		  "element_1gp", NULL, 0, NULL);
  for ( i = 0; i < 8; i++ ) {
    for ( j = 0; j < 1; j++ )
      GiD_WriteScalar(elem[ i].id, ( double)elem[ i].id);
    // second mesh has last_node_id and last_elem_id offset
    for ( j = 0; j < 1; j++ )
      GiD_WriteScalar(elem[ i].id + last_elem_id, ( double)( elem[ i].id + last_elem_id));
  }
  GiD_EndResult();

  GiD_BeginResult("almost node_scalar (r)", analysis,
                  1.0, GiD_Scalar, GiD_OnNodes,
		  NULL, range_table_name, 0, NULL);
  for ( i = 0; i < 9; i++ ) { 
    GiD_WriteScalar(nodes[i].id, Random());
    // second mesh has last_node_id and last_elem_id offset
    GiD_WriteScalar(nodes[i].id + last_node_id, Random());
  }
  GiD_EndResult();
  
  /* scalar result over nodes */
#if defined(USE_GP_IN_GROUP)
  GiD_BeginResultGroup( analysis, 1.0, GiD_OnGaussPoints,"element_3gp" );
#else
  GiD_BeginResultGroup( analysis, 1.0, GiD_OnNodes, NULL);
#endif
  GiD_ResultDescription("EscalarNodosG", GiD_Scalar);
  GiD_ResultUnit( "m" );
  GiD_ResultRange( range_table_name );
  GiD_ResultDescriptionDim("VectorNodosG", GiD_Vector, 4);
  GiD_ResultUnit( "m/s" );
  GiD_ResultDescription("MatrixG", GiD_Matrix);
  GiD_ResultDescription("Local AxesG", GiD_LocalAxes);
  GiD_ResultUnit( "eugen" );
  for ( i = 0; i < sizeof(LOC)/sizeof(LOC[0]); i++ ) {
    for (j = 0; j < LPL; j++) {
      GiD_WriteScalar(LOC[i].id, Random());
      GiD_WriteVectorModule(LOC[i].id, Random(), Random(), Random(),-1);
      GiD_Write3DMatrix(LOC[i].id, Random(), Random(), Random(),
                        Random(), Random(), Random());
      GiD_WriteLocalAxes(LOC[i].id, Random(), Random(), Random());
    }
    // second mesh has last_node_id and last_elem_id offset
    for (j = 0; j < LPL; j++) {
      GiD_WriteScalar(LOC[i].id + LAST_ID, Random());
      GiD_WriteVectorModule(LOC[i].id + LAST_ID, Random(), Random(), Random(),-1);
      GiD_Write3DMatrix(LOC[i].id + LAST_ID, Random(), Random(), Random(),
                        Random(), Random(), Random());
      GiD_WriteLocalAxes(LOC[i].id + LAST_ID, Random(), Random(), Random());
    }
  }
  GiD_EndResult();
  
  /* #else */

  /* scalar results over gauss points: Triangle + 3 GP */
  
  GiD_BeginResult("Sobre_GP Triangle3", analysis, 1.0, GiD_Scalar, GiD_OnGaussPoints, 
		  "element_3gp", range_table_name, 0, NULL);
  for ( i = 0; i < 8; i++ ) {
    for ( j = 0; j < 3; j++ )
      GiD_WriteScalar(elem[ i].id, Random());
    // second mesh has last_node_id and last_elem_id offset
    for ( j = 0; j < 3; j++ )
      GiD_WriteScalar(elem[ i].id + last_elem_id, Random());
  }
  GiD_EndResult();
  GiD_BeginResult("Nodal_Scalar with spaces ( r2)",
                  analysis, 1.0, GiD_Scalar, GiD_OnNodes,
		  NULL, range_table_name, 0, NULL);
  for ( i = 0; i < 9; i++ ) {
    float rr = Random();
    GiD_WriteScalar(nodes[i].id, rr * rr);
    rr = Random();
    // second mesh has last_node_id and last_elem_id offset
    GiD_WriteScalar(nodes[i].id + last_node_id, rr * rr);
  }
  GiD_EndResult();

  /* vector result over nodes */
  
  /*
  GiD_BeginResult("VectorNodos", "Analysis", 1.0, GiD_Vector, GiD_OnNodes,
  NULL, NULL, 0, NULL);*/
  GiD_BeginResultHeader("VELOCITIES (m/s)", analysis,
                        1.0, GiD_Vector, GiD_OnNodes, NULL);
  for ( i = 0; i < 9; i++ ) { 
    GiD_WriteVector(nodes[i].id, Random(), Random(), Random());
    // second mesh has last_node_id and last_elem_id offset
    GiD_WriteVector(nodes[i].id + last_node_id, Random(), Random(), Random());
  }
  GiD_EndResult();

#ifndef DISABLE_CODE

#ifdef DO_FLUSH
  GiD_FlushPostFile();
#endif

#ifdef DO_FLUSH
  GiD_FlushPostFile();
#endif 
  /* matrix result */

  GiD_BeginResult("Matrix", analysis, 1.0, GiD_Matrix, GiD_OnNodes,
		  NULL, NULL, 0, NULL);
  for ( i = 0; i < 9; i++ ) {
    GiD_Write3DMatrix(nodes[ i].id, Random(), Random(), Random(),
		      Random(), Random(), Random());
    // second mesh has last_node_id and last_elem_id offset
    GiD_Write3DMatrix(nodes[ i].id + last_node_id, Random(), Random(), Random(),
		      Random(), Random(), Random());
  }
  GiD_EndResult();  

  /* local axes result */

  GiD_BeginResult("Local Axes", analysis, 1.0, GiD_LocalAxes, GiD_OnNodes,
		  NULL, NULL, 0, NULL);
  for ( i = 0; i < 9; i++ ) {
    GiD_WriteLocalAxes(nodes[ i].id, Random(), Random(), Random());
    // second mesh has last_node_id and last_elem_id offset
    GiD_WriteLocalAxes(nodes[ i].id + last_node_id, Random(), Random(), Random());
  }
  GiD_EndResult();

  GiD_BeginResult("Complex Scalar", analysis, 1.0, GiD_ComplexScalar, GiD_OnNodes,
		  NULL, NULL, 0, NULL);
  for ( i = 0; i < 9; i++ ) {
    GiD_WriteComplexScalar(nodes[ i].id, Random(), Random());
    // second mesh has last_node_id and last_elem_id offset
    GiD_WriteComplexScalar(nodes[ i].id + last_node_id, Random(), Random());
  }
  GiD_EndResult();

  GiD_BeginResult("Complex Vector (last result)", analysis, 1.0, GiD_ComplexVector, GiD_OnNodes,
		  NULL, NULL, 0, NULL);
  for ( i = 0; i < 9; i++ ) {
    GiD_WriteComplexVector( nodes[ i].id,
                            Random(), Random(), Random(),
                            Random(), Random(), Random());
    // second mesh has last_node_id and last_elem_id offset
    GiD_WriteComplexVector( nodes[ i].id + last_node_id,
                            Random(), Random(), Random(),
                            Random(), Random(), Random());
  }
  GiD_EndResult();
 
#endif /* DISABLE_CODE */

  GiD_FlushPostFile();
  if ( early_terminate ) {
    printf( "aborting...\n" );
    abort();
  }

  if ( !strcasecmp( format, "ascii")) {
    GiD_ClosePostMeshFile();
  }
  GiD_ClosePostResultFile();
  GiD_PostDone();
  printf( "... done.\n");
  
  free( prefix); prefix = NULL;
  free( format); format = NULL;
  
  return 0;
}
