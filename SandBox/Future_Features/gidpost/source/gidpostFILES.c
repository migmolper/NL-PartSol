/* gidpost */
/*
 *  gidpost.c--
 *
 *    This file implement a C interface for generating postprocessing
 *    results in the 'New postprocess format' of GiD. See declaration
 *    in gidpost.h
 *
 */

#include <assert.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <ctype.h>

#include "gidpostInt.h"
#include "gidpost.h"
#include "gidpostHash.h"

#include "gidpostFILES.h"

/* #define DEBUG_STATE */

#if defined( DEBUG_STATE )
#define GP_DUMP_STATE( File )                   \
  printf( "Current State = %s\n", GetStateDesc( CPostFile_TopState( File ) ) );
#else
#define GP_DUMP_STATE( File )
#endif

#ifdef WIN32
#define snprintf _snprintf
#endif
  
/* let's change the double quotes to simple ones, to avoid potential problems*/
static char *change_quotes(char *str) {
  unsigned int in;
  
  if ( str && *str) {
    for ( in = 0; in < strlen( str); in++) {
      if ( str[ in] == '"')
	str[ in] = '\'';
    }
  }
  return str;
}

/* to check if a string must be enclosed in quotes*/
int string_hasspace(const char *str){
  unsigned char *tstr;
  /*
    assert(str);
  */
  if ( !str || !*str) return 0;
  tstr=(unsigned char*)str;
  while(*tstr){
    if(isspace(*tstr))
      return 1;
    tstr++;
  }  
  return 0;
}


/* ---------------------------------------------------------------------------
 *
 *  Global files
 *
 * ---------------------------------------------------------------------------
 */

CPostFile *MeshFile = NULL;
CPostFile *ResultFile = NULL;
CPostFile *outputMesh = NULL;

/*
  hay que encapsular HDF5 en las estructuras internas de GiDPost
*/
GiD_PostMode PostMode;

static GP_CONST char * level_desc[] = {
  "UNDEFINED level",
  "TOP level",
  "MESH header",
  "MESHGROUP header",
  "inside a Coordinate block",
  "after a Coordinate block but inside a MESH",
  "inside an Element block",
  "GAUSS point block: implicit",
  "GAUSS point block: explicit",
  "RANGE table block",
  "OnGroup block",
  "Result block (deprecated)",
  "Result block (single)",
  "Result group block",
  "Result description block",
  "writing values",
  "unknown"
};

static GP_CONST char * GetStateDesc(int i)
{
  static int last = sizeof(level_desc)/sizeof(level_desc[0]) - 1;
  return (i<0 || i >= last) ? level_desc[last] : level_desc[i]; 
} 

CPostFile *_GiDfiles_GetMeshFile()
{
  return (outputMesh = MeshFile ? MeshFile : ResultFile);
}

CPostFile * _GiDfiles_NewFile(GiD_PostMode Mode)
{
  GiD_PostInit();
  switch ( Mode ) {
    case GiD_PostAscii:
      return CPostAscii_Create();
    case GiD_PostAsciiZipped:
      return CPostAsciiZ_Create();      
    case GiD_PostBinary:
      return CPostBinary_Create();
    case GiD_PostHDF5:
      return NULL;
    default:
      return NULL;
  }
}

static GP_CONST char * strElementType[]= {
  "",
  "Point",
  "Linear",
  "Triangle",
  "Quadrilateral",
  "Tetrahedra",
  "Hexahedra",
  "Prism",
  "Pyramid",
  "Sphere",
  "Circle"
};

GP_CONST char * GetElementTypeName( GiD_ElementType type )
{
  return strElementType[(int)(type)];
}

int _GiDfiles_CheckState( post_state s_req, CPostFile *file )
{
  post_state s_cur = CPostFile_TopState( file );
  if (s_req!=s_cur) {
    printf("invalid state '%s' should be '%s'\n", GetStateDesc(s_cur), GetStateDesc(s_req));
    return 0;
  }
  return 1;
}

static
int ValidateConnectivity(GiD_ElementType etype , int NNode)
{
  int error;

  switch (etype) {
  case GiD_Point:
  case GiD_Sphere:
  case GiD_Circle:
    error = (NNode != 1);
    break;
  case GiD_Linear:
    error = (NNode != 2 && NNode != 3);
    break;
  case GiD_Triangle:
    error = (NNode != 3 && NNode != 6);
    break;
  case GiD_Quadrilateral:
    error = (NNode != 4 && NNode != 8 && NNode != 9);
    break;
  case GiD_Tetrahedra:
    error = (NNode != 4 && NNode != 10);
    break;
  case GiD_Hexahedra:
    error = (NNode != 8 && NNode != 20 && NNode != 27);
    break;
  case GiD_Prism:
    error = (NNode != 6 && NNode != 15);
    break;
  case GiD_Pyramid:
    error = (NNode != 5 && NNode != 13);
    break;
   default:
    printf("invalid type of element code %d", etype);
    return 0;
  }
  if (error) {
    printf("invalid number of nodes '%d' for element type %s\n",
	   NNode, GetElementTypeName(etype));
    return 0;
  }
  return 1;
}


/*
 *  Begin a new mesh --
 *
 *    Start a mesh block, if MeshFile is opened write using this file
 *    in other case write using ResultFile
 */

int _GiDfiles_BeginMesh( CPostFile *File,
                    GP_CONST char * MeshName, GiD_Dimension Dim,
                    GiD_ElementType EType, int NNode )
{
  int fail = 1;
  char line[LINE_SIZE];
  char *mesh_name;
  post_state cur_state;
  /* here we sould validate EType & NNode */
  assert(ValidateConnectivity(EType,NNode));
  assert(File);
  //GP_DUMP_STATE( File );
  cur_state = CPostFile_TopState( File );
  if ( cur_state != POST_S0 && cur_state != POST_MESHGROUP_S0 )
    {
    return GP_ERROR_INVALID_STATE;
    }
    
  /* TODO:
     - return error code is validation fail
     - add check MESHGROUP vs MESH (!MESHGROUP)
   */
  /*
    assert(_GiDfiles_CheckState(POST_S0,File->level_mesh));
  */

  mesh_name = change_quotes(strdup(MeshName));
  snprintf(line, LINE_SIZE-1,
	   "MESH \"%s\" dimension %d ElemType %s Nnode %d",
	   mesh_name, (int)(Dim), GetElementTypeName(EType), NNode);
  free(mesh_name);
  if (!(fail = CPostFile_WriteString(File,line))) 
    {
    CPostFile_SetConnectivity( File, NNode );
    }
  CPostFile_PushState( File, POST_MESH_S0 );
  //File->level_mesh = POST_MESH_S0;
  return fail;
}

/*
 *  Begin a new mesh --
 *
 *    Start a mesh block, if MeshFile is opened write using this file
 *    in other case write using ResultFile. With this function you can
 *    specify a color for the mesh by its RGB components where each
 *    component take values on the interval [0,1].
 */

int _GiDfiles_BeginMeshColor(CPostFile *File,
		        GP_CONST char * MeshName, GiD_Dimension Dim,
		        GiD_ElementType EType, int NNode,
		        double Red, double Green, double Blue)
{
  int fail ;
  char line[LINE_SIZE];

  if ( !( fail = _GiDfiles_BeginMesh(File, MeshName, Dim, EType, NNode ) ) ) 
    {
    snprintf(line, LINE_SIZE-1, "# color %f %f %f", Red, Green, Blue);
    fail = CPostFile_WriteString( File, line );
    }
  return fail;
}

int _GiDfiles_MeshUnit(CPostFile *File,GP_CONST char * UnitName)
{
  if (UnitName) {
    char line[LINE_SIZE];
    char *tmp_name;
    tmp_name = change_quotes( strdup( UnitName));
    snprintf(line, LINE_SIZE-1, "Unit \"");
    strcat(line, tmp_name);
    strcat(line, "\"");
    free( tmp_name);
    return CPostFile_WriteString(File, line);
  }  
  return 1;
}

/*
 *  Start a coordinate block in the current mesh
 */


int _GiDfiles_BeginCoordinates( CPostFile *File )
{
  /* state checking */
  assert(File);
  assert( _GiDfiles_CheckState( POST_MESH_S0, File ) );
  CPostFile_PushState( File, POST_MESH_COORD0 );
  /* File->level_mesh = POST_MESH_COORD0; */
  return CPostFile_BeginCoordinates( File );
}

/*
 *  Close the current coordinate block
 */

int _GiDfiles_EndCoordinates( CPostFile *File )
{
  assert( PostMode!=GiD_PostHDF5);

  /* state checking */
  assert( _GiDfiles_CheckState( POST_MESH_COORD0, File ) );
  CPostFile_PopState( File );
  CPostFile_PushState( File, POST_MESH_COORD1 );
  // File->level_mesh = POST_MESH_COORD1;
  CPostFile_ResetLastID( File ); /* m_lastID is also checked when writting coordinates */
  return CPostFile_WriteString( File, "End Coordinates" );
}

/*
 * This function open a group of mesh. This makes possible specifying
 * multiples meshes within the group.
 */

int _GiDfiles_BeginMeshGroup( CPostFile *File, GP_CONST char* Name )
{
  int fail = 1;
  char line[LINE_SIZE];
  char *name;
    
  assert(File);
  assert(_GiDfiles_CheckState( POST_S0, File ) );
  
  name = change_quotes( strdup( Name ) );

  snprintf(line, LINE_SIZE-1, "Group \"%s\"", name);
  free(name);
  fail = CPostFile_WriteString(File, line);
  CPostFile_PushState( File, POST_MESHGROUP_S0 );
  //File->level_mesh = POST_S0;
  return fail;
}

int _GiDfiles_EndMeshGroup( CPostFile *fileMesh )
{
  assert( _GiDfiles_CheckState( POST_MESHGROUP_S0, fileMesh ) );  
  CPostFile_PopState( fileMesh );
  assert( _GiDfiles_CheckState( POST_S0, fileMesh ) );  
  return CPostFile_WriteString( fileMesh, "End Group" );
}


/*
 * This function close the previously opened group of mesh. See
 * GiD_BeginMeshGroup.
 */

int _GiDfiles_WriteCoordinates(CPostFile *File, int id, 
                          double x, double y, double z)
{
  int res;
  /* state checking */
  assert(_GiDfiles_CheckState( POST_MESH_COORD0, File));
  /* keep on the same level */
  res = CPostFile_WriteValuesVA( File, id, 3, x, y, z);
  CPostFile_ResetLastID( File ); /* inside _WriteValuesVA there is some trick to write gaussian  results */
  return res;
}

/*
 *  Write a coordinate member at the current Coordinates Block 
 */

int _GiDfiles_WriteCoordinates2D(CPostFile *File, int id, double x, double y)
{
  int res = 0;
  /* state checking */
  assert( _GiDfiles_CheckState( POST_MESH_COORD0, File ) );
  /* keep in the same level */
  res = CPostFile_IsBinary(File)
    ? CPostFile_WriteValuesVA(File,id, 3, x, y, 0.0)
    : CPostFile_WriteValuesVA(File, id, 2, x, y);
  CPostFile_ResetLastID(File); /* inside _WriteValuesVA there is some trick to write gaussian  results */
  return res;
}

/*
 *  Start a elements block in the current mesh
 */

int _GiDfiles_BeginElements( CPostFile *File )
{
  /* state checking */
  assert( _GiDfiles_CheckState( POST_MESH_COORD1, File ) );
  CPostFile_PopState( File );
  CPostFile_PushState( File, POST_MESH_ELEM );
  // File->level_mesh = POST_MESH_ELEM;
  return CPostFile_BeginElements(File);
}

/*
 *  Close the current elements block
 */

int _GiDfiles_EndElements(CPostFile *File)
{
  /* state checking */
  assert( _GiDfiles_CheckState( POST_MESH_ELEM, File ) );
  CPostFile_PopState( File );
  assert( _GiDfiles_CheckState( POST_MESH_S0, File ) );
  CPostFile_PopState( File );
  // File->level_mesh = POST_S0;
  return CPostFile_WriteString(File, "End Elements");
}

/*
 *  Write an element member at the current Elements Block. The rest of
 *  the arguments define the connectivity of the element and its
 *  number depends on the NNode parameter given in the previous call
 *  to BeginPostMesh.
 *  
 */

int _GiDfiles_WriteElement(CPostFile *File, int id, int nid[])
{
  /* state checking */
  assert( _GiDfiles_CheckState( POST_MESH_ELEM, File ) );
  /* keep on the same state */
  return CPostFile_WriteElement(File,
		                id, CPostFile_GetConnectivity(File), nid);
}

/*
 *  Write an element member at the current Elements Block. The rest of
 *  the arguments define the connectivity of the element and its
 *  number depends on the NNode parameter given in the previous call
 *  to BeginPostMesh. The last id correspond to a material number.
 *  
 */

int _GiDfiles_WriteElementMat(CPostFile *File, int id, int nid[])
{
  /* state checking */
  assert(_GiDfiles_CheckState(POST_MESH_ELEM, File));    
  /* keep on the same state */
  return CPostFile_WriteElement(File,
		                id, CPostFile_GetConnectivity(File)+1,
		                nid);
}

/*
 *  Write an sphere element member at the current Elements Block.
 *  An sphere element is defined by:
 *
 *     id: element id
 *
 *     nid: node center given by the node id specified previously in
 *          the coordinate block.
 *
 *     r : radius of the sphere element.
 *  
 */

int _GiDfiles_WriteSphere(CPostFile *File, int id, int nid, double r)
{
  /* state checking */
  assert(_GiDfiles_CheckState(POST_MESH_ELEM, File));    
  /* keep on the same state */
  CPostFile_WriteInteger(File, id, 0);
  CPostFile_WriteInteger(File, nid, 1);
  CPostFile_WriteDouble (File, r,  2);
  if (CPostFile_IsBinary(File)) {
    CPostFile_WriteInteger(File, 1, 1);    
  }
  return 0;
}

/*
 *  Write an sphere element member at the current Elements
 *  Block. Providing also a material identification.
 *  
 *  An sphere element is defined by:
 *
 *     id: element id
 *
 *     nid: node center given by the node id specified previously in
 *          the coordinate block.
 *
 *     r : radius of the sphere element.
 *
 *     mat: material identification.
 *  
 */

int _GiDfiles_WriteSphereMat(CPostFile * File, int id, int nid, double r, int mat)
{
  /* state checking */
  assert( _GiDfiles_CheckState( POST_MESH_ELEM, File ) );    
  /* keep on the same state */
  CPostFile_WriteInteger(File, id,  0);
  CPostFile_WriteInteger(File, nid, 1);
  CPostFile_WriteDouble (File, r,   1);
  CPostFile_WriteInteger(File, mat, 2);
  return 0;
}

/*
 *  Write a circle element member at the current Elements Block.
 *  A circle element is defined by:
 *
 *     id: element id
 *
 *     nid: node center given by the node id specified previously in
 *          the coordinate block.
 *
 *     r : radius of the circle element.      
 *
 *     nx, ny, nz : normal to the plane containing the circle.
 *  
 */

int _GiDfiles_WriteCircle(CPostFile *File,
		     int id, int nid, double r,
		     double nx, double ny, double nz)
{
  /* state checking */
  assert( _GiDfiles_CheckState( POST_MESH_ELEM, File ) );    
  /* keep on the same state */
  CPostFile_WriteInteger(File, id,  0);
  CPostFile_WriteInteger(File, nid, 1);
  CPostFile_WriteDouble (File, r,   1);
  CPostFile_WriteDouble (File, nx,  1);
  CPostFile_WriteDouble (File, ny,  1);
  CPostFile_WriteDouble (File, nz,  2);
  if (CPostFile_IsBinary(File)) {
    CPostFile_WriteInteger(File, 1, 1);    
  }
  return 0;
}

/*
 *  Write a circle element member at the current Elements
 *  Block. Providing also a material identification.
 *  
 *  A circle element is defined by:
 *
 *     id: element id
 *
 *     nid: node center given by the node id specified previously in
 *          the coordinate block.
 *
 *     r : radius of the circle element.      
 *
 *     nx, ny, nz : normal to the plane containing the circle.
 *  
 */

int _GiDfiles_WriteCircleMat(CPostFile *File, int id, int nid, double r,
		        double nx, double ny, double nz, int mat)
{
  /* state checking */
  assert( _GiDfiles_CheckState( POST_MESH_ELEM, File ) );    
  /* keep on the same state */
  CPostFile_WriteInteger(File, id,  0);
  CPostFile_WriteInteger(File, nid, 1);
  CPostFile_WriteDouble (File, r,   1);
  CPostFile_WriteDouble (File, nx,  1);
  CPostFile_WriteDouble (File, ny,  1);
  CPostFile_WriteDouble (File, nz,  1);
  CPostFile_WriteInteger(File, mat, 2);
  return 0;
}

/* ---------------------------------------------------------------------------
 *
 *  Post Result Interface
 *
 * ---------------------------------------------------------------------------
 */

/*
 *  Open a new post result file
 */

/*
 *  Close the current post result file
 */

/*
 *  Begin Gauss Points definition
 */

int _GiDfiles_BeginGaussPoint(CPostFile *File,
		         GP_CONST char * name, GiD_ElementType EType,
		         GP_CONST char * MeshName,
		         int GP_number, int NodesIncluded, int InternalCoord)
{
  char line[LINE_SIZE];
  char *gp_name;
  char *mesh_name;
  
  /* check state & validation */
  post_state st = CPostFile_TopState( File );
  if ( st != POST_S0 && st != POST_RESULT_ONGROUP )
    {
    return GP_ERROR_SCOPE;
    }

  gp_name = change_quotes(strdup( name ) );

  snprintf( line, LINE_SIZE-1,
	    "GaussPoints \"%s\" ElemType %s",
	    gp_name,
	    GetElementTypeName( EType ) );
  if ( MeshName && *MeshName ) 
    {
    mesh_name = change_quotes( strdup( MeshName ) );
    strcat( line, " \"") ;
    strcat( line, mesh_name );
    strcat( line, "\"" );
    free( mesh_name );
  }
  free( gp_name ); 
  gp_name = NULL;
  if ( CPostFile_WriteString( File, line ) )
    {
    return GP_ERROR_WRITESTRING;
    }
  snprintf( line, LINE_SIZE, "Number Of Gauss Points: %d", GP_number );
  if ( CPostFile_WriteString( File, line ) )
    {
    return GP_ERROR_WRITESTRING;
    }
  /* here we could save the number of GP in order to check at
     EndGaussPoint */
  File->GP_number_check = GP_number;
  if ( EType == GiD_Linear )
    {
    if ( NodesIncluded ) 
      {
      if ( CPostFile_WriteString(File, "  Nodes included" ) )
        {
	return GP_ERROR_WRITESTRING;
        }
      } 
    else if (CPostFile_WriteString( File, "  Nodes not included" ) )
      {
      return GP_ERROR_WRITESTRING;
      }
    }
  if ( InternalCoord ) 
    {
    if ( CPostFile_WriteString( File, "Natural Coordinates: Internal" ) )
      {
      return GP_ERROR_WRITESTRING;
      }
    CPostFile_PushState( File, POST_GAUSS_S0 );
    } 
  else 
    {
    if ( CPostFile_WriteString(File, "Natural Coordinates: Given" ) )
      {
      return GP_ERROR_WRITESTRING;
      }
    CPostFile_PushState( File, POST_GAUSS_GIVEN );
    /* here we can save the size of the coordinates to check later
       in WriteGaussPointXX*/
    }
  return GP_OK;
}

/*
 *  End current Gauss Points definition
 */

static int CheckGaussPointEnd(CPostFile* File)
{
  post_state st = CPostFile_TopState( File );
  if (st != POST_GAUSS_S0 && st != POST_GAUSS_GIVEN)
    {
    printf("Invalid call of GiD_EndGaussPoint. Current state is '%s' and should be '%s' or '%s'\n",
	   GetStateDesc( st ),
	   GetStateDesc( POST_GAUSS_S0 ),
	   GetStateDesc( POST_GAUSS_GIVEN ) );
    return 0;
    }
  return 1;
}

static int CheckGaussPointGiven( int written, int check )
{
  if (written !=  check)
    {
    printf( "missmatch in gauss point given, written %d and %d were required",
	    written, check );
    return 0;
  }
  return 1;
}

int _GiDfiles_EndGaussPoint( CPostFile *File )
{
  /* check state */
  post_state st;
  assert(CheckGaussPointEnd(File));
  st = CPostFile_TopState( File );
  if (st != POST_GAUSS_S0 && st != POST_GAUSS_GIVEN)
    {
      return GP_ERROR_SCOPE;
    }
#ifndef NDEBUG
  if ( CPostFile_TopState( File ) == POST_GAUSS_GIVEN )
    { 
    assert( CheckGaussPointGiven( File->gauss_written, File->GP_number_check ) );
    }
#endif
  CPostFile_PopState( File );
  st = CPostFile_TopState( File );
  if ( st != POST_S0 && st != POST_RESULT_ONGROUP )
    {
    return GP_ERROR_SCOPE;
    }
  File->GP_number_check = File->gauss_written = 0;
  return CPostFile_WriteString( File, "End GaussPoints" );
}

/*
 *  Write internal gauss point coordinate.
 */

int _GiDfiles_WriteGaussPoint2D( CPostFile *File, double x, double y )
{
  /* check state */
  assert(_GiDfiles_CheckState( POST_GAUSS_GIVEN, File ) );    
  if ( !CPostFile_Write2D( File, x, y ) ) 
    {
    ++File->gauss_written;
    return GP_OK;
    }
  return GP_ERROR_WRITEPOINT;
}

int _GiDfiles_WriteGaussPoint3D(CPostFile *File, double x, double y, double z)
{
  /* check state */
  assert(_GiDfiles_CheckState( POST_GAUSS_GIVEN, File ) );
  if (!CPostFile_Write3D( File, x, y, z ) ) 
    {
    ++File->gauss_written;
    return GP_OK;    
    }
  return GP_ERROR_WRITEPOINT;
}

/*
 *  Begin a Range Table definition
 */

int _GiDfiles_BeginRangeTable( CPostFile *File, GP_CONST char * name )
{
  char line[LINE_SIZE];
  char *rt_name;
  
  /* check & update state */
  assert( _GiDfiles_CheckState( POST_S0, File ) );
  
  rt_name = change_quotes( strdup( name ) );
  snprintf( line, LINE_SIZE-1, "ResultRangesTable \"%s\"", rt_name );
  free( rt_name );
  if (!CPostFile_WriteString( File, line ) ) 
    {
    CPostFile_PushState( File, POST_RANGE_S0 );
    return GP_OK;
    }
  return GP_ERROR_WRITESTRING;
}

/*
 *  End a Range Table definition
 */

int _GiDfiles_EndRangeTable( CPostFile *File )
{
  /* check & update state */
  assert( _GiDfiles_CheckState( POST_RANGE_S0, File ) );

  if ( !CPostFile_WriteString( File, "End ResultRangesTable" ) )
    {
    CPostFile_PopState( File );
    assert( _GiDfiles_CheckState( POST_S0, File ) );
    return GP_OK;
    }
  return GP_ERROR_WRITESTRING;
}

/*
 *  Write Range functions --
 *
 *   WriteMinRange : write a range with an implicit minimum value, the
 *   minimum absolute in the result set.
 *
 *   WriteRange : write an explicit range.
 *
 *   WritemaxRange: write a range with an implicit maximum value, the
 *   maximum absolute in the result set.
 */

int _GiDfiles_WriteMinRange( CPostFile* File, double max, GP_CONST char * name )
{
  char line[LINE_SIZE];
  char *tmp_name;
  static char local_format[100];
  static int create_format=1;
  if(create_format){
    sprintf(local_format," - %s : \"%%s\"",GiD_PostGetFormatReal());
    create_format=0;
  }
  /* check state */
  assert( _GiDfiles_CheckState( POST_RANGE_S0, File ) );

  tmp_name = change_quotes(strdup( name ) );
  snprintf( line, LINE_SIZE-1,local_format, max, tmp_name );
  free( tmp_name );
  return CPostFile_WriteString( File, line );
}

int _GiDfiles_WriteRange( CPostFile* File,
		     double min, double max, GP_CONST char * name )
{
  char line[LINE_SIZE];
  char *tmp_name;
  static char local_format[100];
  static int create_format=1;
  if(create_format){
    sprintf(local_format," %s - %s : \"%%s\"",GiD_PostGetFormatReal(),GiD_PostGetFormatReal());
    create_format=0;
  }
/* check state */
  assert(_GiDfiles_CheckState( POST_RANGE_S0, File ) );
  
  tmp_name = change_quotes(strdup( name ) );
  snprintf( line, LINE_SIZE-1,local_format, min, max, tmp_name );
  free( tmp_name );
  return CPostFile_WriteString( File, line );
}

int _GiDfiles_WriteMaxRange(CPostFile *File, double min, GP_CONST char * name)
{
  char line[LINE_SIZE];
  char *tmp_name;
  static char local_format[100];
  static int create_format=1;
  if(create_format){
    sprintf(local_format,"%s - : \"%%s\"",GiD_PostGetFormatReal());
    create_format=0;
  }
  /* check state */
  assert( _GiDfiles_CheckState( POST_RANGE_S0, File ) );
  
  tmp_name = change_quotes( strdup( name));
  snprintf(line, LINE_SIZE-1,local_format, min, tmp_name);
  free( tmp_name);
  return CPostFile_WriteString(File, line);
}

/*
 *  Begin Result Block
 */

int _GiDfiles_BeginResult(CPostFile *File,
		     GP_CONST char     *Result,
		     GP_CONST char     *Analysis,
		     double             step,
		     GiD_ResultType     Type,
		     GiD_ResultLocation Where,
		     GP_CONST char     *GaussPointsName,
		     GP_CONST char     *RangeTable, 
		     int                compc,
		     GP_CONST char      *compv[])
{
  char line[LINE_SIZE];
  const char * loc;
  char *res_name, *analysis_name, *tmp_name;
  int i;
  char step_string[100];

  /* check & change state */
  post_state st = CPostFile_TopState( File );
  if ( st != POST_S0 && st != POST_RESULT_ONGROUP )
    {
    return GP_ERROR_SCOPE;
    }

  //loc = ( Where == GiD_OnNodes ) ? "OnNodes" : "OnGaussPoints";
  switch (Where) {
  case GiD_OnNodes: 
    loc = "OnNodes";
    break;
  case GiD_OnGaussPoints: 
    loc = "OnGaussPoints";
    break;
  case GiD_OnNurbsLine: 
    loc = "OnNurbsLine";
    break;
  case GiD_OnNurbsSurface: 
    loc = "OnNurbsSurface";
    break;
  case GiD_OnNurbsVolume: 
    loc = "OnNurbsVolume";
    break;
  default:
    loc = "OnNodes";
    break;
  }
  res_name = change_quotes( strdup( Result ) );
  analysis_name = change_quotes( strdup( Analysis ) );
  sprintf(step_string,GiD_PostGetFormatStep(),step);
  snprintf(line,LINE_SIZE-1,"Result \"%s\" \"%s\" %s %s %s",res_name,analysis_name,step_string,GetResultTypeName(Type,0 ),loc);
  free( res_name );
  free( analysis_name );
  if ( Where == GiD_OnGaussPoints )
    {
    assert( GaussPointsName );
    tmp_name = change_quotes( strdup( GaussPointsName ) );
    strcat( line, " \"" );
    strcat( line, tmp_name );
    strcat( line, "\"" );
    free( tmp_name );
    }
  if ( CPostFile_WriteString( File, line ) )
    {
    return GP_ERROR_WRITESTRING;
    }
  if ( RangeTable ) 
    {
    tmp_name = change_quotes( strdup( RangeTable ) );
    snprintf( line, LINE_SIZE-1, "ResultRangesTable \"%s\"", tmp_name );
    free( tmp_name );
    if (CPostFile_WriteString( File, line ) )
      {
      return GP_ERROR_WRITESTRING;
      }
    }
  if (compc > 0) 
    {
    snprintf(line, LINE_SIZE-1, "ComponentNames");
    for (i = 0; i < compc; i++) 
      {  
      tmp_name = change_quotes(strdup(compv[ i]));
      strcat(line, " ");
      strcat(line, "\"");
      strcat(line, tmp_name);
      strcat(line, "\"");
      free(tmp_name);
      }
    if ( CPostFile_WriteString( File, line ) )
      {
      return GP_ERROR_WRITESTRING;
      }
    }
  CPostFile_PushState( File, POST_RESULT_DEPRECATED );
  File->flag_isgroup = 0;
  File->flag_begin_values = 0;
  return GP_OK;
}

int _GiDfiles_BeginResultHeader( CPostFile         *File,
                            GP_CONST char     *Result,
                            GP_CONST char     *Analysis,
                            double             step,
                            GiD_ResultType     Type,
                            GiD_ResultLocation Where,
                            GP_CONST char     *GaussPointsName)
{
  char lineHeader[LINE_SIZE];
  const char * loc;
  char *res_name, *analysis_name, *tmp_name;
  post_state st;
  char step_string[100];
  
  /* check & change state */
  assert( File );
  assert( Result );
  assert( Analysis );

  st = CPostFile_TopState( File );
  if ( st != POST_S0 && st != POST_RESULT_ONGROUP )
    {
    return GP_ERROR_SCOPE;
    }

  //loc = ( Where == GiD_OnNodes ) ? "OnNodes" : "OnGaussPoints";
  switch (Where) {
  case GiD_OnNodes: 
    loc = "OnNodes";
    break;
  case GiD_OnGaussPoints: 
    loc = "OnGaussPoints";
    break;
  case GiD_OnNurbsLine: 
    loc = "OnNurbsLine";
    break;
  case GiD_OnNurbsSurface: 
    loc = "OnNurbsSurface";
    break;
  case GiD_OnNurbsVolume: 
    loc = "OnNurbsVolume";
    break;
  default:
    loc = "OnNodes";
    break;
  }
  res_name = change_quotes(strdup( Result ) );
  analysis_name = change_quotes( strdup(Analysis ) );
  sprintf(step_string,GiD_PostGetFormatStep(),step);
  snprintf(lineHeader,LINE_SIZE-1,"Result \"%s\" \"%s\" %s %s %s",res_name,
    analysis_name,step_string,GetResultTypeName(Type,0),loc);
  free( res_name );
  free( analysis_name );
  if ( Where == GiD_OnGaussPoints )
    {
    assert( GaussPointsName );
    tmp_name = change_quotes( strdup( GaussPointsName ) );
    strcat( lineHeader, " \"" );
    strcat( lineHeader, tmp_name );
    strcat( lineHeader, "\"" );
    free( tmp_name);
  }
  if ( CPostFile_WriteString( File, lineHeader ) ) 
    {
    /* could not write result header */
    return GP_ERROR_WRITESTRING;
    }
  CPostFile_PushState( File, POST_RESULT_SINGLE );
  File->flag_isgroup = 0;
  File->flag_begin_values = 0;
  return 0;
}

static int CheckResultHeaderState( CPostFile* File )
{
  post_state st = CPostFile_TopState( File );
  if ( st == POST_RESULT_SINGLE || st == POST_RESULT_GROUP )
    return 1;
  printf("Invalid result state '%s'. Should be: '%s' or '%s'\n",
	 GetStateDesc( st ),
	 GetStateDesc( POST_RESULT_SINGLE ),
	 GetStateDesc( POST_RESULT_GROUP ) );
  return 0;
}

int _GiDfiles_ResultRange( CPostFile *File, GP_CONST char * RangeTable )
{
  char line[LINE_SIZE];
  char *tmp_name;

  /* check state */
  assert( File );
  assert( CheckResultHeaderState( File ) );
  assert( RangeTable );
  
  if ( RangeTable ) 
    {
    tmp_name = change_quotes( strdup(RangeTable ) );
    snprintf( line, LINE_SIZE-1, "ResultRangesTable \"%s\"", tmp_name );
    free( tmp_name );
    return CPostFile_WriteString( File, line );
    }
  return GP_ERROR_NULLSTRING;
}

int _GiDfiles_ResultComponents( CPostFile *File, int compc, GP_CONST char * compv[] )
{
  char line[LINE_SIZE];
  char *tmp_name;
  int i;
  
  /* check state */
  assert(File);
  assert( CheckResultHeaderState( File ) );
  assert(compc>0);
  
  if ( compc > 0 )
    {
    snprintf( line, LINE_SIZE-1, "ComponentNames" );
    for (i = 0; i < compc; i++) 
      {
      tmp_name = change_quotes( strdup( compv[ i ] ) );
      strcat( line, " " );
      strcat( line, "\"" );
      strcat( line, tmp_name );
      strcat( line, "\"" );
      free( tmp_name );
    }
    return CPostFile_WriteString(File, line);
  }
  return GP_ERROR_ZEROCOMPONENTS;
}

int _GiDfiles_ResultUnit( CPostFile *File, GP_CONST char * UnitName )
{
  /* check state */
  assert( CheckResultHeaderState( File ) );
  assert( UnitName );

  if ( UnitName )
    {
    char line[LINE_SIZE];
    char *tmp_name;
    tmp_name = change_quotes( strdup( UnitName ) );
    snprintf( line, LINE_SIZE-1, "Unit \"");
    strcat( line, tmp_name);
    strcat( line, "\"");
    free( tmp_name );
    return CPostFile_WriteString( File, line );
    }
  return GP_ERROR_NULLSTRING;
}

int _GiDfiles_ResultUserDefined( CPostFile *File, GP_CONST char * Name,GP_CONST char * Value )
  {  
  char line[LINE_SIZE];
  char* tmp_name=change_quotes( strdup( Name ) );
  char* tmp_value=change_quotes( strdup( Value ) );
  if( string_hasspace( Name ) )
    {
    if(string_hasspace(Value))
      {
      snprintf(line, LINE_SIZE-1, "# ResultUserDefined \"%s\" %s", tmp_name, tmp_value );  
      } 
    else 
      {
      snprintf(line, LINE_SIZE-1, "# ResultUserDefined \"%s\" \"%s\"", tmp_name, tmp_value );  
      }
    }
  else 
    {
    if( string_hasspace( Value ) )
      {
      snprintf(line, LINE_SIZE-1, "# ResultUserDefined %s %s", tmp_name, tmp_value );  
      } 
    else 
      {
      snprintf(line, LINE_SIZE-1, "# ResultUserDefined %s \"%s\"", tmp_name, tmp_value );  
      }
    }
  free(tmp_name);
  free(tmp_value);
  return CPostFile_WriteString(File, line);
}

int _GiDfiles_BeginResultGroup( CPostFile         *File,
                           GP_CONST char     *Analysis,
                           double             step,
                           GiD_ResultLocation Where,
                           GP_CONST char     *GaussPointsName )
{
  char line[LINE_SIZE];
  const char * loc;
  char *analysis_name, *tmp_name;
  char step_string[100];
  /* check & change state */
  assert(File);
  assert( _GiDfiles_CheckState( POST_S0, File ) );
  assert( Analysis );

  //loc = ( Where == GiD_OnNodes) ? "OnNodes" : "OnGaussPoints";
  switch (Where) {
  case GiD_OnNodes: 
    loc = "OnNodes";
    break;
  case GiD_OnGaussPoints: 
    loc = "OnGaussPoints";
    break;
  case GiD_OnNurbsLine: 
    loc = "OnNurbsLine";
    break;
  case GiD_OnNurbsSurface: 
    loc = "OnNurbsSurface";
    break;
  case GiD_OnNurbsVolume: 
    loc = "OnNurbsVolume";
    break;
  default:
    loc = "OnNodes";
    break;
  }
  analysis_name = change_quotes( strdup(Analysis ) );
  
  sprintf(step_string,GiD_PostGetFormatStep(),step);
  snprintf(line,LINE_SIZE-1,"ResultGroup \"%s\" %s %s",analysis_name,step_string,loc);
  free( analysis_name );
  if ( Where == GiD_OnGaussPoints )
    {
    assert( GaussPointsName );
    tmp_name = change_quotes( strdup( GaussPointsName ) );
    strcat( line, " \"" );
    strcat( line, tmp_name );
    strcat( line, "\"" );
    free( tmp_name );
    }
  if ( CPostFile_WriteString( File, line ) )
    {
    /* could not write result header */
    return GP_ERROR_WRITESTRING;
    }
  CPostFile_PushState( File, POST_RESULT_GROUP );
  /* initialize values counter */
  File->flag_isgroup = 1;
  File->flag_begin_values = 0;
  CPostFile_ResultGroupOnBegin( File );
  return 0;
}

static int CheckStateDesc( CPostFile* File )
{
  post_state st = CPostFile_TopState( File );
  if ( st == POST_RESULT_GROUP || st == POST_RESULT_DESC)
    return 1;
  printf( "Invalid result state '%s'. Should be: '%s' or '%s'\n",
          GetStateDesc( st ),
          GetStateDesc( POST_RESULT_GROUP ),
          GetStateDesc( POST_RESULT_DESC ) );
  return 0;
}

int _GiDfiles_ResultDescription_( CPostFile     *File,
		             GP_CONST char *Result,
		             GiD_ResultType Type,
		             size_t         s )
{
  char line[LINE_SIZE];
  char *tmp_name;

  /* check & change state */
  assert( File );
  //assert(CheckStateDesc(File));
  assert( _GiDfiles_CheckState( POST_RESULT_GROUP, File ) );
  assert( Result );

  if ( CPostFile_TopState( File ) != POST_RESULT_GROUP )
    {
    return GP_ERROR_NOTINGROUP;
    }
  line[0] = '\0';
  tmp_name = change_quotes( strdup( Result ) );
  snprintf( line, LINE_SIZE-1, "ResultDescription \"%s\" %s",
	   tmp_name, GetResultTypeName( Type, s ) );
  free( tmp_name );
  if (CPostFile_WriteString( File, line ) ) 
    {
    /* could not write result description */
    return GP_ERROR_WRITESTRING;
    }
  //File->level_res = POST_RESULT_DESC;
  /* update number of values to check per location */
  CPostFile_ResultGroupOnNewType( File, Type );
  return GP_OK;
}

/*
 *  End Result Block
 */

int _GiDfiles_EndResult(CPostFile *File)
{
  int _fail = 0;
  int status;
  post_state cur_state;

  /* check & change state */
  assert( File );
  assert( _GiDfiles_CheckState( POST_RESULT_VALUES, File ) );

  if ( File->flag_isgroup )
    {
    status = CPostFile_ResultGroupIsEmpty( File );
    assert( status );
    }

  _fail = CPostFile_EndValues( File );
  CPostFile_ResetLastID( File );
  GP_DUMP_STATE( File );
  CPostFile_PopState( File );
  GP_DUMP_STATE( File );
  CPostFile_PopState( File );
  GP_DUMP_STATE( File );
  cur_state = CPostFile_TopState( File );
  if ( cur_state != POST_S0 && cur_state != POST_RESULT_ONGROUP )
    {
    assert( cur_state == POST_S0 || cur_state == POST_RESULT_ONGROUP );
    return GP_ERROR_SCOPE;
    }
  return _fail;
}

/*
 * Results which are referred to a mesh group (see GiD_BeginMeshGroup)
 * should be written between a call to this function and
 * GiD_EndOnMeshGroup.
 */

int _GiDfiles_BeginOnMeshGroup( CPostFile *File, char *Name )
{
  char line[LINE_SIZE];
  
  /* check state */
  assert( _GiDfiles_CheckState( POST_S0, File ) );

  if ( CPostFile_TopState( File ) != POST_S0 )
    {
    return GP_ERROR_SCOPE;
    }
  snprintf(line, LINE_SIZE-1, "OnGroup \"%s\"", Name);
  CPostFile_PushState( File, POST_RESULT_ONGROUP );
  return CPostFile_WriteString( File, line );
}

int _GiDfiles_EndOnMeshGroup( CPostFile *File )
{
  /* check state */
  assert( _GiDfiles_CheckState( POST_RESULT_ONGROUP, File ) );

  CPostFile_PopState( File );
  assert( _GiDfiles_CheckState( POST_S0, File ) );
  return CPostFile_WriteString( File, "End OnGroup" );
}

/*
 * This function close a previously opened block of result over a mesh
 * group.
 */

/*
 *  Write result functions
 */

static
int _GiDfiles_internal_EnsureBeginValues( CPostFile *File )
{
  post_state st;

  assert( File );
  
  if ( !File->flag_begin_values ) 
    {
    if ( !CPostFile_BeginValues( File ) ) 
      {
      st = CPostFile_TopState( File );
      
      CPostFile_PushState( File, POST_RESULT_VALUES );
      if ( File->flag_isgroup )
        {
	CPostFile_ResultGroupOnBeginValues( File );
        }
      File->flag_begin_values = 1;
      return 0;
      }
    }
  else
    {
    assert( _GiDfiles_CheckState( POST_RESULT_VALUES, File ) );
    }
  return 1;
}

int _GiDfiles_WriteScalar(CPostFile *File, int id, double v )
{
  _GiDfiles_internal_EnsureBeginValues( File );
  /* check state */
  assert( _GiDfiles_CheckState( POST_RESULT_VALUES, File ) );

  return File->flag_isgroup ?
    CPostFile_ResultGroupWriteValues( File, GiD_Scalar, id, 1, v) :
    CPostFile_WriteValues(File, id, 1, &v);
}

int _GiDfiles_Write2DVector(CPostFile *File, int id, double x, double y)
{
  _GiDfiles_internal_EnsureBeginValues(File);
  /* check state */
  assert(_GiDfiles_CheckState( POST_RESULT_VALUES, File ) );

  if (File->flag_isgroup || !CPostFile_IsBinary(File)) {
    return File->flag_isgroup ?
      CPostFile_ResultGroupWriteValues(File,
		                       GiD_Vector, id, 2, x, y) :
      CPostFile_WriteValuesVA(File, id, 2, x, y);
  } else {
    /* single result & binary */
    double mod;
    if(((float)x)!=GP_UNKNOWN){
      mod=sqrt(x*x + y*y);
    } else {
      assert(((float)y)==GP_UNKNOWN);
      mod=GP_UNKNOWN;
    }
    return CPostFile_WriteValuesVA( File, id, 4, x, y, 0.0, mod );
  }  
}

int _GiDfiles_WriteVector( CPostFile *File, int id, double x, double y, double z )
{
  _GiDfiles_internal_EnsureBeginValues( File );
  /* check state */
  assert(_GiDfiles_CheckState( POST_RESULT_VALUES, File ) );

  if ( File->flag_isgroup || !CPostFile_IsBinary( File ) ) {
    return File->flag_isgroup
      ?
      CPostFile_ResultGroupWriteValues( File,
                                        GiD_Vector, id, 3, x, y, z )
      :
      CPostFile_WriteValuesVA( File, id, 3, x, y, z );
  } else {
    /* single result & binary */
    double mod;
    if(((float)x)!=GP_UNKNOWN){
      mod=sqrt(x*x + y*y + z*z);
    } else {
      assert(((float)y)==GP_UNKNOWN && ((float)z)==GP_UNKNOWN);      
      mod=GP_UNKNOWN;
    }
    return CPostFile_WriteValuesVA( File, id, 4, x, y, z, mod );
  }  
}

int _GiDfiles_WriteVectorModule( CPostFile *File,
                            int id, double x, double y, double z, double mod )
{
  _GiDfiles_internal_EnsureBeginValues(File);
  /* check state */
  assert(_GiDfiles_CheckState( POST_RESULT_VALUES, File ) );

  /* 4-vectors can not be written on RG-ASCII */
  return File->flag_isgroup
    ?
    CPostFile_ResultGroupWriteValues(File,
		                     GiD_Vector, id, 4, x, y, z, mod)
    :
    CPostFile_WriteValuesVA(File, id, 4, x, y, z, mod );
}

int _GiDfiles_Write2DMatrix(CPostFile *File,
		       int id, double Sxx, double Syy, double Sxy)
{
  assert(File);
  if (CPostFile_IsBinary(File)) {
    return _GiDfiles_Write3DMatrix(File,
		              id,
		              Sxx, Syy,         0.0 /*Szz*/,
		              Sxy, 0.0 /*Syz*/, 0.0 /*Sxz*/);
  }
  _GiDfiles_internal_EnsureBeginValues(File);
  /* check state */
  assert(_GiDfiles_CheckState(POST_RESULT_VALUES, File ) );
  return File->flag_isgroup
    ?
    CPostFile_ResultGroupWriteValues(File,
		                     GiD_Matrix, id, 3, Sxx, Syy, Sxy)
    :
    CPostFile_WriteValuesVA(File, id, 3, Sxx, Syy, Sxy);
}

int _GiDfiles_Write3DMatrix(CPostFile *File,
		       int id, double Sxx, double Syy, double Szz,
		       double Sxy, double Syz, double Sxz)
{
  _GiDfiles_internal_EnsureBeginValues(File);
  /* check state */
  assert( _GiDfiles_CheckState( POST_RESULT_VALUES, File ) );
  
  return File->flag_isgroup
    ?
    CPostFile_ResultGroupWriteValues(File,
		                     GiD_Matrix,
		                     id, 6, Sxx, Syy, Szz, Sxy, Syz, Sxz)
    :
    CPostFile_WriteValuesVA(File, id, 6, Sxx, Syy, Szz, Sxy, Syz, Sxz);
}

int _GiDfiles_WritePlainDefMatrix(CPostFile *File, int id,
		             double Sxx, double Syy, double Sxy, double Szz )
{
#if 0
    if (CPostFile_IsBinary(File)) {
      return _GiDfiles_Write3DMatrix(File, id,
                                Sxx, Syy,         Szz,
                                Sxy, 0.0 /*Syz*/, 0.0 /*Sxz*/);
    }
#endif
  _GiDfiles_internal_EnsureBeginValues(File);
  /* check state */
  assert( _GiDfiles_CheckState( POST_RESULT_VALUES, File ) );
   
  return File->flag_isgroup
    ? 
    CPostFile_ResultGroupWriteValues(File,
		                     GiD_PlainDeformationMatrix,
		                     id, 4, Sxx, Syy, Sxy, Szz)
    :
    CPostFile_WriteValuesVA(File, id, 4, Sxx, Syy, Sxy, Szz);
}

int _GiDfiles_WriteMainMatrix(CPostFile *File, int id,
		         double Si, double Sii, double Siii,
		         double Vix, double Viy, double Viz,
		         double Viix, double Viiy, double Viiz,
		         double Viiix, double Viiiy, double Viiiz)
{
  _GiDfiles_internal_EnsureBeginValues(File);
  /* check state */
  assert( _GiDfiles_CheckState( POST_RESULT_VALUES, File ) );
  
  return File->flag_isgroup
    ? 
    CPostFile_ResultGroupWriteValues(File,
		                     GiD_MainMatrix, id, 12,
		                     Si,    Sii,   Siii,
		                     Vix,   Viy,   Viz,
		                     Viix,  Viiy,  Viiz,
		                     Viiix, Viiiy, Viiiz)
    :
    CPostFile_WriteValuesVA(File, id, 12,
		            Si,    Sii,   Siii,
		            Vix,   Viy,   Viz,
		            Viix,  Viiy,  Viiz,
		            Viiix, Viiiy, Viiiz);
}

int _GiDfiles_WriteLocalAxes(CPostFile *File,
		        int id, double euler_1, double euler_2, double euler_3)
{
  _GiDfiles_internal_EnsureBeginValues(File);
  /* check state */
  assert(_GiDfiles_CheckState(POST_RESULT_VALUES, File ) );
  
  return File->flag_isgroup
    ?
    CPostFile_ResultGroupWriteValues(File,
		                     GiD_LocalAxes, id, 3,
		                     euler_1, euler_2, euler_3)
    :
    CPostFile_WriteValuesVA(File, id, 3, euler_1, euler_2, euler_3);
}

/*
 * Complex numbers
 */

int _GiDfiles_WriteComplexScalar( CPostFile *File,
                             int id, double complex_real, double complex_imag) {
  _GiDfiles_internal_EnsureBeginValues( File);
  /* check state */
  assert( _GiDfiles_CheckState( POST_RESULT_VALUES, File ) );
  
  return File->flag_isgroup
    ?
    CPostFile_ResultGroupWriteValues( File,
				      GiD_ComplexScalar, id, 2,
				      complex_real, complex_imag)
    :
    CPostFile_WriteValuesVA(File, id, 2, complex_real, complex_imag);
}

int _GiDfiles_Write2DComplexVector( CPostFile *File, int id,
			       double x_real, double x_imag,
			       double y_real, double y_imag) {  
  double mod_r,mod_i,mod;
  if(((float)x_real)!=GP_UNKNOWN){
    double mod2_r=x_real * x_real + y_real * y_real;
    double mod2_i=x_imag * x_imag + y_imag * y_imag;
    mod_r=sqrt(mod2_r);
    mod_i=sqrt(mod2_i);
    mod=sqrt(mod2_r + mod2_i);
  } else {
    assert(((float)y_real)==GP_UNKNOWN && ((float)x_imag)==GP_UNKNOWN && ((float)y_imag)==GP_UNKNOWN);
    mod_r=mod_i=mod=GP_UNKNOWN;
  }
  
  _GiDfiles_internal_EnsureBeginValues( File);
  /* check state */
  assert( _GiDfiles_CheckState( POST_RESULT_VALUES, File ) );
  
  return File->flag_isgroup
    ?
    CPostFile_ResultGroupWriteValues( File,
				      GiD_ComplexVector, id,  9,
				      x_real, x_imag, y_real, y_imag, 0.0, 0.0,
				      mod_r, mod_i, mod)
    :
    CPostFile_WriteValuesVA(File, id, 9, 
			    x_real, x_imag, y_real, y_imag, 0.0, 0.0,
			    mod_r, mod_i, mod);
}

int _GiDfiles_WriteComplexVector( CPostFile *File, int id,
			     double x_real, double x_imag,
			     double y_real, double y_imag,
			     double z_real, double z_imag) {
  double mod_r,mod_i,mod;
  if(((float)x_real)!=GP_UNKNOWN){
    double mod2_r=x_real * x_real + y_real * y_real + z_real * z_real;
    double mod2_i=x_imag * x_imag + y_imag * y_imag + z_imag * z_imag;
    mod_r=sqrt(mod2_r);
    mod_i=sqrt(mod2_i);
    mod=sqrt(mod2_r + mod2_i);
  } else {
    assert(((float)y_real)==GP_UNKNOWN && ((float)z_real)==GP_UNKNOWN && ((float)x_imag)==GP_UNKNOWN && ((float)y_imag)==GP_UNKNOWN && ((float)z_imag)==GP_UNKNOWN);
    mod_r=mod_i=mod=GP_UNKNOWN;
  }  
  
  _GiDfiles_internal_EnsureBeginValues( File);
  /* check state */
  assert( _GiDfiles_CheckState( POST_RESULT_VALUES, File ) );

  return File->flag_isgroup
    ?
    CPostFile_ResultGroupWriteValues( File,
				      GiD_ComplexVector, id, 9,
				      x_real, x_imag, y_real, y_imag, z_real, z_imag,
				      mod_r, mod_i, mod)
    :
    CPostFile_WriteValuesVA(File, id, 9, 
			    x_real, x_imag, y_real, y_imag, z_real, z_imag,
			    mod_r, mod_i, mod);
}

int _GiDfiles_Write2DComplexMatrix(CPostFile *File, int id,
                                   double Sxx_real, double Syy_real, double Sxy_real,
                                   double Sxx_imag, double Syy_imag, double Sxy_imag) {
  assert(File);
  if (CPostFile_IsBinary(File)) {
    return _GiDfiles_WriteComplexMatrix(File,
                                        id,
                                        Sxx_real, Syy_real,         0.0 /*Szz_real*/,
                                        Sxy_real, 0.0 /*Syz_real*/, 0.0 /*Sxz_real*/,
                                        Sxx_imag, Syy_imag,         0.0 /*Szz_imag*/,
                                        Sxy_imag, 0.0 /*Syz_imag*/, 0.0 /*Sxz_imag*/);
  }
  _GiDfiles_internal_EnsureBeginValues(File);
  /* check state */
  assert(_GiDfiles_CheckState(POST_RESULT_VALUES, File ) );
  return File->flag_isgroup
      ?
      CPostFile_ResultGroupWriteValues(File,
                                       GiD_ComplexMatrix, id, 6,
                                       Sxx_real, Syy_real, Sxy_real, Sxx_imag, Syy_imag, Sxy_imag)
      :
      CPostFile_WriteValuesVA(File, id, 6,
                              Sxx_real, Syy_real, Sxy_real, Sxx_imag, Syy_imag, Sxy_imag);
}

int _GiDfiles_WriteComplexMatrix(CPostFile *File, int id,
                                 double Sxx_real, double Syy_real, double Szz_real,
                                 double Sxy_real, double Syz_real, double Sxz_real,
                                 double Sxx_imag, double Syy_imag, double Szz_imag,
                                 double Sxy_imag, double Syz_imag, double Sxz_imag) {
  _GiDfiles_internal_EnsureBeginValues(File);
  /* check state */
  assert( _GiDfiles_CheckState( POST_RESULT_VALUES, File ) );
  
  return File->flag_isgroup
      ?
      CPostFile_ResultGroupWriteValues(File, GiD_ComplexMatrix, id, 12, 
                                       Sxx_real, Syy_real, Szz_real,
                                       Sxy_real, Syz_real, Sxz_real,
                                       Sxx_imag, Syy_imag, Szz_imag,
                                       Sxy_imag, Syz_imag, Sxz_imag)
      :
      CPostFile_WriteValuesVA(File, id, 12, 
                              Sxx_real, Syy_real, Szz_real,
                              Sxy_real, Syz_real, Sxz_real,
                              Sxx_imag, Syy_imag, Szz_imag,
                              Sxy_imag, Syz_imag, Sxz_imag);
}

/* 
* Nurbs 
*/

int _GiDfiles_WriteNurbsSurface( CPostFile *File, int id, int n, double* v )
{
  _GiDfiles_internal_EnsureBeginValues( File );
  /* check state */
  assert(_GiDfiles_CheckState( POST_RESULT_VALUES, File ) );
  
  return File->flag_isgroup ?
    CPostFile_ResultGroupWriteValues( File, GiD_Scalar, id, 1, v) :
    CPostFile_WriteValuesNS(File, id, n, v);
}

int _GiDfiles_WriteNurbsSurfaceVector( CPostFile *File, int id, int n, int num_comp, double* v )
{
  _GiDfiles_internal_EnsureBeginValues( File );
  /* check state */
  assert(_GiDfiles_CheckState( POST_RESULT_VALUES, File ) );
  
  return File->flag_isgroup ?
    CPostFile_ResultGroupWriteValues( File, GiD_Vector, id, 3, v) :
    CPostFile_WriteValuesNSV(File, id, n, num_comp, v);
}
