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
#include <assert.h>

#include "gidpostInt.h"
#include "gidpost.h"
#include "gidpostHash.h"
#include "gidpostFILES.h"

#ifdef HDF5
  #include "gidpostHDF5.h"
#endif

/* defined in gidpostFILES.c */
extern CPostFile *MeshFile;
extern CPostFile *ResultFile;
extern CPostFile *outputMesh;
extern GiD_PostMode PostMode;


#define FD2FILE(_fd_, _File_)                     \
  do {                                            \
    assert(_fd_);                                 \
    _File_ = (CPostFile*)GiD_HashFind(_fd_);      \
    if (!_File_) {                                \
      /* file handler does not exists */          \
      return -8;                                  \
    }                                             \
  } while (0);


int GiD_PostInit()
{
  return GiD_HashInit();
}

int GiD_PostDone()
{
  return GiD_HashDone();
}

/* ---------------------------------------------------------------------------
 *
 *  Post Mesh Interface
 *
 * ---------------------------------------------------------------------------
 */

/*
 *  Open a new post mesh file
 */

int GiD_OpenPostMeshFile(GP_CONST char * FileName, GiD_PostMode Mode )
{
  PostMode=Mode;

  if(PostMode==GiD_PostHDF5){
#ifdef HDF5
      return GiD_OpenPostMeshFile_HDF5(FileName);
#else
      return 5;
#endif
    }
  /* Binary mode not allowed in mesh file!!! */
  assert(Mode!=GiD_PostBinary);
  assert(!MeshFile);
  if (MeshFile) {
    /* must be closed */
    return GP_ERROR_FILEOPENED;
  }
  if (!(MeshFile = _GiDfiles_NewFile(Mode)))
    /* not enough memory */
    return GP_ERROR_NOMEM;
  _GiDfiles_GetMeshFile();
  /* now outputMesh points to MeshFile */
  if (CPostFile_Open(MeshFile, FileName)) {
    /* Open failed */
    CPostFile_Release(MeshFile);
    MeshFile = NULL;
    return GP_ERROR_OPENFAILED;
  }
  
  //MeshFile->level_mesh = POST_S0;
  CPostFile_PushState( MeshFile, POST_S0 );
  return GP_OK;
}

GiD_FILE GiD_fOpenPostMeshFile( GP_CONST char * FileName,
		                GiD_PostMode Mode )
{
  CPostFile *File = NULL;
  GiD_FILE fd;
  
#ifdef HDF5
  if ( PostMode == GiD_PostHDF5 ) {
    return -1;
  }
#endif
  
  /* Binary mode not allowed in mesh file!!! */
  assert( Mode != GiD_PostBinary );

  if ( !( File = _GiDfiles_NewFile( Mode ) ) )
    /* not enough memory = -2 */
    return 0;
  /* now open the File */
  if ( CPostFile_Open( File, FileName ) ) {
    /* Open failed = -4 */
    return 0;
  }
  fd = GiD_HashAdd( File );
  if (!fd) {
    /* could not create a file handler = -6 */
    CPostFile_Release(File);
    return 0;
  }
  CPostFile_PushState( File, POST_S0 );
  return fd;
}

/*
 *  Close the current post mesh file
 */
int GiD_ClosePostMeshFile()
{
  int fail = 1;

#ifdef HDF5
  if(PostMode==GiD_PostHDF5){
    return GiD_ClosePostMeshFile_HDF5( );
  }
#endif
  assert( MeshFile );
  assert( _GiDfiles_CheckState( POST_S0, MeshFile ) );    
  
  if ( MeshFile ) {
    fail = CPostFile_Release( MeshFile );
    MeshFile = NULL;
    /* reset outpuMesh */
    _GiDfiles_GetMeshFile( );
  }
  return fail;
}

int GiD_fClosePostMeshFile( GiD_FILE fd )
{
  int fail = 1;
  CPostFile *File = NULL;
  
#ifdef HDF5
  if ( PostMode == GiD_PostHDF5 ) {
    return -1;
  }
#endif

  FD2FILE(fd,File);
  
  fail = CPostFile_Release( File );
  GiD_HashRemove( fd );
  
  return fail;
}

/*
 *  Begin a new mesh --
 *
 *    Start a mesh block, if MeshFile is opened write using this file
 *    in other case write using ResultFile
 */

int GiD_BeginMesh(GP_CONST char * MeshName, GiD_Dimension Dim,
		  GiD_ElementType EType, int NNode)
{
  CPostFile *mesh;
  
#ifdef HDF5
  if(PostMode==GiD_PostHDF5){
    return GiD_BeginMesh_HDF5(MeshName, Dim,EType,NNode);
  }
#endif

  mesh = _GiDfiles_GetMeshFile();
  return _GiDfiles_BeginMesh(mesh, MeshName, Dim, EType, NNode);
  
}

int GiD_fBeginMesh(GiD_FILE fd, GP_CONST char * MeshName,
		   GiD_Dimension Dim, GiD_ElementType EType, int NNode)
{
  CPostFile *mesh;
  
#ifdef HDF5
  if ( PostMode == GiD_PostHDF5 ) {
    return -1;
  }
#endif

  FD2FILE(fd,mesh);
  return _GiDfiles_BeginMesh(mesh, MeshName, Dim, EType, NNode);
}
/*
 *  Begin a new mesh --
 *
 *    Start a mesh block, if MeshFile is opened write using this file
 *    in other case write using ResultFile. With this function you can
 *    specify a color for the mesh by its RGB components where each
 *    component take values on the interval [0,1].
 */

int GiD_BeginMeshColor(GP_CONST char * MeshName, GiD_Dimension Dim,
		       GiD_ElementType EType, int NNode,
		       double Red, double Green, double Blue)
{
  CPostFile * mesh;

#ifdef HDF5
  if( PostMode == GiD_PostHDF5 ){
    return GiD_BeginMeshColor_HDF5( MeshName, Dim, EType, NNode, Red, Green, Blue );
  }
#endif

  mesh = _GiDfiles_GetMeshFile();
  return _GiDfiles_BeginMeshColor( mesh, MeshName, Dim, EType, NNode,
		             Red, Green, Blue );
}

int GiD_fBeginMeshColor(GiD_FILE fd, GP_CONST char * MeshName,
		        GiD_Dimension Dim, GiD_ElementType EType,
		        int NNode,
		        double Red, double Green, double Blue)
{
  CPostFile * mesh;
  
#ifdef HDF5
  if ( PostMode == GiD_PostHDF5 ) {
    return -1;
  }
#endif

  FD2FILE(fd,mesh);
    
  return _GiDfiles_BeginMeshColor(mesh, MeshName, Dim, EType, NNode,
		             Red, Green, Blue);
}

/*
 *  End current mesh
 */

int GiD_EndMesh()
{
  /* check & change state */
  /*
  assert(_GiDfiles_CheckState(POST_MESH_ELEM,level_mesh));
  level_mesh = POST_S0;
  */
  
#ifdef HDF5
  if(PostMode==GiD_PostHDF5){
    return GiD_EndMesh_HDF5();
  }
#endif

  
  return 0;
}

int GiD_fEndMesh(GiD_FILE fd)
{  
#ifdef HDF5
  if ( PostMode == GiD_PostHDF5 ) {
    return -1;
  }
#endif

  return 0;
}

int GiD_MeshUnit(GP_CONST char * UnitName)
{
#ifdef HDF5
  if(PostMode==GiD_PostHDF5){
    return  GiD_MeshUnit_HDF5(UnitName);
  }
#endif
  return _GiDfiles_MeshUnit( _GiDfiles_GetMeshFile( ), UnitName );
  return 1;
}


int GiD_fMeshUnit(GiD_FILE fd,GP_CONST char * UnitName)
{
  CPostFile *File = NULL;  
#ifdef HDF5
  if ( PostMode == GiD_PostHDF5 ) {
    return -1;
  }
#endif

  FD2FILE(fd,File);
  return _GiDfiles_MeshUnit(File,UnitName);
}

int GiD_MeshLocalAxes(GP_CONST char * Result, GP_CONST char * Analysis,double step)
{
#ifdef HDF5
    if(PostMode==GiD_PostHDF5){
    return GiD_MeshLocalAxes_HDF5(Result,Analysis,step);
  }
#endif
  return 1;
}

/*
 *  Start a coordinate block in the current mesh
 */

int GiD_BeginCoordinates()
{

#ifdef HDF5
  if(PostMode==GiD_PostHDF5){
    return GiD_BeginCoordinates_HDF5();
  }
#endif

  return _GiDfiles_BeginCoordinates( _GiDfiles_GetMeshFile( ) );
}

int GiD_fBeginCoordinates(GiD_FILE fd)
{
  CPostFile *File = NULL;
  
#ifdef HDF5
  if ( PostMode == GiD_PostHDF5 ) {
    return -1;
  }
#endif

  FD2FILE(fd,File);
  return _GiDfiles_BeginCoordinates( File );
}

/*
 *  Close the current coordinate block
 */

int GiD_EndCoordinates()
{
#ifdef HDF5
  if(PostMode==GiD_PostHDF5){
    return GiD_EndCoordinates_HDF5();
  }
#endif

  return _GiDfiles_EndCoordinates( _GiDfiles_GetMeshFile( ) );
}

int GiD_fEndCoordinates(GiD_FILE fd)
{
  CPostFile *File = NULL;
  
#ifdef HDF5
  if ( PostMode == GiD_PostHDF5 ) {
    return -1;
  }
#endif

  FD2FILE(fd,File);
  
  return _GiDfiles_EndCoordinates(File);
}

/*
 * This function open a group of mesh. This makes possible specifying
 * multiples meshes within the group.
 */

int GiD_BeginMeshGroup( GP_CONST char* Name )
{
#ifdef HDF5
  if(PostMode==GiD_PostHDF5){
    return GiD_BeginMeshGroup_HDF5( Name);
  }
#endif

  return _GiDfiles_BeginMeshGroup( _GiDfiles_GetMeshFile( ), Name );
}

int GiD_fBeginMeshGroup(GiD_FILE fd, GP_CONST char* Name)
{
  CPostFile *File = NULL;
  
#ifdef HDF5
  if ( PostMode == GiD_PostHDF5 ) {
    return -1;
  }
#endif

  FD2FILE(fd,File);
  
  return _GiDfiles_BeginMeshGroup(File, Name);
}

/*
 * This function close the previously opened group of mesh. See
 * GiD_BeginMeshGroup.
 */

int GiD_EndMeshGroup()
{
#ifdef HDF5
  if( PostMode == GiD_PostHDF5 ) {
    return GiD_EndMeshGroup_HDF5();
    }
#endif
  return _GiDfiles_EndMeshGroup( _GiDfiles_GetMeshFile( ) );
}

int GiD_fEndMeshGroup(GiD_FILE fd)
{
  CPostFile *File = NULL;
  
#ifdef HDF5
  if ( PostMode == GiD_PostHDF5 ) {
    return -1;
  }
#endif

  FD2FILE(fd,File);
  
  return _GiDfiles_EndMeshGroup( File );
}


/*
 *  Write a coordinate member at the current Coordinates Block 
 */

int GiD_WriteCoordinates( int id, double x, double y, double z )
{
  int res = 0;
#ifdef HDF5
  if(PostMode==GiD_PostHDF5){
    res = GiD_WriteCoordinates_HDF5(id,x,y,z);
    return res;
  }
#endif

  return _GiDfiles_WriteCoordinates( _GiDfiles_GetMeshFile(), id, x, y, z );
}

int GiD_fWriteCoordinates(GiD_FILE fd, int id, double x, double y, double z)
{
  CPostFile *File = NULL;
  int res = 0;
  
#ifdef HDF5
  if ( PostMode == GiD_PostHDF5 ) {
    return -1;
  }
#endif

  FD2FILE( fd, File );
  
  return _GiDfiles_WriteCoordinates( File, id, x, y, z );
}

int GiD_WriteCoordinates2D(int id, double x, double y)
{
#ifdef HDF5
  if(PostMode==GiD_PostHDF5){
    int res = GiD_WriteCoordinates2D_HDF5(id,x,y);
    return res;
  }
#endif

  return _GiDfiles_WriteCoordinates2D( _GiDfiles_GetMeshFile(), id, x, y );
}

int GiD_fWriteCoordinates2D(GiD_FILE fd, int id, double x, double y)
{
  CPostFile *File = NULL;
  
#ifdef HDF5
  if ( PostMode == GiD_PostHDF5 ) {
    return -1;
  }
#endif

  FD2FILE( fd, File );
  
  return _GiDfiles_WriteCoordinates2D( File, id, x, y );
}

/*
 *  Start a elements block in the current mesh
 */

int GiD_BeginElements()
{
#ifdef HDF5
  if(PostMode==GiD_PostHDF5){
    return GiD_BeginElements_HDF5();
  }
#endif

  return _GiDfiles_BeginElements( _GiDfiles_GetMeshFile() );
}

int GiD_fBeginElements(GiD_FILE fd)
{
  CPostFile *File = NULL;
  
#ifdef HDF5
  if ( PostMode == GiD_PostHDF5 ) {
    return -1;
  }
#endif

  FD2FILE( fd, File );
  return _GiDfiles_BeginElements( File );
}

/*
 *  Close the current elements block
 */

int GiD_EndElements()
{
#ifdef HDF5
  if(PostMode==GiD_PostHDF5){
    return GiD_EndElements_HDF5();
  }
#endif

  return _GiDfiles_EndElements( _GiDfiles_GetMeshFile() );
}

int GiD_fEndElements(GiD_FILE fd)
{
  CPostFile *File = NULL;
    
#ifdef HDF5
  if ( PostMode == GiD_PostHDF5 ) {
    return -1;
  }
#endif

  FD2FILE(fd,File);
  
  return _GiDfiles_EndElements(File);
}

/*
 *  Write an element member at the current Elements Block. The rest of
 *  the arguments define the connectivity of the element and its
 *  number depends on the NNode parameter given in the previous call
 *  to BeginPostMesh.
 *  
 */

int GiD_WriteElement( int id, int nid[] )
{
#ifdef HDF5
  if(PostMode==GiD_PostHDF5){
    return GiD_WriteElement_HDF5(id,nid);
  }
#endif

  return _GiDfiles_WriteElement( _GiDfiles_GetMeshFile(), id, nid);
}

int GiD_fWriteElement(GiD_FILE fd, int id, int nid[])
{
  CPostFile *File = NULL;
  
#ifdef HDF5
  if ( PostMode == GiD_PostHDF5 ) {
    return -1;
  }
#endif

  FD2FILE(fd,File);

  return _GiDfiles_WriteElement(File, id, nid);
}

/*
 *  Write an element member at the current Elements Block. The rest of
 *  the arguments define the connectivity of the element and its
 *  number depends on the NNode parameter given in the previous call
 *  to BeginPostMesh. The last id correspond to a material number.
 *  
 */

int GiD_WriteElementMat(int id, int nid[])
{
#ifdef HDF5
  if(PostMode==GiD_PostHDF5){
    return GiD_WriteElementMat_HDF5(id,nid);
  }
#endif

  return _GiDfiles_WriteElementMat( _GiDfiles_GetMeshFile(), id, nid);
}

int GiD_fWriteElementMat(GiD_FILE fd, int id, int nid[])
{
  CPostFile *File = NULL;
  
#ifdef HDF5
  if ( PostMode == GiD_PostHDF5 ) {
    return -1;
  }
#endif

  FD2FILE(fd,File);

  return _GiDfiles_WriteElementMat(File, id, nid);
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

int GiD_WriteSphere(int id, int nid, double r)
{
#ifdef HDF5
  if(PostMode==GiD_PostHDF5){
    return GiD_WriteSphere_HDF5(id,nid,r);
  }
#endif

  return _GiDfiles_WriteSphere( _GiDfiles_GetMeshFile(), id, nid, r);
}

int GiD_fWriteSphere(GiD_FILE fd, int id, int nid, double r)
{
  CPostFile *File = NULL;
  
#ifdef HDF5
  if ( PostMode == GiD_PostHDF5 ) {
    return -1;
  }
#endif

  FD2FILE(fd,File);

  return _GiDfiles_WriteSphere(File, id, nid, r);
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

int GiD_WriteSphereMat(int id, int nid, double r, int mat)
{
#ifdef HDF5
  if(PostMode==GiD_PostHDF5){
    return GiD_WriteSphereMat_HDF5(id,nid,r,mat);
  }
#endif

  return _GiDfiles_WriteSphereMat( _GiDfiles_GetMeshFile(), id, nid, r, mat);
}

int GiD_fWriteSphereMat(GiD_FILE fd, int id, int nid, double r, int mat)
{
  CPostFile *File = NULL;
  
#ifdef HDF5
  if ( PostMode == GiD_PostHDF5 ) {
    return -1;
  }
#endif

  FD2FILE(fd,File);

  return _GiDfiles_WriteSphereMat(File, id, nid, r, mat);
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

int GiD_WriteCircle(int id, int nid, double r,
		    double nx, double ny, double nz)
{
#ifdef HDF5
  if(PostMode==GiD_PostHDF5){
    return GiD_WriteCircle_HDF5(id,nid,r,nx,ny,nz);
  }
#endif

  return _GiDfiles_WriteCircle( _GiDfiles_GetMeshFile() , id, nid, r, nx, ny, nz);
}

int GiD_fWriteCircle(GiD_FILE fd, int id, int nid, double r,
		     double nx, double ny, double nz)
{
  CPostFile *File = NULL;
  
#ifdef HDF5
  if ( PostMode == GiD_PostHDF5 ) {
    return -1;
  }
#endif

  FD2FILE(fd,File);

  return _GiDfiles_WriteCircle(File, id, nid, r, nx, ny, nz);
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

int GiD_WriteCircleMat(int id, int nid, double r,
		       double nx, double ny, double nz, int mat)
{
#ifdef HDF5
  if(PostMode==GiD_PostHDF5){
    return GiD_WriteCircleMat_HDF5(id,nid,r,nx,ny,nz,mat);
  }
#endif
  
  return _GiDfiles_WriteCircleMat( _GiDfiles_GetMeshFile(), id, nid, r, nx, ny, nz, mat);
}

int GiD_fWriteCircleMat(GiD_FILE fd, int id, int nid, double r,
		        double nx, double ny, double nz, int mat)
{
  CPostFile *File = NULL;
  
#ifdef HDF5
  if ( PostMode == GiD_PostHDF5 ) {
    return -1;
  }
#endif

  FD2FILE(fd,File);

  return _GiDfiles_WriteCircleMat(File, id, nid, r, nx, ny, nz, mat);
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

int GiD_OpenPostResultFile( GP_CONST char * FileName, GiD_PostMode Mode )
{
  PostMode=Mode;

  if(PostMode==GiD_PostHDF5){
#ifdef HDF5
    return GiD_OpenPostResultFile_HDF5(FileName);
#else
    return 5;
#endif
  }

  if ( ResultFile ) 
    {
    /* must be closed */
    return  GP_ERROR_FILEOPENED;
    }

  if ( !( ResultFile = _GiDfiles_NewFile( Mode) ) )
    {
    /* not enough memory */
    return GP_ERROR_NOMEM;
    }

  if ( CPostFile_Open( ResultFile, FileName ) )
    {
    /* could not open file */
    CPostFile_Release( ResultFile );
    return GP_ERROR_OPENFAILED;
    }

  if ( CPostFile_WritePostHeader( ResultFile ) ) 
    {
    /* WritePostHeader failed */
    GiD_ClosePostResultFile();
    return GP_ERROR_WRITESTRING;
    }
  CPostFile_PushState( ResultFile, POST_S0 );
  return 0;
}

GiD_FILE GiD_fOpenPostResultFile(GP_CONST char * FileName, GiD_PostMode Mode)
{
  CPostFile *File = NULL;
  GiD_FILE fd;
  
#ifdef HDF5
  if ( PostMode == GiD_PostHDF5 ) {
    return -1;
  }
#endif

  if ( !( File = _GiDfiles_NewFile( Mode ) ) )
    {
    /* not enough memory = -2 GP_ERROR_NOMEM */
    return 0;
    }

  /* now open the File */
  if ( CPostFile_Open( File, FileName ) ) 
    {
    /* Open failed = -4 GP_ERROR_OPENFAILED */
    return 0;
    }
  fd = GiD_HashAdd( File );
  if ( !fd ) 
    {
    /* could not create a file handler = -5 GP_ERROR_HANDLEFAIL */
    CPostFile_Release( File );
    return 0;
    }

  if ( CPostFile_WritePostHeader( File ) )
    {
    /* WritePostHeader failed = -6 GP_ERROR_WRITESTRING */
    GiD_ClosePostResultFile( );
    return 0;
    } 
  else 
    {
    CPostFile_PushState( File, POST_S0 );
    }
  return fd;
}

/*
 *  Close the current post result file
 */

int GiD_ClosePostResultFile( )
{
  int fail;

#ifdef HDF5
  if(PostMode==GiD_PostHDF5){
    return GiD_ClosePostResultFile_HDF5();
  }
#endif

  assert( ResultFile != NULL );
  assert( _GiDfiles_CheckState( POST_S0, ResultFile ) );    

  if ( ResultFile ) 
    {
    fail = CPostFile_Release(ResultFile);
    ResultFile = NULL;
    /* reset outputMesh pointer */
    _GiDfiles_GetMeshFile();
    return fail;
    }
  return GP_ERROR_NULLFILE;
}

int GiD_fClosePostResultFile(GiD_FILE fd)
{
  int fail = 1;
  CPostFile *File = NULL;
  
#ifdef HDF5
  if ( PostMode == GiD_PostHDF5 ) {
    return -1;
  }
#endif

  FD2FILE( fd, File );
  
  assert( _GiDfiles_CheckState( POST_S0, File ) );    
  fail = CPostFile_Release( File );
  GiD_HashRemove( fd );
  
  return fail;
}

/*
 *  Begin Gauss Points definition
 */

int GiD_BeginGaussPoint( GP_CONST char * name, GiD_ElementType EType,
                         GP_CONST char * MeshName,
                         int GP_number, int NodesIncluded, int InternalCoord )
{
#ifdef HDF5
  if(PostMode==GiD_PostHDF5){
    return GiD_BeginGaussPoint_HDF5(name,EType,MeshName,GP_number,NodesIncluded,InternalCoord);
  }
#endif
  
  return _GiDfiles_BeginGaussPoint( ResultFile,
                               name, EType, MeshName, GP_number,
                               NodesIncluded, InternalCoord );
}

int GiD_fBeginGaussPoint( GiD_FILE fd, GP_CONST char * name,
                          GiD_ElementType EType,
                          GP_CONST char * MeshName,
                          int GP_number, int NodesIncluded, int InternalCoord )
{
  CPostFile *File = NULL;
  
#ifdef HDF5
  if ( PostMode == GiD_PostHDF5 ) {
    return -1;
  }
#endif

  FD2FILE( fd, File );
  
  return _GiDfiles_BeginGaussPoint( File,
                               name, EType, MeshName, GP_number,
                               NodesIncluded, InternalCoord );
}


/*
 *  End current Gauss Points definition
 */

int GiD_EndGaussPoint( )
{
#ifdef HDF5
  if(PostMode==GiD_PostHDF5){
    return GiD_EndGaussPoint_HDF5( );
  }
#endif

  return _GiDfiles_EndGaussPoint( ResultFile );
}

int GiD_fEndGaussPoint( GiD_FILE fd )
{
  CPostFile *File = NULL;
  
#ifdef HDF5
  if ( PostMode == GiD_PostHDF5 ) {
    return -1;
  }
#endif

  FD2FILE( fd,File );

  return _GiDfiles_EndGaussPoint( File );
}

/*
 *  Write internal gauss point coordinate.
 */

int GiD_WriteGaussPoint2D(double x, double y)
{
#ifdef HDF5
  if(PostMode==GiD_PostHDF5){
    return GiD_WriteGaussPoint2D_HDF5(x,y);
  }
#endif

  return _GiDfiles_WriteGaussPoint2D( ResultFile, x, y );
}

int GiD_fWriteGaussPoint2D(GiD_FILE fd, double x, double y)
{
  CPostFile *File = NULL;
  
#ifdef HDF5
  if ( PostMode == GiD_PostHDF5 ) {
    return -1;
  }
#endif

  FD2FILE( fd, File );

  return _GiDfiles_WriteGaussPoint2D( File, x, y );
}

int GiD_WriteGaussPoint3D( double x, double y, double z )
{
#ifdef HDF5
  if(PostMode==GiD_PostHDF5){
    return GiD_WriteGaussPoint3D_HDF5(x,y,z);
  }
#endif
  
  return _GiDfiles_WriteGaussPoint3D( ResultFile, x, y, z );
}

int GiD_fWriteGaussPoint3D(GiD_FILE fd, double x, double y, double z)
{
  CPostFile *File = NULL;
  
#ifdef HDF5
  if ( PostMode == GiD_PostHDF5 ) {
    return -1;
  }
#endif

  FD2FILE(fd,File);

  return _GiDfiles_WriteGaussPoint3D( File, x, y, z );
}

/*
 *  Begin a Range Table definition
 */

int GiD_BeginRangeTable( GP_CONST char * name )
{
#ifdef HDF5
  if(PostMode==GiD_PostHDF5){
    return GiD_BeginRangeTable_HDF5(name);
  }
#endif
  
  return _GiDfiles_BeginRangeTable( ResultFile, name ); 
}

int GiD_fBeginRangeTable( GiD_FILE fd, GP_CONST char * name )
{
  CPostFile *File = NULL;
  
#ifdef HDF5
  if ( PostMode == GiD_PostHDF5 ) {
    return -1;
  }
#endif

  FD2FILE( fd, File );
  
  return _GiDfiles_BeginRangeTable( File, name );
}

/*
 *  End a Range Table definition
 */

int GiD_EndRangeTable()
{
#ifdef HDF5
  if(PostMode==GiD_PostHDF5){
    return GiD_EndRangeTable_HDF5();
  }
#endif

  return _GiDfiles_EndRangeTable( ResultFile );
}

int GiD_fEndRangeTable(GiD_FILE fd)
{
  CPostFile *File = NULL;
  
#ifdef HDF5
  if ( PostMode == GiD_PostHDF5 ) {
    return -1;
  }
#endif

  FD2FILE( fd, File );

  return _GiDfiles_EndRangeTable( File );
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

int GiD_WriteMinRange( double max, GP_CONST char * name )
{
#ifdef HDF5
  if(PostMode==GiD_PostHDF5){
    return GiD_WriteMinRange_HDF5(max,name);
  }
#endif

  return _GiDfiles_WriteMinRange( ResultFile, max, name );
}

int GiD_fWriteMinRange( GiD_FILE fd, double max, GP_CONST char * name )
{
  CPostFile *File = NULL;
    
#ifdef HDF5
  if ( PostMode == GiD_PostHDF5 ) {
    return -1;
  }
#endif

  FD2FILE( fd, File );

  return _GiDfiles_WriteMinRange( File, max, name );
}

int GiD_WriteRange( double min, double max, GP_CONST char * name )
{
#ifdef HDF5
  if(PostMode==GiD_PostHDF5){
    return GiD_WriteRange_HDF5(min,max,name);
  }
#endif
  
  return _GiDfiles_WriteRange( ResultFile, min, max, name );
}

int GiD_fWriteRange( GiD_FILE fd, double min, double max, GP_CONST char * name )
{
  CPostFile *File = NULL;
  
#ifdef HDF5
  if ( PostMode == GiD_PostHDF5 ) {
    return -1;
  }
#endif

  FD2FILE( fd, File );

  return _GiDfiles_WriteRange(File, min, max, name);
}

int GiD_WriteMaxRange( double min, GP_CONST char * name )
{
#ifdef HDF5
  if(PostMode==GiD_PostHDF5){
    return GiD_WriteMaxRange_HDF5(min,name);
  }
#endif

  return _GiDfiles_WriteMaxRange(ResultFile, min, name);
}

int GiD_fWriteMaxRange(GiD_FILE fd, double min, GP_CONST char * name)
{
  CPostFile *File = NULL;
  
#ifdef HDF5
  if ( PostMode == GiD_PostHDF5 ) {
    return -1;
  }
#endif

  FD2FILE( fd, File );
  
  return _GiDfiles_WriteMaxRange(File, min, name);
}

/*
 *  Begin Result Block
 */

int GiD_BeginResult( GP_CONST char     *Result,
                     GP_CONST char     *Analysis,
                     double             step,
                     GiD_ResultType     Type,
                     GiD_ResultLocation Where,
                     GP_CONST char     *GaussPointsName,
                     GP_CONST char     *RangeTable, 
                     int                compc,
                     GP_CONST char     *compv[] )
{
#ifdef HDF5
  if(PostMode==GiD_PostHDF5){
    return GiD_BeginResult_HDF5(Result,Analysis,step,Type,Where,GaussPointsName,RangeTable,compc,compv);
  }
#endif
  
  return _GiDfiles_BeginResult( ResultFile, Result, Analysis, step, Type, Where,
                           GaussPointsName, RangeTable, compc, compv );
}

int GiD_fBeginResult( GiD_FILE fd, GP_CONST char * Result,
                      GP_CONST char * Analysis,
                      double step,
                      GiD_ResultType Type, GiD_ResultLocation Where,
                      GP_CONST char * GaussPointsName,
                      GP_CONST char * RangeTable, 
		     int compc, GP_CONST char * compv[] )
{
  CPostFile *File = NULL;
  
#ifdef HDF5
  if ( PostMode == GiD_PostHDF5 ) {
    return -1;
  }
#endif

  FD2FILE( fd, File );

  return _GiDfiles_BeginResult( File, Result, Analysis, step, Type, Where,
                           GaussPointsName, RangeTable, compc, compv );
}

int GiD_BeginResultHeader( GP_CONST char     *Result,
                           GP_CONST char     *Analysis,
                           double             step,
                           GiD_ResultType     Type,
                           GiD_ResultLocation Where,
                           GP_CONST char     *GaussPointsName )
{
#ifdef HDF5
  if(PostMode==GiD_PostHDF5){
    return GiD_BeginResultHeader_HDF5(Result,Analysis,step,Type,Where,GaussPointsName);
  }
#endif
  
  return _GiDfiles_BeginResultHeader( ResultFile, Result, Analysis, step, Type,
                                 Where, GaussPointsName );
}

int GiD_fBeginResultHeader( GiD_FILE fd, GP_CONST char * Result,
                            GP_CONST char * Analysis, 
                            double step,
                            GiD_ResultType Type, GiD_ResultLocation Where,
                            GP_CONST char * GaussPointsName )
{
  CPostFile *File = NULL;
  
#ifdef HDF5
  if ( PostMode == GiD_PostHDF5 ) {
    return -1;
  }
#endif

  FD2FILE( fd, File );
  
  return _GiDfiles_BeginResultHeader( File, Result, Analysis, step, Type,
                                 Where, GaussPointsName );
}

int GiD_ResultRange( GP_CONST char * RangeTable )
{
#ifdef HDF5
  if(PostMode==GiD_PostHDF5){
    return GiD_ResultRange_HDF5(RangeTable);
  }
#endif

  return _GiDfiles_ResultRange( ResultFile, RangeTable );
}

int GiD_fResultRange( GiD_FILE fd, GP_CONST char * RangeTable )
{
  CPostFile *File = NULL;
  
#ifdef HDF5
  if ( PostMode == GiD_PostHDF5 ) {
    return -1;
  }
#endif

  FD2FILE( fd, File );

  return _GiDfiles_ResultRange( File, RangeTable );
}

int GiD_ResultComponents( int compc, GP_CONST char * compv[] )
{
#ifdef HDF5
  if(PostMode==GiD_PostHDF5){
    return GiD_ResultComponents_HDF5(compc,compv);
  }
#endif

  return _GiDfiles_ResultComponents( ResultFile, compc, compv );
}

int GiD_fResultComponents( GiD_FILE fd, int compc, GP_CONST char * compv[] )
{
  CPostFile *File = NULL;
  
#ifdef HDF5
  if ( PostMode == GiD_PostHDF5 ) {
    return -1;
  }
#endif

  FD2FILE( fd, File );
  
  return _GiDfiles_ResultComponents( File, compc, compv );
}

int GiD_ResultUnit( GP_CONST char * UnitName )
{
#ifdef HDF5
  if(PostMode==GiD_PostHDF5){
    return GiD_ResultUnit_HDF5(UnitName);
  }
#endif
  return _GiDfiles_ResultUnit( ResultFile, UnitName );  
}

int GiD_fResultUnit(GiD_FILE fd, GP_CONST char * UnitName)
{
  CPostFile *File = NULL;  
#ifdef HDF5
  if ( PostMode == GiD_PostHDF5 ) {
    return -1;
  }
#endif

  FD2FILE( fd, File );
  return _GiDfiles_ResultUnit( File, UnitName );  
}

int GiD_ResultUserDefined( GP_CONST char * Name,GP_CONST char * Value )
{
#ifdef HDF5
  if(PostMode==GiD_PostHDF5){
    return GiD_ResultUserDefined_HDF5(Name,Value);
  }
#endif
  return _GiDfiles_ResultUserDefined( ResultFile, Name, Value );
}

int GiD_fResultUserDefined( GiD_FILE fd, GP_CONST char * Name,GP_CONST char * Value )
{
  CPostFile *File = NULL;  
#ifdef HDF5
  if ( PostMode == GiD_PostHDF5 ) {
    return -1;
  }
#endif

  FD2FILE( fd, File );
  return _GiDfiles_ResultUserDefined( File, Name, Value );
}

int GiD_BeginResultGroup( GP_CONST char     *Analysis,
                          double             step,
                          GiD_ResultLocation Where,
                          GP_CONST char     *GaussPointsName )
{
#ifdef HDF5
  if(PostMode==GiD_PostHDF5){
    return GiD_BeginResultGroup_HDF5(Analysis,step,Where,GaussPointsName);
  }
#endif

  return _GiDfiles_BeginResultGroup( ResultFile, Analysis, step, Where,
                                GaussPointsName );
}

int GiD_fBeginResultGroup( GiD_FILE fd, GP_CONST char * Analysis, double step,
		           GiD_ResultLocation Where,
                           GP_CONST char * GaussPointsName )
{
  CPostFile *File = NULL;
  
#ifdef HDF5
  if ( PostMode == GiD_PostHDF5 ) {
    return -1;
  }
#endif

  FD2FILE( fd, File );

  return _GiDfiles_BeginResultGroup( File, Analysis, step, Where,
                                GaussPointsName );
}

int GiD_ResultDescription( GP_CONST char * Result, GiD_ResultType Type )
{
#ifdef HDF5
  if(PostMode==GiD_PostHDF5){
    return GiD_ResultDescription_HDF5(Result,Type);
  }
#endif

  return _GiDfiles_ResultDescription_( ResultFile, Result, Type, 0 );
}

int GiD_fResultDescription( GiD_FILE fd,
		            GP_CONST char * Result, GiD_ResultType Type )
{
  CPostFile *File = NULL;
  
#ifdef HDF5
  if ( PostMode == GiD_PostHDF5 ) {
    return -1;
  }
#endif

  FD2FILE( fd, File );
  
  return _GiDfiles_ResultDescription_( File, Result, Type, 0 );
}

int GiD_ResultDescriptionDim( GP_CONST char * Result, GiD_ResultType Type,
                              size_t s )
{
#ifdef HDF5
  if(PostMode==GiD_PostHDF5){
    /* disregard s in hdf5, because ResultGroup is really written in HDF5 as independent Results*/
    return GiD_ResultDescription_HDF5(Result,Type);
  }
#endif

  return _GiDfiles_ResultDescription_( ResultFile, Result, Type, s );
}

int GiD_fResultDescriptionDim( GiD_FILE fd,
                               GP_CONST char * Result, GiD_ResultType Type,
                               size_t s )
{
  CPostFile *File = NULL;
  
#ifdef HDF5
  if ( PostMode == GiD_PostHDF5 ) {
    return -1;
  }
#endif

  FD2FILE( fd, File );
  
  return _GiDfiles_ResultDescription_( File, Result, Type, s );
}

int GiD_ResultLocalAxes( GP_CONST char * Result, GP_CONST char * Analysis,
                         double step,double vx,double vy,double vz )
{
#ifdef HDF5
    if(PostMode==GiD_PostHDF5){
    return GiD_ResultLocalAxes_HDF5(Result,Analysis,step,vx,vy,vz);
  }
#endif
  return 1;
}

int GiD_ResultIsDeformationVector( int boolean )
{
#ifdef HDF5
    if(PostMode==GiD_PostHDF5){
    return GiD_ResultIsDeformationVector_HDF5(boolean);
  }
#endif
  return 1;
}

int GiD_ResultValues()
{
#ifdef HDF5
  if(PostMode==GiD_PostHDF5){
    return GiD_ResultValues_HDF5();
  }
#endif

  return GP_OK;
}

int GiD_fResultValues( GiD_FILE fd )
{  
#ifdef HDF5
  if ( PostMode == GiD_PostHDF5 ) {
    return -1;
  }
#endif

  return GP_OK;
}

/*
 *  End Result Block
 */

int GiD_EndResult( )
{
#ifdef HDF5
  if(PostMode==GiD_PostHDF5){
    return GiD_EndResult_HDF5();
  }
#endif

  return _GiDfiles_EndResult( ResultFile );
}

int GiD_fEndResult( GiD_FILE fd )
{
  CPostFile *File = NULL;
  
#ifdef HDF5
  if ( PostMode == GiD_PostHDF5 ) {
    return -1;
  }
#endif

  FD2FILE( fd, File );

  return _GiDfiles_EndResult( File );
}

/*
 * Results which are referred to a mesh group (see GiD_BeginMeshGroup)
 * should be written between a call to this function and
 * GiD_EndOnMeshGroup.
 */

int GiD_BeginOnMeshGroup( char * Name )
{
#ifdef HDF5
  if(PostMode==GiD_PostHDF5){
    return GiD_BeginOnMeshGroup_HDF5( Name);
  }
#endif

  return _GiDfiles_BeginOnMeshGroup( ResultFile, Name );
}

int GiD_fBeginOnMeshGroup(GiD_FILE fd, char * Name)
{
  CPostFile *File = NULL;
  
#ifdef HDF5
  if ( PostMode == GiD_PostHDF5 ) {
    return -1;
  }
#endif

  FD2FILE( fd, File );

  return _GiDfiles_BeginOnMeshGroup( File, Name );
}

/*
 * This function close a previously opened block of result over a mesh
 * group.
 */

int GiD_EndOnMeshGroup()
{
#ifdef HDF5
  if(PostMode==GiD_PostHDF5){
    return GiD_EndOnMeshGroup_HDF5();
  }
#endif
  return _GiDfiles_EndOnMeshGroup( ResultFile );
}

int GiD_fEndOnMeshGroup(GiD_FILE fd)
{
  CPostFile *File = NULL;
  
#ifdef HDF5
  if ( PostMode == GiD_PostHDF5 ) {
    return -1;
  }
#endif

  FD2FILE( fd, File );
  
  return _GiDfiles_EndOnMeshGroup( File );
}

/*
 * Flushes all pending output into the compressed file.
 */

int GiD_FlushPostFile()
{
  int res1 = 0, res2 = 0;
  
#ifdef HDF5
  if( PostMode == GiD_PostHDF5 ) {
    return GiD_FlushPostFile_HDF5();
  }
#endif
  if( MeshFile ) res1 = CPostFile_Flush( MeshFile );
  if( !res1 && ResultFile ) res2 = CPostFile_Flush( ResultFile );
  return res1 || res2;
}

int GiD_fFlushPostFile(GiD_FILE fd)
{
  CPostFile *File = NULL;
  
#ifdef HDF5
  if ( PostMode == GiD_PostHDF5 ) {
    return -1;
  }
#endif

  FD2FILE( fd, File );

  return CPostFile_Flush( File );
}

/*
 *  Write result functions
 */

static
int GiD_EnsureBeginValues( CPostFile *File )
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

int GiD_WriteScalar(int id, double v)
{
#ifdef HDF5
  if(PostMode==GiD_PostHDF5){
    return GiD_WriteScalar_HDF5(id,v);
  }
#endif

  return _GiDfiles_WriteScalar( ResultFile, id, v );
}

int GiD_fWriteScalar(GiD_FILE fd, int id, double v)
{
  CPostFile *File = NULL;
  
#ifdef HDF5
  if ( PostMode == GiD_PostHDF5 ) {
    return -1;
  }
#endif

  FD2FILE( fd, File );
  
  return _GiDfiles_WriteScalar( File, id, v );
}

int GiD_Write2DVector( int id, double x, double y )
{
#ifdef HDF5
  if(PostMode==GiD_PostHDF5){
    return GiD_Write2DVector_HDF5( id, x, y );
  }
#endif

  return _GiDfiles_Write2DVector( ResultFile, id, x, y );
}

int GiD_fWrite2DVector( GiD_FILE fd, int id, double x, double y )
{
  CPostFile *File = NULL;
  
#ifdef HDF5
  if ( PostMode == GiD_PostHDF5 ) {
    return -1;
  }
#endif

  FD2FILE( fd, File );
  
  return _GiDfiles_Write2DVector( File, id, x, y );
}

int GiD_WriteVector( int id, double x, double y, double z )
{
#ifdef HDF5
  if(PostMode==GiD_PostHDF5){
    return GiD_WriteVector_HDF5(id,x,y,z);
  }
#endif
  
  return _GiDfiles_WriteVector( ResultFile, id, x, y, z );
}

int GiD_fWriteVector( GiD_FILE fd, int id, double x, double y, double z )
{
  CPostFile *File = NULL;
  
#ifdef HDF5
  if ( PostMode == GiD_PostHDF5 ) {
    return -1;
  }
#endif

  FD2FILE( fd, File );
  
  return _GiDfiles_WriteVector( File, id, x, y, z );
}

int GiD_WriteVectorModule(int id, double x, double y, double z, double mod)
{
#ifdef HDF5
  if(PostMode==GiD_PostHDF5){
    return GiD_WriteVectorModule_HDF5(id,x,y,z,mod);
  }
#endif

  return _GiDfiles_WriteVectorModule(ResultFile, id, x, y, z, mod);
}

int GiD_fWriteVectorModule(GiD_FILE fd,
		           int id, double x, double y, double z, double mod)
{
  CPostFile *File = NULL;
  
#ifdef HDF5
  if ( PostMode == GiD_PostHDF5 ) {
    return -1;
  }
#endif

  FD2FILE( fd, File );
  
  return _GiDfiles_WriteVectorModule(File, id, x, y, z, mod);
}

int GiD_Write2DMatrix(int id, double Sxx, double Syy, double Sxy)
{
#ifdef HDF5
  if(PostMode==GiD_PostHDF5){
    return GiD_Write2DMatrix_HDF5(id,Sxx,Syy,Sxy);
  }
#endif

  return _GiDfiles_Write2DMatrix(ResultFile, id, Sxx, Syy, Sxy);
}

int GiD_fWrite2DMatrix(GiD_FILE fd, int id, double Sxx, double Syy, double Sxy)
{
  CPostFile *File = NULL;
  
#ifdef HDF5
  if ( PostMode == GiD_PostHDF5 ) {
    return -1;
  }
#endif

  FD2FILE( fd, File );
  
  return _GiDfiles_Write2DMatrix(File, id, Sxx, Syy, Sxy);
}

int GiD_Write3DMatrix(int id,
		      double Sxx, double Syy, double Szz,
		      double Sxy, double Syz, double Sxz)
{
#ifdef HDF5
  if(PostMode==GiD_PostHDF5){
    return GiD_Write3DMatrix_HDF5(id,Sxx,Syy,Szz,Sxy,Syz,Sxz);
  }
#endif

  return _GiDfiles_Write3DMatrix(ResultFile, id,
		            Sxx, Syy, Szz,
		            Sxy, Syz, Sxz);
}

int GiD_fWrite3DMatrix(GiD_FILE fd, int id,
		       double Sxx, double Syy, double Szz,
		       double Sxy, double Syz, double Sxz)
{
  CPostFile *File = NULL;
  
#ifdef HDF5
  if ( PostMode == GiD_PostHDF5 ) {
    return -1;
  }
#endif

  FD2FILE( fd, File );
  
  return _GiDfiles_Write3DMatrix(File, id,
		            Sxx, Syy, Szz,
		            Sxy, Syz, Sxz);
}

int GiD_WritePlainDefMatrix(int id,
		            double Sxx, double Syy, double Sxy, double Szz )
{
#ifdef HDF5
  if(PostMode==GiD_PostHDF5){
    return GiD_WritePlainDefMatrix_HDF5(id,Sxx,Syy,Sxy,Szz);
  }
#endif
  
  return _GiDfiles_WritePlainDefMatrix(ResultFile,
		                  id, Sxx, Syy, Sxy, Szz);
}

int GiD_fWritePlainDefMatrix(GiD_FILE fd, int id,
		             double Sxx, double Syy, double Sxy, double Szz )
{
  CPostFile *File = NULL;
  
#ifdef HDF5
  if ( PostMode == GiD_PostHDF5 ) {
    return -1;
  }
#endif

  FD2FILE( fd, File );
  
  return _GiDfiles_WritePlainDefMatrix(File,
		                  id, Sxx, Syy, Sxy, Szz);
}

int GiD_WriteMainMatrix(int id,
		        double Si, double Sii, double Siii,
		        double Vix, double Viy, double Viz,
		        double Viix, double Viiy, double Viiz,
		        double Viiix, double Viiiy, double Viiiz)
{
#ifdef HDF5
  if(PostMode==GiD_PostHDF5){
    return GiD_WriteMainMatrix_HDF5(id,Si,Sii,Siii,Vix,Viy,Viz,Viix,Viiy,Viiz,Viiix,Viiiy,Viiiz);
  }
#endif

  return _GiDfiles_WriteMainMatrix(ResultFile, id, 
		              Si,  Sii,  Siii,
		              Vix,  Viy,  Viz,
		              Viix,  Viiy,  Viiz,
		              Viiix,  Viiiy,  Viiiz);
}

int GiD_fWriteMainMatrix(GiD_FILE fd, int id,
		         double Si, double Sii, double Siii,
		         double Vix, double Viy, double Viz,
		         double Viix, double Viiy, double Viiz,
		         double Viiix, double Viiiy, double Viiiz)
{
  CPostFile *File = NULL;
  
#ifdef HDF5
  if ( PostMode == GiD_PostHDF5 ) {
    return -1;
  }
#endif

  FD2FILE( fd, File );
  
  return _GiDfiles_WriteMainMatrix(File, id, 
		              Si,  Sii,  Siii,
		              Vix,  Viy,  Viz,
		              Viix,  Viiy,  Viiz,
		              Viiix,  Viiiy,  Viiiz);
}

int GiD_WriteLocalAxes(int id, double euler_1, double euler_2, double euler_3)
{
#ifdef HDF5
  if(PostMode==GiD_PostHDF5){
    return GiD_WriteLocalAxes_HDF5(id,euler_1,euler_2,euler_3);
  }
#endif

  return _GiDfiles_WriteLocalAxes(ResultFile, id, euler_1, euler_2, euler_3);
}

int GiD_fWriteLocalAxes(GiD_FILE fd,
		        int id, double euler_1, double euler_2, double euler_3)
{
  CPostFile *File = NULL;
  
#ifdef HDF5
  if ( PostMode == GiD_PostHDF5 ) {
    return -1;
  }
#endif

  FD2FILE( fd, File );
  
  return _GiDfiles_WriteLocalAxes(File, id, euler_1, euler_2, euler_3);
}

/*
 * Complex numbers
 */

int GiD_WriteComplexScalar( int id, double complex_real, double complex_imag) {
#ifdef HDF5
  if ( PostMode == GiD_PostHDF5) {
    return GiD_WriteComplexScalar_HDF5(id, complex_real, complex_imag);
  }
#endif
  return _GiDfiles_WriteComplexScalar(ResultFile, id, complex_real, complex_imag);
}

int GiD_fWriteComplexScalar(GiD_FILE fd, int id, double complex_real, double complex_imag) {
  CPostFile *File = NULL;  
#ifdef HDF5
  if ( PostMode == GiD_PostHDF5 ) {
    return -1;
  }
#endif

  FD2FILE( fd, File);
  return _GiDfiles_WriteComplexScalar(File, id, complex_real, complex_imag);
}

int GiD_Write2DComplexVector( int id, 
			      double x_real, double x_imag,
			      double y_real, double y_imag) {
#ifdef HDF5
  if ( PostMode == GiD_PostHDF5) {
    return GiD_Write2DComplexVector_HDF5(id, x_real, x_imag, y_real, y_imag);
  }
#endif
  return _GiDfiles_Write2DComplexVector(ResultFile, id, x_real, x_imag, y_real, y_imag);
}

int GiD_fWrite2DComplexVector( GiD_FILE fd, int id,
			       double x_real, double x_imag,
			       double y_real, double y_imag) {
  CPostFile *File = NULL;  
#ifdef HDF5
  if ( PostMode == GiD_PostHDF5 ) {
    return -1;
  }
#endif

  FD2FILE( fd, File);
  return _GiDfiles_Write2DComplexVector(File, id, x_real, x_imag, y_real, y_imag);
}

int GiD_WriteComplexVector( int id, 
			    double x_real, double x_imag,
			    double y_real, double y_imag,
			    double z_real, double z_imag) {
#ifdef HDF5
  if ( PostMode == GiD_PostHDF5) {
    return GiD_WriteComplexVector_HDF5(id, x_real, x_imag, y_real, y_imag, z_real, z_imag);
  }
#endif
  return _GiDfiles_WriteComplexVector(ResultFile, id, x_real, x_imag, y_real, y_imag, z_real, z_imag);
}

int GiD_fWriteComplexVector( GiD_FILE fd, int id,
			     double x_real, double x_imag,
			     double y_real, double y_imag,
			     double z_real, double z_imag) {
  CPostFile *File = NULL;  
#ifdef HDF5
  if ( PostMode == GiD_PostHDF5 ) {
    return -1;
  }
#endif

  FD2FILE( fd, File);
  return _GiDfiles_WriteComplexVector(File, id, x_real, x_imag, y_real, y_imag, z_real, z_imag);
}


int GiD_Write2DComplexMatrix(int id, double Sxx_real, double Syy_real, double Sxy_real,
                             double Sxx_imag, double Syy_imag, double Sxy_imag) {
#ifdef HDF5
  if ( PostMode == GiD_PostHDF5) {
    return GiD_Write2DComplexMatrix( id, Sxx_real, Syy_real, Sxy_real,
                                     Sxx_imag, Syy_imag, Sxy_imag);
  }
#endif
  return _GiDfiles_Write2DComplexMatrix(ResultFile, id, Sxx_real, Syy_real, Sxy_real,
                                        Sxx_imag, Syy_imag, Sxy_imag);
}

int GiD_fWrite2DComplexMatrix(GiD_FILE fd, int id,
                              double Sxx_real, double Syy_real, double Sxy_real,
                              double Sxx_imag, double Syy_imag, double Sxy_imag) {
  CPostFile *File = NULL;  
#ifdef HDF5
  if ( PostMode == GiD_PostHDF5 ) {
    return -1;
  }
#endif

  FD2FILE( fd, File);
  return _GiDfiles_Write2DComplexMatrix(File, id, Sxx_real, Syy_real, Sxy_real,
                                        Sxx_imag, Syy_imag, Sxy_imag);
}

int GiD_Write3DComplexMatrix(int id,
                             double Sxx_real, double Syy_real, double Szz_real,
                             double Sxy_real, double Syz_real, double Sxz_real,
                             double Sxx_imag, double Syy_imag, double Szz_imag,
                             double Sxy_imag, double Syz_imag, double Sxz_imag) {
#ifdef HDF5
  if ( PostMode == GiD_PostHDF5) {
    return GiD_WriteComplexMatrix_HDF5(id, 
                                       Sxx_real, Syy_real, Szz_real,
                                       Sxy_real, Syz_real, Sxz_real,
                                       Sxx_imag, Syy_imag, Szz_imag,
                                       Sxy_imag, Syz_imag, Sxz_imag);
  }
#endif
  return _GiDfiles_WriteComplexMatrix(ResultFile, id, 
                                      Sxx_real, Syy_real, Szz_real,
                                      Sxy_real, Syz_real, Sxz_real,
                                      Sxx_imag, Syy_imag, Szz_imag,
                                      Sxy_imag, Syz_imag, Sxz_imag);
}

int GiD_fWrite3DComplexMatrix(GiD_FILE fd, int id,
                              double Sxx_real, double Syy_real, double Szz_real,
                              double Sxy_real, double Syz_real, double Sxz_real,
                              double Sxx_imag, double Syy_imag, double Szz_imag,
                              double Sxy_imag, double Syz_imag, double Sxz_imag) {
  CPostFile *File = NULL;  
#ifdef HDF5
  if ( PostMode == GiD_PostHDF5 ) {
    return -1;
  }
#endif

  FD2FILE( fd, File);
  return _GiDfiles_WriteComplexMatrix(File, id, 
                                      Sxx_real, Syy_real, Szz_real,
                                      Sxy_real, Syz_real, Sxz_real,
                                      Sxx_imag, Syy_imag, Szz_imag,
                                      Sxy_imag, Syz_imag, Sxz_imag);
}


/* 
* Nurbs 
*/

int GiD_WriteNurbsSurface( int id, int n, double* v )
{
#ifdef HDF5
  if(PostMode==GiD_PostHDF5){
    return GiD_WriteNurbsSurface_HDF5(id,n,v);
  }
#endif
  
  return _GiDfiles_WriteNurbsSurface( ResultFile, id, n, v );
}

int GiD_fWriteNurbsSurface( GiD_FILE fd, int id, int n, double* v )
{
  CPostFile *File = NULL;
  
#ifdef HDF5
  if ( PostMode == GiD_PostHDF5 ) {
    return -1;
  }
#endif

  FD2FILE( fd, File );
  
  return _GiDfiles_WriteNurbsSurface( File, id, n, v );
}

int GiD_WriteNurbsSurfaceVector( int id, int n, int num_comp, double* v )
{
#ifdef HDF5
  if(PostMode==GiD_PostHDF5){
    return GiD_WriteNurbsSurfaceVector_HDF5(id,n,num_comp,v);
  }
#endif
  
  return _GiDfiles_WriteNurbsSurfaceVector( ResultFile, id, n, num_comp, v );
}

int GiD_fWriteNurbsSurfaceVector( GiD_FILE fd, int id, int n, int num_comp, double* v )
{
  CPostFile *File = NULL;
  
#ifdef HDF5
  if ( PostMode == GiD_PostHDF5 ) {
    return -1;
  }
#endif

  FD2FILE( fd, File );
  
  return _GiDfiles_WriteNurbsSurfaceVector( File, id, n, num_comp, v );
}
