PROGRAM EXAMPLE
  use gidpost

  REAL(8) :: nodes(1:9,1:3), rx, ry, value
  INTEGER :: elems(1:4,1:4)
  INTEGER :: idx, x, y, nid(1:4)
  CHARACTER(LEN=1024) :: nombre
  CHARACTER(LEN=1024) :: component_names(2)

  idx = 1
  do x=0,2
    do y=0,2
      nodes(idx,1) = x
      nodes(idx,2) = y
      nodes(idx,3) = 0
      idx = idx + 1
    end do
  end do
  idx = 1
  do x=1,2
    do y=1,2
      elems(idx,1) = x + (y-1)*3
      elems(idx,2) = elems(idx,1) + 1
      elems(idx,4) = elems(idx,1) + 3
      elems(idx,3) = elems(idx,4) + 1
      idx = idx + 1
    end do
  end do

  !CALL GID_OPENPOSTMESHFILE('testfortran90.post.msh',GiD_PostAscii)      
  !CALL GiD_OpenPostResultFile('testfortran90.post.res',GiD_PostBinary)      
  CALL GiD_OpenPostResultFile('testfortran90.post.res',GiD_PostHDF5)
  CALL GiD_BeginMeshColor('quadmesh',GiD_2D,GiD_Quadrilateral,4,0.7d0,0.7d0,0.4d0)
  CALL GiD_MeshUnit('mm');
  CALL GiD_BeginCoordinates
  idx = 1
  do x=0,2
    do y=0,2
      rx = x
      ry = y
      CALL GiD_WriteCoordinates2D(idx, rx, ry)
      idx = idx + 1
    end do
  end do
  CALL GiD_EndCoordinates
  CALL GiD_BeginElements
  idx = 1

  do x=1,2
    do y=1,2
      nid(1) = x + (y-1)*3
      nid(2) = nid(1) + 1
      nid(4) = nid(1) + 3
      nid(3) = nid(4) + 1
      CALL GiD_WriteElement(idx, nid)
      idx = idx + 1
    end do
  end do      
  CALL GiD_EndElements
  CALL GiD_EndMesh
  !Scalar result  
  CALL GiD_BeginResultHeader('Pressure','Analysis',1.0d0,GiD_Scalar,GiD_onNodes,GiD_NULL)  
  !CALL GiD_ScalarComp('Comp')
  CALL GiD_ResultUnit('Pa');
  CALL GiD_ResultValues
  !CALL GiD_BeginScalarResult('Pressure','Analysis',1.0d0,GiD_onNodes,GiD_NULL,GiD_NULL,GiD_NULL) 
  do idx=1,9
    value=idx*1.5
    CALL GiD_WriteScalar(idx,value)
  end do
  CALL GiD_EndResult
  CALL GiD_FlushPostFile
  !Vectorial result
  CALL GiD_BeginResultHeader('Velocity','Analysis',1.0d0,GiD_Vector,GiD_onNodes,GiD_NULL)
  CALL GiD_VectorComp('Comp X','Comp Y','Comp Z','Module')
  CALL GiD_ResultUnit('m/s');
  CALL GiD_ResultValues
  do idx=1,9
    value=idx*1.5
    CALL GiD_Write2DVector(idx,value,value*3)
  end do
  CALL GiD_EndResult

  !Complex Scalar result
  CALL GiD_BeginResultHeader('Potential','Analysis',1.0d0,GiD_ComplexScalar,GiD_onNodes,GiD_NULL)    
  component_names(1)='real_part'
  component_names(2)='imaginary_part'
  CALL GiD_ResultComponents(2,component_names)
  CALL GiD_ResultUnit('Volt');
  CALL GiD_ResultValues
  !CALL GiD_BeginComplexScalarResult('Complex scalar','Analysis',1.0d0,GiD_onNodes,GiD_NULL,GiD_NULL,'real_part','imaginary_part')
  do idx=1,9
    value=idx*1.5
    CALL GiD_WriteComplexScalar(idx,value,value*1.2)
  end do
  CALL GiD_EndResult
  
  !CALL GiD_ClosePostMeshFile
  CALL GiD_ClosePostResultFile()

END PROGRAM
