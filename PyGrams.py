def Write_GramsBox(FileName, TypeMesh, FileMesh, Field, Values):

    """
    GramsBox (Type=GID,File=FEM_Mesh.msh) {
    left=GramsBoundary {
        BcDirichlet V.x test/Dyka_MPM/CurveConstat.txt
 	BcDirichlet V.y NULL	
    }
    right=GramsBoundary {
  	BcDirichlet V.x NULL
   	BcDirichlet V.y NULL	
    }
    top=GramsBoundary {
      	BcDirichlet V.x NULL
    	BcDirichlet V.y NULL
    }
    bottom=GramsBoundary {
  	BcDirichlet V.x NULL
    	BcDirichlet V.y NULL
    }
    }
    
    """
    GDF_file = open('%s.gdf' %(FileName),'a')

    # Open
    GDF_file.write('GramsBox (Type=%s,File=%s) { \n' %(TypeMesh,FileMesh))

    if(Field == 'V'):
        # Left boundary
        GDF_file.write('\t %s=%s { \n' %('left','GramsBoundary'))
        GDF_file.write('\t \t %s %s %s \n' %('BcDirichlet','V.x',Values[0][0]))
        GDF_file.write('\t \t %s %s %s \n' %('BcDirichlet','V.y',Values[0][1]))
        GDF_file.write('\t } \n')
        
        # Right boundary
        GDF_file.write('\t %s=%s { \n' %('right','GramsBoundary'))
        GDF_file.write('\t \t %s %s %s \n' %('BcDirichlet','V.x',Values[1][0]))
        GDF_file.write('\t \t %s %s %s \n' %('BcDirichlet','V.y',Values[1][1]))
        GDF_file.write('\t } \n')
        
        # Top boundary
        GDF_file.write('\t %s=%s { \n' %('top','GramsBoundary'))
        GDF_file.write('\t \t %s %s %s \n' %('BcDirichlet','V.x',Values[2][0]))
        GDF_file.write('\t \t %s %s %s \n' %('BcDirichlet','V.y',Values[2][1]))
        GDF_file.write('\t } \n')
        
        # Bottom boundary
        GDF_file.write('\t %s=%s { \n' %('bottom','GramsBoundary'))
        GDF_file.write('\t \t %s %s %s \n' %('BcDirichlet','V.x',Values[3][0]))
        GDF_file.write('\t \t %s %s %s \n' %('BcDirichlet','V.y',Values[3][1]))
        GDF_file.write('\t } \n')


    # Close
    GDF_file.write('} \n')
    GDF_file.close()

def Write_GramsTime(FileName, TypeAnalysis, Properties):
    """
    GramsTime (Type=FE) {
	CEL=5053.02
	CFL=0.2
	N=800000
    }
    """
    GDF_file = open('%s.gdf' %(FileName),'a')

    # Open
    GDF_file.write('\t GramsTime (Type=%s) { \n' %(TypeAnalysis))

    # Properties
    GDF_file.write('\t \t CEL=%f \n' %(Properties[0]))
    GDF_file.write('\t \t CFL=%f \n' %(Properties[1]))
    GDF_file.write('\t \t N=%i \n' %(Properties[2]))

    # Close
    GDF_file.write('\t } \n')
    GDF_file.close()

def Write_GramsMaterials(FileName, IdMaterial, Properties):
    """
    GramsMaterials (id=0) {
	Type=LE
	rho=7833
	E=200.e9
	mu=0.0
    }
    """
    GDF_file = open('%s.gdf' %(FileName),'a')

    # Open
    GDF_file.write('\t GramsMaterials (id=%i) { \n' %(IdMaterial))

    # Properties
    GDF_file.write('\t \t Type=%s \n' %(Properties[0]))
    GDF_file.write('\t \t rho=%f \n' %(Properties[1]))
    GDF_file.write('\t \t E=%f \n' %(Properties[2]))
    GDF_file.write('\t \t mu=%f \n' %(Properties[3]))

    # Close
    GDF_file.write('\t } \n')
    GDF_file.close()

    
def Write_GramsShapeFun(FileName, TypeShapeFun):
    """
    GramsShapeFun (Type=MPMQ4) {
    }
    """
    GDF_file = open('%s.gdf' %(FileName),'a')

    # Open
    GDF_file.write('\t GramsShapeFun (type=%s) { \n' %(TypeShapeFun))

    # Close
    GDF_file.write('\t } \n')
    GDF_file.close()

def Write_GramsInitials(FileName,FileNodes,Value):
    """
    GramsInitials (Nodes=ListInit.txt) {
	Value=[5.0,0.0]
    }
    """
    GDF_file = open('%s.gdf' %(FileName),'a')

    # Open
    GDF_file.write('\t GramsInitials (Nodes=%s) { \n' %(FileNodes))

    # Properties
    GDF_file.write('\t \t Value=[%f,%f] \n' %(Value[0],Value[1]))    

    # Close
    GDF_file.write('\t } \n')
    GDF_file.close()

def Write_GramsSolid2D(FileName,FileMesh):

    # Open
    GDF_file = open('%s.gdf' %(FileName),'a')
    GDF_file.write('GramsSolid2D (File=%s) { \n' %(FileMesh))
    GDF_file.close()

    # Properties
    Write_GramsTime(FileName, 'FE', [5053.02,0.2,800000])
    Write_GramsMaterials(FileName, 0, ['LE',7833,200.e9,0.0])
    Write_GramsShapeFun(FileName, 'MPMQ4')
    Write_GramsInitials(FileName,'ListInit.txt',[5.0,0.0])

    # Close
    GDF_file = open('%s.gdf' %(FileName),'a')
    GDF_file.write('} \n')
    GDF_file.close()

    
    
# # Material points definitions
# GramsSolid2D (File=MPM_Mesh.msh) {

#     # Loads over the GPs
#     # GramsBodyForces (Nodes=) {}
#     # GramsContactForces (Nodes=) {}
  
# }

# # Outputs 
# GramsOutputs (i=10000) {
# 	     DIR=test/Dyka_MPM	
# }



