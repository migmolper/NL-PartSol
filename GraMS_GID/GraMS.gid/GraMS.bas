//-----------------------------------------------------------------------
FILE          
*GenData(1)
//-----------------------------------------------------------------------
DIMENSIONS
NNODE= *npoin
NCONE= *nnode
NDIME= *ndime
NELEM= *nelem
NMATS= *nmats
//-----------------------------------------------------------------------
GEOMETRY
CONNECTIVITY
*loop elems
  *elemsnum *elemsmat *elemsConec(swap)
*end elems
END_CON
NODES
*loop nodes
*format "%6i%15.5f%15.5f%15.5f"
  *NodesNum *NodesCoord
*end
END_NODES
END_GEOMETRY
//-----------------------------------------------------------------------
BOUNDARY_CONDITIONS
*Set Cond BC_NSET *nodes
*Add Cond BC_LSET *nodes
*Set var valor1 = 1
*loop nodes *OnlyInCond
*Set var valor = Cond(1,Int)
*if(valor > valor1)
*Set var valor1 = valor
*endif
*end nodes
*for(i=1;i<=valor1;i=i+1)
*Set var total1 = 0
*loop nodes *OnlyInCond
*Set var valor = Cond(1,Int)
*if(valor == i)
*Set var total1 = total1 + 1
*endif
*end
BC_SET *i *total1
*loop nodes *OnlyInCond
*Set var valor = Cond(1,Int)
*if(valor == i)
*NodesNum
*endif
*end
*end
END_BC
//-----------------------------------------------------------------------
LOADS
*Set Cond POINT_LOAD *nodes
*Add Cond POINT_LOAD_LINE *nodes
*Set var valor1 = 1
*loop nodes *OnlyInCond
*Set var valor = Cond(1,Int)
*if(valor > valor1)
*Set var valor1 = valor
*endif
*end nodes
*for(i=1;i<=valor1;i=i+1)
*Set var total1 = 0
*loop nodes *OnlyInCond
*Set var valor = Cond(1,Int)
*if(valor == i)
*Set var total1 = total1 + 1
*endif
*end
POINT_LOAD *i *total1
*loop nodes *OnlyInCond
*Set var valor = Cond(1,Int)
*if(valor == i)
*NodesNum
*endif
*end
*end
*Set Cond SURFACE_LOAD *elems
*Set var valor1 = 1
*loop elems *OnlyInCond
*Set var valor = Cond(1,Int)
*if(valor > valor1)
*Set var valor1 = valor
*endif
*end elemes
*for(i=1;i<=valor1;i=i+1)
*Set var total1 = 0
*loop elems *OnlyInCond
*Set var valor = Cond(1,Int)
*if(valor == i)
*Set var total1 = total1 + 1
*endif
*end
LINE_LOAD *i *total1
*loop elems *OnlyInCond
*Set var valor = Cond(1,Int)
*if(valor == i)
*globalnodes
*endif
*end
*end
*Set Cond VOLUME_LOAD *nodes
*Set var valor1 = 1
*loop nodes *OnlyInCond
*Set var valor = Cond(1,Int)
*if(valor > valor1)
*Set var valor1 = valor
*endif
*end nodes
*for(i=1;i<=valor1;i=i+1)
*Set var total1 = 0
*loop nodes *OnlyInCond
*Set var valor = Cond(1,Int)
*if(valor == i)
*Set var total1 = total1 + 1
*endif
*end
VOLUME_LOAD *i *total1
*loop nodes *OnlyInCond
*Set var valor = Cond(1,Int)
*if(valor == i)
*NodesNum
*endif
*end
*end
END_LDS
