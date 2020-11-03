# Name of the program
EXE = nl-partsol
ROUTE_EXE = ../nl-partsol/build/./$(EXE)
TESTDIR = test

# Rule for benchmarks
.PHONY: Element Consolidation Debug_MMS Debug_MPM Debug_LME Sulsky_MPM Sulsky_GIMP Sulsky_LME Sulsky_SV Sulsky_NH Sulsky_NPC_FS Cuadrado_MPM Disco_MPM Disco_GIMP Taylor_bar Traccion_MPM Traccion_GIMP Traccion_LME Dyka_MPM Dyka_GIMP Dyka_LME Dyka_GA Brazilian_MPM Three_Point_Bar

# Caso para debugear I
Element:
	@echo "**************** Run test case ******************"
	$(ROUTE_EXE) -V Unix/Element/Element.gdf 
	@echo " "
	@echo "*************** Test case passed ****************"

# Caso para debugear II
Debug_MPM:
	@echo "**************** Run test case ******************"
	$(ROUTE_EXE) -2D -V Unix/Debug_MPM/ejemplo.gdf 
	@echo " "
	@echo "*************** Test case passed ****************"
Debug_LME:
	@echo "**************** Run test case ******************"
	$(ROUTE_EXE) -2D -V Unix/Debug_LME/ejemplo.gdf
	@echo " "
	@echo "*************** Test case passed ****************"

Debug_MMS:
	@echo "**************** Run test case ******************"
	$(ROUTE_EXE) -V Unix/Debug_MMS/ejemplo.gdf
	@echo " "
	@echo "*************** Test case passed ****************"

# Caso impacto elastico discos Sulsky
Sulsky_MPM:
	@echo "**************** Run test case ******************"
	$(ROUTE_EXE) -V Unix/Sulsky_MPM/ejemplo.gdf
	@echo " "
	@echo "*************** Test case passed ****************"
Sulsky_GIMP:
	@echo "**************** Run test case ******************"
	$(ROUTE_EXE) -V Unix/Sulsky_GIMP/ejemplo.gdf
	@echo " "
	@echo "*************** Test case passed ****************"
Sulsky_LME:
	@echo "**************** Run test case ******************"
	$(ROUTE_EXE) -V Unix/Sulsky_LME/ejemplo.gdf
	@echo " "
	@echo "*************** Test case passed ****************"
Sulsky_SV:
	@echo "**************** Run test case ******************"
	$(ROUTE_EXE) -V Unix/Sulsky_Saint_Venant/ejemplo.gdf
	@echo " "
	@echo "*************** Test case passed ****************"
Sulsky_NH:
	@echo "**************** Run test case ******************"
	$(ROUTE_EXE) -V Unix/Sulsky_Neo_Hookeano/ejemplo.gdf
	@echo " "
	@echo "*************** Test case passed ****************"
Sulsky_NPC_FS:
	@echo "**************** Run test case ******************"
	$(ROUTE_EXE) -V Unix/Sulsky_NPC_FS/ejemplo.gdf
	@echo " "
	@echo "*************** Test case passed ****************"

# Consolidacion Neo-Hookeno
Consolidation:
	@echo "**************** Run test case ******************"
	$(ROUTE_EXE) -V Unix/Consolidation/ejemplo.gdf
	@echo " "
	@echo "*************** Test case passed ****************"

# Caso cuadrado caida
Cuadrado_MPM:
	@echo "**************** Run test case ******************"
	$(ROUTE_EXE) -2D -V Unix/Cuadrado_MPM/ejemplo.gdf
	@echo " "
	@echo "*************** Test case passed ****************"

# Caso disco caida
Disco_MPM:
	@echo "**************** Run test case ******************"
	$(ROUTE_EXE) -2D -V Unix/Disco_MPM/ejemplo.gdf
	@echo " "
	@echo "*************** Test case passed ****************"
Disco_GIMP:
	@echo "**************** Run test case ******************"
	$(ROUTE_EXE) -2D -V Unix/Disco_GIMP/ejemplo.gdf
	@echo " "
	@echo "*************** Test case passed ****************"

# Barra de Taylor
Taylor_bar:
	@echo "**************** Run test case ******************"
	$(ROUTE_EXE) -V Unix/Taylor_bar/ejemplo.gdf
	@echo " "
	@echo "*************** Test case passed ****************"

# Caso barra traccion
Traccion_MPM:
	@echo "**************** Run test case ******************"
	$(ROUTE_EXE) -2D -V Unix/Traccion_MPM/ejemplo.gdf
	@echo " "
	@echo "*************** Test case passed ****************"
Traccion_GIMP:
	@echo "**************** Run test case ******************"
	$(ROUTE_EXE) -2D -V Unix/Traccion_GIMP/ejemplo.gdf
	@echo " "
	@echo "*************** Test case passed ****************"
Traccion_LME:
	@echo "**************** Run test case ******************"
	$(ROUTE_EXE) -2D -V Unix/Traccion_LME/ejemplo.gdf
	@echo " "
	@echo "*************** Test case passed ****************"

# Caso barra Dyka Thomas Blanc thesis
Dyka_MPM:
	@echo "**************** Run test case ******************"
	$(ROUTE_EXE) -2D -V Unix/Dyka_MPM/ejemplo.gdf
	@echo " "
	@echo "*************** Test case passed ****************"
Dyka_GIMP:
	@echo "**************** Run test case ******************"
	$(ROUTE_EXE) -2D -V Unix/Dyka_GIMP/ejemplo.gdf
	@echo " "
	@echo "*************** Test case passed ****************"
Dyka_LME:
	@echo "**************** Run test case ******************"
	$(ROUTE_EXE) -2D -V Unix/Dyka_LME/ejemplo.gdf
	@echo " "
	@echo "*************** Test case passed ****************"
Dyka_GenAlpha:
	@echo "**************** Run test case ******************"
	$(ROUTE_EXE) -2D -V Unix/Dyka_GenAlpha/ejemplo.gdf
	@echo " "
	@echo "*************** Test case passed ****************"

# Brazilian test for fracture
Brazilian_MPM:
	@echo "**************** Run test case ******************"
	$(ROUTE_EXE) -2D -V Unix/Brazilian_MPM/ejemplo.gdf
	@echo " "
	@echo "*************** Test case passed ****************"
Brazilian_LME:
	@echo "**************** Run test case ******************"
	$(ROUTE_EXE) -2D -V Unix/Brazilian_LME/ejemplo.gdf
	@echo " "
	@echo "*************** Test case passed ****************"
# Three point test fracture 
Three_Point_Bar:
	@echo "**************** Run test case ******************"
	$(ROUTE_EXE) -2D -V Unix/Three_Point_Bar/ejemplo.gdf
	@echo " "
	@echo "*************** Test case passed ****************"

#####################################################################
