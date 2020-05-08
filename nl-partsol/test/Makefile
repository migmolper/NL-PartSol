# Name of the program
EXE = nl-partsol
TESTDIR = test

# Rule for benchmarks
.PHONY: Debug_MPM Debug_LME Sulsky_MPM Sulsky_GIMP Sulsky_LME Cuadrado_MPM Disco_MPM Disco_GIMP Traccion_MPM Traccion_GIMP Traccion_LME Dyka_MPM Dyka_GIMP Dyka_LME Dyka_GA Brazilian_MPM Three_Point_Bar

# Caso para debugear
Debug_MPM:
	@echo "**************** Run test case ******************"
	$ .././$(EXE) -2D -V Debug_MPM/ejemplo.gdf 
	@echo " "
	@echo "*************** Test case passed ****************"
Debug_LME:
	@echo "**************** Run test case ******************"
	$ .././$(EXE) -2D -V Debug_LME/ejemplo.gdf
	@echo " "
	@echo "*************** Test case passed ****************"

# Caso impacto elastico discos Sulsky
Sulsky_MPM:
	@echo "**************** Run test case ******************"
	$ .././$(EXE) -2D -V Sulsky_MPM/ejemplo.gdf
	@echo " "
	@echo "*************** Test case passed ****************"
Sulsky_GIMP:
	@echo "**************** Run test case ******************"
	$ .././$(EXE) -2D -V Sulsky_GIMP/ejemplo.gdf
	@echo " "
	@echo "*************** Test case passed ****************"
Sulsky_LME:
	@echo "**************** Run test case ******************"
	$ .././$(EXE) -2D -V Sulsky_LME/ejemplo.gdf
	@echo " "
	@echo "*************** Test case passed ****************"



# Caso cuadrado caida
Cuadrado_MPM:
	@echo "**************** Run test case ******************"
	$ .././$(EXE) -2D -V Cuadrado_MPM/ejemplo.gdf
	@echo " "
	@echo "*************** Test case passed ****************"

# Caso disco caida
Disco_MPM:
	@echo "**************** Run test case ******************"
	$ .././$(EXE) -2D -V Disco_MPM/ejemplo.gdf
	@echo " "
	@echo "*************** Test case passed ****************"
Disco_GIMP:
	@echo "**************** Run test case ******************"
	$ .././$(EXE) -2D -V Disco_GIMP/ejemplo.gdf
	@echo " "
	@echo "*************** Test case passed ****************"

# Caso barra traccion
Traccion_MPM:
	@echo "**************** Run test case ******************"
	$ .././$(EXE) -2D -V Traccion_MPM/ejemplo.gdf
	@echo " "
	@echo "*************** Test case passed ****************"
Traccion_GIMP:
	@echo "**************** Run test case ******************"
	$ .././$(EXE) -2D -V Traccion_GIMP/ejemplo.gdf
	@echo " "
	@echo "*************** Test case passed ****************"
Traccion_LME:
	@echo "**************** Run test case ******************"
	$ .././$(EXE) -2D -V Traccion_LME/ejemplo.gdf
	@echo " "
	@echo "*************** Test case passed ****************"

# Caso barra Dyka Thomas Blanc thesis
Dyka_MPM:
	@echo "**************** Run test case ******************"
	$ .././$(EXE) -2D -V Dyka_MPM/ejemplo.gdf
	@echo " "
	@echo "*************** Test case passed ****************"
Dyka_GIMP:
	@echo "**************** Run test case ******************"
	$ .././$(EXE) -2D -V Dyka_GIMP/ejemplo.gdf
	@echo " "
	@echo "*************** Test case passed ****************"
Dyka_LME:
	@echo "**************** Run test case ******************"
	$ .././$(EXE) -2D -V Dyka_LME/ejemplo.gdf
	@echo " "
	@echo "*************** Test case passed ****************"
Dyka_GenAlpha:
	@echo "**************** Run test case ******************"
	$ .././$(EXE) -2D -V Dyka_GenAlpha/ejemplo.gdf
	$ python $(TESTDIR)/Dyka_GenAlpha/Velocity.py
	@echo " "
	@echo "*************** Test case passed ****************"

# Brazilian test for fracture
Brazilian_MPM:
	@echo "**************** Run test case ******************"
	$ .././$(EXE) -2D -V Brazilian_MPM/ejemplo.gdf
	@echo " "
	@echo "*************** Test case passed ****************"
Brazilian_LME:
	@echo "**************** Run test case ******************"
	$ .././$(EXE) -2D -V Brazilian_LME/ejemplo.gdf
	@echo " "
	@echo "*************** Test case passed ****************"
# Three point test fracture 
Three_Point_Bar:
	@echo "**************** Run test case ******************"
	$ .././$(EXE) -2D -V Three_Point_Bar/ejemplo.gdf
	@echo " "
	@echo "*************** Test case passed ****************"

#####################################################################
