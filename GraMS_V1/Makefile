# Makefile for the GraMS

# Compiler root
GCC_ROOT = /usr/bin

# Name of the directories
SRCDIR = src
OBJDIR = obj
TESTDIR = test
INSTDIR = /usr/bin

# Name of the program
EXE = GraMS

# Compiler and preprocesor
CC = $(GCC_ROOT)/gcc

#####################################################################

# Source code
SOURCES = \
	$(wildcard src/GRAMS/*.c)\
	$(wildcard src/ToolsLib/*.c)\
	$(wildcard src/InOutFun/*.c)\
	$(wildcard src/Constitutive/*.c)\
	$(wildcard src/Matlib/*.c)\
	$(wildcard src/ShapeFun/*.c)\
	$(wildcard src/MPM/*.c)\
	$(wildcard src/Boundary_Conditions/*.c)\
	$(wildcard src/Formulations/*.c)

#####################################################################

# Transforms the contents of the src variable, changing all file suffixes
# from .c to .o, thus constructing the object file list we need.
OBJS = $(SOURCES:.c=.c.o)

#####################################################################

# List of libraries required during linking
LDFLAGS = -lm
DBFLAGS = -g -Wall -Wpedantic -Wextra 
# RLFLAGS = -O2

#####################################################################

# Compile recipie
# debug: $(DEBDIR)/$(EXE)
debug: $(EXE)

# Link process
# $(DEBDIR)/$(EXE): $(OBJS)
$(EXE): $(OBJS)
	@echo "*************************************************"
	$(CC) -o $@ $^ $(LDFLAGS)
	@echo " "
	@echo "***** GRAnular Mechanics Solver compiled !! *****"
	@echo " "
# Compile	
%.c.o: %.c
	@echo "****************** Debug mode *******************"
	$(CC) -o $@ -c $(DBFLAGS) $<
	@echo "\t Compiled "$<" successfully !!"

#####################################################################
# Rule for benchmarks
.PHONY: Debug_MPM Debug_LME Sulsky_MPM Sulsky_GIMP Sulsky_LME Cuadrado_MPM Disco_MPM Disco_GIMP Traccion_MPM Traccion_GIMP Dyka_MPM Dyka_GIMP Brazilian_MPM

# Caso para debugear
Debug_MPM:
	@echo "**************** Run test case ******************"
	$ ./$(EXE) $(TESTDIR)/Debug_MPM/ejemplo.gdf
	@echo " "
	@echo "*************** Test case passed ****************"
Debug_LME:
	@echo "**************** Run test case ******************"
	$ ./$(EXE) $(TESTDIR)/Debug_LME/ejemplo.gdf
	@echo " "
	@echo "*************** Test case passed ****************"

# Caso impacto elastico discos Sulsky
Sulsky_MPM:
	@echo "**************** Run test case ******************"
	$ ./$(EXE) $(TESTDIR)/Sulsky_MPM/ejemplo.gdf
	@echo " "
	@echo "*************** Test case passed ****************"
Sulsky_GIMP:
	@echo "**************** Run test case ******************"
	$ ./$(EXE) $(TESTDIR)/Sulsky_GIMP/ejemplo.gdf
	@echo " "
	@echo "*************** Test case passed ****************"
Sulsky_LME:
	@echo "**************** Run test case ******************"
	$ ./$(EXE) $(TESTDIR)/Sulsky_LME/ejemplo.gdf
	@echo " "
	@echo "*************** Test case passed ****************"



# Caso cuadrado caida
Cuadrado_MPM:
	@echo "**************** Run test case ******************"
	$ ./$(EXE) $(TESTDIR)/Cuadrado_MPM/ejemplo.gdf
	@echo " "
	@echo "*************** Test case passed ****************"

# Caso disco caida
Disco_MPM:
	@echo "**************** Run test case ******************"
	$ ./$(EXE) $(TESTDIR)/Disco_MPM/ejemplo.gdf
	@echo " "
	@echo "*************** Test case passed ****************"
Disco_GIMP:
	@echo "**************** Run test case ******************"
	$ ./$(EXE) $(TESTDIR)/Disco_GIMP/ejemplo.gdf
	@echo " "
	@echo "*************** Test case passed ****************"

# Caso barra traccion
Traccion_MPM:
	@echo "**************** Run test case ******************"
	$ ./$(EXE) $(TESTDIR)/Traccion_MPM/ejemplo.gdf
	@echo " "
	@echo "*************** Test case passed ****************"
Traccion_GIMP:
	@echo "**************** Run test case ******************"
	$ ./$(EXE) $(TESTDIR)/Traccion_GIMP/ejemplo.gdf
	@echo " "
	@echo "*************** Test case passed ****************"

# Caso barra Dyka Thomas Blanc thesis
Dyka_MPM:
	@echo "**************** Run test case ******************"
	$ ./$(EXE) $(TESTDIR)/Dyka_MPM/ejemplo.gdf
	@echo " "
	@echo "*************** Test case passed ****************"
Dyka_GIMP:
	@echo "**************** Run test case ******************"
	$ ./$(EXE) $(TESTDIR)/Dyka_GIMP/ejemplo.gdf
	@echo " "
	@echo "*************** Test case passed ****************"

# Brazilian test for fracture
Brazilian_MPM:
	@echo "**************** Run test case ******************"
	$ ./$(EXE) $(TESTDIR)/Brazilian_MPM/ejemplo.gdf
	@echo " "
	@echo "*************** Test case passed ****************"	

#####################################################################
# Rule for installing GeHoJerez
# .PHONY: install
# install: $(RELDIR)/$(EXE)
# 	@echo "************** Install GeHoJerez ****************"
# 	mkdir -p $(INSTDIR)/$(EXE)
# 	cp $< $(INSTDIR)/$(EXE)
# 	@echo "************ Installed successfully *************"

# .PHONY: uninstall
# uninstall:
# 	@echo "************* Uninstall GeHoJerez ***************"
# 	rm -rf $(INSTDIR)/$(EXE)
# 	@echo "*********** Uninstalled successfully ************"

#####################################################################

# Rule for cleaning up every target, in order to rebuild
# the whole program from scratch
.PHONY: clean
clean:
	rm -f $(OBJS)
	rm -f $(EXE)
	# rm -f $(BINDIR)/$(DEB)/$(EXE)
	# rm -f $(BINDIR)/$(REL)/$(EXE)
	@echo "Cleanup complete!"
