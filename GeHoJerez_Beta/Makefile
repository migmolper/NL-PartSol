# Makefile for the GeHoJerez

# Compiler root
GCC_ROOT = /usr/bin

# Name of the directories
SRCDIR = src
OBJDIR = obj
DEBDIR = bin/Debug
RELDIR = bin/Release
TESTDIR = test
INSTDIR = /usr/bin

# Name of the program
EXE = GeHoJerez_V1_Linux

# Compiler and preprocesor
CC = $(GCC_ROOT)/gcc

#####################################################################

# Source code
SOURCES = \
	$(wildcard src/*.c)\
	$(wildcard src/ToolsLib/*.c)\
	$(wildcard src/InOutFun/*.c)\
	$(wildcard src/Constitutive/*.c)\
	$(wildcard src/Solvers/*.c)\
	$(wildcard src/ElementsFunctions/*.c)\
	$(wildcard src/GaussPointsFunctions/*.c)\
	$(wildcard src/MPM_Subroutines/*.c)\
	$(wildcard src/Boundary_Conditions/*.c)

#####################################################################

# Transforms the contents of the src variable, changing all file suffixes
# from .c to .o, thus constructing the object file list we need.
OBJS = $(SOURCES:.c=.o)

#####################################################################

# List of libraries required during linking
LDFLAGS = -lm
DBFLAGS = -g -Wall -Wpedantic -Wextra 
RLFLAGS = -O2

#####################################################################

# Debug recipie
debug: $(DEBDIR)/$(EXE)

# Link process
$(DEBDIR)/$(EXE): $(OBJS)
	@echo "*************************************************"
	$(CC) -o $@ $^ $(LDFLAGS)
	@echo " "
	@echo "*************  GeHoJerez compiled !! ************"
	@echo " "
# Compile	
%.o: %.c
	@echo "****************** Debug mode *******************"
	$(CC) -o $@ -c $(DBFLAGS) $<
	@echo "\t Compiled "$<" successfully !!"

#####################################################################

# Release recipie
release: $(RELDIR)/$(EXE)

# Link process
$(RELDIR)/$(EXE): $(OBJS)
	@echo "*************************************************"
	$(CC) -o $@ $^ $(LDFLAGS)
	@echo " "
	@echo "*************  GeHoJerez compiled !! ************"
	@echo " "
# Compile	
%.o: %.c
	@echo "***************** Release mode ******************"
	$(CC) -o $@ -c $(DBFLAGS) $<
	@echo "\t Compiled "$<" successfully !!"

#####################################################################
# Rule for benchmarks
.PHONY: Benchmarks Disco2D Cuadrado2D Barra1D

# Run al the Benchmarks
Benchmarks: Disco2D Cuadrado2D Barra1D

# Caso Disco
Disco2D:
	@echo "**************** Run test case ******************"
	$(DEBDIR)/./$(EXE) $(TESTDIR)/CasoDisco/ejemplo.dat
	@echo " "
	@echo "*************** Test case passed ****************"

# Caso cuadrado
Cuadrado2D:
	@echo "**************** Run test case ******************"
	$(DEBDIR)/./$(EXE) $(TESTDIR)/CasoCuadrado/ejemplo.dat
	@echo " "
	@echo "*************** Test case passed ****************"

# Caso barra 1D
Barra1D:
	@echo "**************** Run test case ******************"
	$(DEBDIR)/./$(EXE) $(TESTDIR)/CasoBarraVertical/ejemplo.dat
	@echo " "
	@echo "*************** Test case passed ****************"

#####################################################################
# Rule for installing GeHoJerez
.PHONY: install
install: $(RELDIR)/$(EXE)
	@echo "************** Install GeHoJerez ****************"
	mkdir -p $(INSTDIR)/$(EXE)
	cp $< $(INSTDIR)/$(EXE)
	@echo "************ Installed successfully *************"

.PHONY: uninstall
uninstall:
	@echo "************* Uninstall GeHoJerez ***************"
	rm -rf $(INSTDIR)/$(EXE)
	@echo "*********** Uninstalled successfully ************"

#####################################################################

# Rule for cleaning up every target, in order to rebuild
# the whole program from scratch
.PHONY: clean
clean:
	rm -f $(OBJS)
	rm -f $(BINDIR)/$(DEB)/$(EXE)
	rm -f $(BINDIR)/$(REL)/$(EXE)
	@echo "Cleanup complete!"
