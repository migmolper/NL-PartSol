# Makefile for the GeHoJerez

# Compiler root
GCC_ROOT = /usr/bin

# Name of the directories
SRCDIR = src
OBJDIR = obj
# DEBDIR = bin/Debug
# RELDIR = bin/Release
TESTDIR = test
INSTDIR = /usr/bin

# Name of the program
EXE = GraMS

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
OBJS = $(SOURCES:.c=.c.o)

#####################################################################

# List of libraries required during linking
LDFLAGS = -lm
DBFLAGS = -g -Wall -Wpedantic -Wextra 
# RLFLAGS = -O2

#####################################################################

# Debug recipie
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

# # Release recipie
# # release: $(RELDIR)/$(EXE)
# release: $(EXE)

# # Link process
# # $(RELDIR)/$(EXE): $(OBJS)
# $(EXE): $(OBJS)
# 	@echo "*************************************************"
# 	$(CC) -o $@ $^ $(LDFLAGS)
# 	@echo " "
# 	@echo "*************  GeHoJerez compiled !! ************"
# 	@echo " "
# # Compile	
# %.o: %.c
# 	@echo "***************** Release mode ******************"
# 	$(CC) -o $@ -c $(DBFLAGS) $<
# 	@echo "\t Compiled "$<" successfully !!"

#####################################################################
# Rule for benchmarks
.PHONY: Benchmarks Disco2D Cuadrado2D Barra1D_caida Barra1D_traccion Dyka1D Dyka1D_Blanc

# Run al the Benchmarks
Benchmarks: Disco2D Cuadrado2D Barra1D_caida Barra1D_traccion Dyka1D_Blanc

# Caso Disco
Disco2D:
	@echo "**************** Run test case ******************"
	$ ./$(EXE) $(TESTDIR)/CasoDisco/ejemplo.dat
	# $(DEBDIR)/./$(EXE) $(TESTDIR)/CasoDisco/ejemplo.dat
	@echo " "
	@echo "*************** Test case passed ****************"

# Caso cuadrado
Cuadrado2D:
	@echo "**************** Run test case ******************"
	$ ./$(EXE) $(TESTDIR)/CasoCuadrado/ejemplo.dat
	# $(DEBDIR)/./$(EXE) $(TESTDIR)/CasoCuadrado/ejemplo.dat
	@echo " "
	@echo "*************** Test case passed ****************"

# Caso barra caida 1D
Barra1D_caida:
	@echo "**************** Run test case ******************"
	# $(DEBDIR)/./$(EXE) $(TESTDIR)/CasoBarraVertical/ejemplo.dat
	$ ./$(EXE) $(TESTDIR)/CasoBarraVertical/ejemplo.dat
	@echo " "
	@echo "*************** Test case passed ****************"

# Caso barra traccion 1D
Barra1D_traccion:
	@echo "**************** Run test case ******************"
	# $(DEBDIR)/./$(EXE) $(TESTDIR)/CasoBarraVertical/ejemplo.dat
	$ ./$(EXE) $(TESTDIR)/BarraTraccion1D/ejemplo.dat
	@echo " "
	@echo "*************** Test case passed ****************"
# Caso barra Dyka
Dyka1D:
	@echo "**************** Run test case ******************"
	$ ./$(EXE) $(TESTDIR)/CasoDyka/ejemplo.dat
	@echo " "
	@echo "*************** Test case passed ****************"

# Caso barra Dyka Thomas Blanc thesis
Dyka1D_Blanc:
	@echo "**************** Run test case ******************"
	$ ./$(EXE) $(TESTDIR)/CasoDyka_Pastor_Blanc/ejemplo.dat
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
