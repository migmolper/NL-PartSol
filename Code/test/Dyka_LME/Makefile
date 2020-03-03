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
# Rule for benchmarks
.PHONY: Dyka


# Caso barra Dyka Thomas Blanc thesis
Dyka:
	@echo "**************** Run test case ******************"
	$ ./$(EXE) -2D -V Dyka.gdf
	@echo " "
	@echo "*************** Test case passed ****************"
