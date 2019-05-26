# Makefile for the FEM_Miguel

# Name of the program
myprog = GeHoJerez

# Compiler and preprocesor
CC = gcc

# Source code
csrc = \
	$(wildcard src/*.c)\
	$(wildcard src/ToolsLib/*.c)\
	$(wildcard src/InOutFun/*.c)\
	$(wildcard src/Constitutive/*.c)\
	$(wildcard src/Solvers/*.c)\
	$(wildcard src/ElementsFunctions/*.c)\
	$(wildcard src/GaussPointsFunctions/*.c)\
	$(wildcard src/Num_Schemes/*.c)\
	$(wildcard src/MPM_Subroutines/*.c)\
	$(wildcard src/Boundary_Conditions/*.c)

# Transforms the contents of the src variable, changing all file suffixes
# from .c to .o, thus constructing the object file list we need.
obj = $(csrc:.c=.o)

# List of libraries required during linking
LDFLAGS = -lm

# compile
$(myprog): $(obj)
	$(CC) -o $@  $^ $(LDFLAGS)
	@echo " "
	@echo "*************  GeHoJerez compiled !! ************"
	@echo " "
	@echo "**************** Run test case ******************"
	./$(myprog) test/MPM_test/CasoSimple/ejemplo.dat
	@echo " "
	@echo "*************** Test case passed ****************"
# run
test1D:
	./$(myprog) test/FEM_test/ejemplo.dat
test2D:
	./$(myprog) test/FEM_test/ejemplo.dat

# New rule for cleaning up every target, in order to rebuild
# the whole program from scratch
.PHONY: clean
clean:
	rm -f $(obj) $(myprog)

