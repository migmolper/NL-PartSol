# Makefile for the FEM_Miguel

# Name of the program
myprog = GIMP

# Compiler and preprocesor
CC = gcc

# Source code
csrc = \
	$(wildcard src/*.c)\
	$(wildcard src/ToolsLib/*.c)\
	$(wildcard src/ReadData/*.c)\
	$(wildcard src/Constitutive/*.c)\
	$(wildcard src/Solvers/*.c)\
	$(wildcard src/ElementsFunctions/*.c)\
	$(wildcard src/GaussPointsFunctions/*.c)\
	$(wildcard src/Forces/*.c)

# Transforms the contents of the src variable, changing all file suffixes
# from .c to .o, thus constructing the object file list we need.
obj = $(csrc:.c=.o)

# List of libraries required during linking
LDFLAGS = -lm

# compile
$(myprog): $(obj)
	$(CC) -o $@  $^ $(LDFLAGS)
	./$(myprog)

# run
run:
	./$(myprog)
# New rule for cleaning up every target, in order to rebuild
# the whole program from scratch
.PHONY: clean
clean:
	rm -f $(obj) $(myprog)

