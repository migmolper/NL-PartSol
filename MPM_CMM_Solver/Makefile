# Makefile for the FEM_Miguel

# Name of the program
myprog = FEM_Miguel

# Compiler and preprocesor
CC = gcc

# Source code
csrc = $(wildcard src/utilities/*.c) \
	$(wildcard src/mesh/*.c) \
	$(wildcard src/inputs/*.c)\
	$(wildcard src/*.c) 

# Collects all source files in the current directory in the src variable
src = $(wildcard *.c)

# Transforms the contents of the src variable, changing all file suffixes
# from .c to .o, thus constructing the object file list we need.
obj = $(csrc:.c=.o)

# List of libraries required during linking
LDFLAGS = -lm


# 
$(myprog): $(obj)
	$(CC) -o $@ $^ $(LDFLAGS)

# New rule for cleaning up every target, in order to rebuild
# the whole program from scratch
.PHONY: clean
clean:
	rm -f $(obj) $(myprog)
