# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 2.8

#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
.SUFFIXES:

# Remove some rules from gmake that .SUFFIXES does not remove.
SUFFIXES =

.SUFFIXES: .hpux_make_needs_suffix_list

# Suppress display of executed commands.
$(VERBOSE).SILENT:

# A target that is always out of date.
cmake_force:
.PHONY : cmake_force

#=============================================================================
# Set environment variables for the build.

# The shell in which to execute make rules.
SHELL = /bin/sh

# The CMake executable.
CMAKE_COMMAND = /usr/bin/cmake

# The command to remove a file.
RM = /usr/bin/cmake -E remove -f

# Escaping for special characters.
EQUALS = =

# The program to use to edit the cache.
CMAKE_EDIT_COMMAND = /usr/bin/cmake-gui

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /home/miguel/GID-git/gidproject/gidpost

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/miguel/GID-git/gidproject/gidpost/buildL

# Include any dependencies generated for this target.
include examples/CMakeFiles/BugMeshGroup03.dir/depend.make

# Include the progress variables for this target.
include examples/CMakeFiles/BugMeshGroup03.dir/progress.make

# Include the compile flags for this target's objects.
include examples/CMakeFiles/BugMeshGroup03.dir/flags.make

examples/CMakeFiles/BugMeshGroup03.dir/BugMeshGroup03.c.o: examples/CMakeFiles/BugMeshGroup03.dir/flags.make
examples/CMakeFiles/BugMeshGroup03.dir/BugMeshGroup03.c.o: ../examples/BugMeshGroup03.c
	$(CMAKE_COMMAND) -E cmake_progress_report /home/miguel/GID-git/gidproject/gidpost/buildL/CMakeFiles $(CMAKE_PROGRESS_1)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Building C object examples/CMakeFiles/BugMeshGroup03.dir/BugMeshGroup03.c.o"
	cd /home/miguel/GID-git/gidproject/gidpost/buildL/examples && /usr/bin/cc  $(C_DEFINES) $(C_FLAGS) -o CMakeFiles/BugMeshGroup03.dir/BugMeshGroup03.c.o   -c /home/miguel/GID-git/gidproject/gidpost/examples/BugMeshGroup03.c

examples/CMakeFiles/BugMeshGroup03.dir/BugMeshGroup03.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/BugMeshGroup03.dir/BugMeshGroup03.c.i"
	cd /home/miguel/GID-git/gidproject/gidpost/buildL/examples && /usr/bin/cc  $(C_DEFINES) $(C_FLAGS) -E /home/miguel/GID-git/gidproject/gidpost/examples/BugMeshGroup03.c > CMakeFiles/BugMeshGroup03.dir/BugMeshGroup03.c.i

examples/CMakeFiles/BugMeshGroup03.dir/BugMeshGroup03.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/BugMeshGroup03.dir/BugMeshGroup03.c.s"
	cd /home/miguel/GID-git/gidproject/gidpost/buildL/examples && /usr/bin/cc  $(C_DEFINES) $(C_FLAGS) -S /home/miguel/GID-git/gidproject/gidpost/examples/BugMeshGroup03.c -o CMakeFiles/BugMeshGroup03.dir/BugMeshGroup03.c.s

examples/CMakeFiles/BugMeshGroup03.dir/BugMeshGroup03.c.o.requires:
.PHONY : examples/CMakeFiles/BugMeshGroup03.dir/BugMeshGroup03.c.o.requires

examples/CMakeFiles/BugMeshGroup03.dir/BugMeshGroup03.c.o.provides: examples/CMakeFiles/BugMeshGroup03.dir/BugMeshGroup03.c.o.requires
	$(MAKE) -f examples/CMakeFiles/BugMeshGroup03.dir/build.make examples/CMakeFiles/BugMeshGroup03.dir/BugMeshGroup03.c.o.provides.build
.PHONY : examples/CMakeFiles/BugMeshGroup03.dir/BugMeshGroup03.c.o.provides

examples/CMakeFiles/BugMeshGroup03.dir/BugMeshGroup03.c.o.provides.build: examples/CMakeFiles/BugMeshGroup03.dir/BugMeshGroup03.c.o

# Object files for target BugMeshGroup03
BugMeshGroup03_OBJECTS = \
"CMakeFiles/BugMeshGroup03.dir/BugMeshGroup03.c.o"

# External object files for target BugMeshGroup03
BugMeshGroup03_EXTERNAL_OBJECTS =

examples/BugMeshGroup03: examples/CMakeFiles/BugMeshGroup03.dir/BugMeshGroup03.c.o
examples/BugMeshGroup03: examples/CMakeFiles/BugMeshGroup03.dir/build.make
examples/BugMeshGroup03: source/libgidpost.a
examples/BugMeshGroup03: /home/miguel/lib/libhdf5.so
examples/BugMeshGroup03: /home/miguel/lib/libz.so
examples/BugMeshGroup03: /usr/lib64/libdl.so
examples/BugMeshGroup03: /usr/lib64/libm.so
examples/BugMeshGroup03: /usr/lib64/libz.so
examples/BugMeshGroup03: examples/CMakeFiles/BugMeshGroup03.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --red --bold "Linking C executable BugMeshGroup03"
	cd /home/miguel/GID-git/gidproject/gidpost/buildL/examples && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/BugMeshGroup03.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
examples/CMakeFiles/BugMeshGroup03.dir/build: examples/BugMeshGroup03
.PHONY : examples/CMakeFiles/BugMeshGroup03.dir/build

examples/CMakeFiles/BugMeshGroup03.dir/requires: examples/CMakeFiles/BugMeshGroup03.dir/BugMeshGroup03.c.o.requires
.PHONY : examples/CMakeFiles/BugMeshGroup03.dir/requires

examples/CMakeFiles/BugMeshGroup03.dir/clean:
	cd /home/miguel/GID-git/gidproject/gidpost/buildL/examples && $(CMAKE_COMMAND) -P CMakeFiles/BugMeshGroup03.dir/cmake_clean.cmake
.PHONY : examples/CMakeFiles/BugMeshGroup03.dir/clean

examples/CMakeFiles/BugMeshGroup03.dir/depend:
	cd /home/miguel/GID-git/gidproject/gidpost/buildL && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/miguel/GID-git/gidproject/gidpost /home/miguel/GID-git/gidproject/gidpost/examples /home/miguel/GID-git/gidproject/gidpost/buildL /home/miguel/GID-git/gidproject/gidpost/buildL/examples /home/miguel/GID-git/gidproject/gidpost/buildL/examples/CMakeFiles/BugMeshGroup03.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : examples/CMakeFiles/BugMeshGroup03.dir/depend
