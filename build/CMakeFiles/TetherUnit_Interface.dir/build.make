# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.18

# Delete rule output on recipe failure.
.DELETE_ON_ERROR:


#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
.SUFFIXES:


# Disable VCS-based implicit rules.
% : %,v


# Disable VCS-based implicit rules.
% : RCS/%


# Disable VCS-based implicit rules.
% : RCS/%,v


# Disable VCS-based implicit rules.
% : SCCS/s.%


# Disable VCS-based implicit rules.
% : s.%


.SUFFIXES: .hpux_make_needs_suffix_list


# Command-line flag to silence nested $(MAKE).
$(VERBOSE)MAKESILENT = -s

#Suppress display of executed commands.
$(VERBOSE).SILENT:

# A target that is always out of date.
cmake_force:

.PHONY : cmake_force

#=============================================================================
# Set environment variables for the build.

# The shell in which to execute make rules.
SHELL = /bin/sh

# The CMake executable.
CMAKE_COMMAND = /usr/local/bin/cmake

# The command to remove a file.
RM = /usr/local/bin/cmake -E rm -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /home/jerluen/atu-crt

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/jerluen/atu-crt/build

# Include any dependencies generated for this target.
include CMakeFiles/TetherUnit_Interface.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/TetherUnit_Interface.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/TetherUnit_Interface.dir/flags.make

CMakeFiles/TetherUnit_Interface.dir/src/tether_units/TetherUnit_Interface.cpp.o: CMakeFiles/TetherUnit_Interface.dir/flags.make
CMakeFiles/TetherUnit_Interface.dir/src/tether_units/TetherUnit_Interface.cpp.o: ../src/tether_units/TetherUnit_Interface.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/jerluen/atu-crt/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/TetherUnit_Interface.dir/src/tether_units/TetherUnit_Interface.cpp.o"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/TetherUnit_Interface.dir/src/tether_units/TetherUnit_Interface.cpp.o -c /home/jerluen/atu-crt/src/tether_units/TetherUnit_Interface.cpp

CMakeFiles/TetherUnit_Interface.dir/src/tether_units/TetherUnit_Interface.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/TetherUnit_Interface.dir/src/tether_units/TetherUnit_Interface.cpp.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/jerluen/atu-crt/src/tether_units/TetherUnit_Interface.cpp > CMakeFiles/TetherUnit_Interface.dir/src/tether_units/TetherUnit_Interface.cpp.i

CMakeFiles/TetherUnit_Interface.dir/src/tether_units/TetherUnit_Interface.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/TetherUnit_Interface.dir/src/tether_units/TetherUnit_Interface.cpp.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/jerluen/atu-crt/src/tether_units/TetherUnit_Interface.cpp -o CMakeFiles/TetherUnit_Interface.dir/src/tether_units/TetherUnit_Interface.cpp.s

# Object files for target TetherUnit_Interface
TetherUnit_Interface_OBJECTS = \
"CMakeFiles/TetherUnit_Interface.dir/src/tether_units/TetherUnit_Interface.cpp.o"

# External object files for target TetherUnit_Interface
TetherUnit_Interface_EXTERNAL_OBJECTS =

../lib/static/libTetherUnit_Interface.a: CMakeFiles/TetherUnit_Interface.dir/src/tether_units/TetherUnit_Interface.cpp.o
../lib/static/libTetherUnit_Interface.a: CMakeFiles/TetherUnit_Interface.dir/build.make
../lib/static/libTetherUnit_Interface.a: CMakeFiles/TetherUnit_Interface.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/jerluen/atu-crt/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX static library ../lib/static/libTetherUnit_Interface.a"
	$(CMAKE_COMMAND) -P CMakeFiles/TetherUnit_Interface.dir/cmake_clean_target.cmake
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/TetherUnit_Interface.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/TetherUnit_Interface.dir/build: ../lib/static/libTetherUnit_Interface.a

.PHONY : CMakeFiles/TetherUnit_Interface.dir/build

CMakeFiles/TetherUnit_Interface.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/TetherUnit_Interface.dir/cmake_clean.cmake
.PHONY : CMakeFiles/TetherUnit_Interface.dir/clean

CMakeFiles/TetherUnit_Interface.dir/depend:
	cd /home/jerluen/atu-crt/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/jerluen/atu-crt /home/jerluen/atu-crt /home/jerluen/atu-crt/build /home/jerluen/atu-crt/build /home/jerluen/atu-crt/build/CMakeFiles/TetherUnit_Interface.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/TetherUnit_Interface.dir/depend

