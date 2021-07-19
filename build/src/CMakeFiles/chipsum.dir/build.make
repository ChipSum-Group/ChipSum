# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.16

# Delete rule output on recipe failure.
.DELETE_ON_ERROR:


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

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /home/lky/code/git/ChipSum

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/lky/code/git/ChipSum/build

# Include any dependencies generated for this target.
include src/CMakeFiles/chipsum.dir/depend.make

# Include the progress variables for this target.
include src/CMakeFiles/chipsum.dir/progress.make

# Include the compile flags for this target's objects.
include src/CMakeFiles/chipsum.dir/flags.make

src/CMakeFiles/chipsum.dir/numeric/operator.cpp.o: src/CMakeFiles/chipsum.dir/flags.make
src/CMakeFiles/chipsum.dir/numeric/operator.cpp.o: ../src/numeric/operator.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/lky/code/git/ChipSum/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object src/CMakeFiles/chipsum.dir/numeric/operator.cpp.o"
	cd /home/lky/code/git/ChipSum/build/src && /usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/chipsum.dir/numeric/operator.cpp.o -c /home/lky/code/git/ChipSum/src/numeric/operator.cpp

src/CMakeFiles/chipsum.dir/numeric/operator.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/chipsum.dir/numeric/operator.cpp.i"
	cd /home/lky/code/git/ChipSum/build/src && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/lky/code/git/ChipSum/src/numeric/operator.cpp > CMakeFiles/chipsum.dir/numeric/operator.cpp.i

src/CMakeFiles/chipsum.dir/numeric/operator.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/chipsum.dir/numeric/operator.cpp.s"
	cd /home/lky/code/git/ChipSum/build/src && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/lky/code/git/ChipSum/src/numeric/operator.cpp -o CMakeFiles/chipsum.dir/numeric/operator.cpp.s

src/CMakeFiles/chipsum.dir/numeric/vector.cpp.o: src/CMakeFiles/chipsum.dir/flags.make
src/CMakeFiles/chipsum.dir/numeric/vector.cpp.o: ../src/numeric/vector.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/lky/code/git/ChipSum/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Building CXX object src/CMakeFiles/chipsum.dir/numeric/vector.cpp.o"
	cd /home/lky/code/git/ChipSum/build/src && /usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/chipsum.dir/numeric/vector.cpp.o -c /home/lky/code/git/ChipSum/src/numeric/vector.cpp

src/CMakeFiles/chipsum.dir/numeric/vector.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/chipsum.dir/numeric/vector.cpp.i"
	cd /home/lky/code/git/ChipSum/build/src && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/lky/code/git/ChipSum/src/numeric/vector.cpp > CMakeFiles/chipsum.dir/numeric/vector.cpp.i

src/CMakeFiles/chipsum.dir/numeric/vector.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/chipsum.dir/numeric/vector.cpp.s"
	cd /home/lky/code/git/ChipSum/build/src && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/lky/code/git/ChipSum/src/numeric/vector.cpp -o CMakeFiles/chipsum.dir/numeric/vector.cpp.s

src/CMakeFiles/chipsum.dir/numeric/impl/vector_impl.cpp.o: src/CMakeFiles/chipsum.dir/flags.make
src/CMakeFiles/chipsum.dir/numeric/impl/vector_impl.cpp.o: ../src/numeric/impl/vector_impl.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/lky/code/git/ChipSum/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_3) "Building CXX object src/CMakeFiles/chipsum.dir/numeric/impl/vector_impl.cpp.o"
	cd /home/lky/code/git/ChipSum/build/src && /usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/chipsum.dir/numeric/impl/vector_impl.cpp.o -c /home/lky/code/git/ChipSum/src/numeric/impl/vector_impl.cpp

src/CMakeFiles/chipsum.dir/numeric/impl/vector_impl.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/chipsum.dir/numeric/impl/vector_impl.cpp.i"
	cd /home/lky/code/git/ChipSum/build/src && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/lky/code/git/ChipSum/src/numeric/impl/vector_impl.cpp > CMakeFiles/chipsum.dir/numeric/impl/vector_impl.cpp.i

src/CMakeFiles/chipsum.dir/numeric/impl/vector_impl.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/chipsum.dir/numeric/impl/vector_impl.cpp.s"
	cd /home/lky/code/git/ChipSum/build/src && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/lky/code/git/ChipSum/src/numeric/impl/vector_impl.cpp -o CMakeFiles/chipsum.dir/numeric/impl/vector_impl.cpp.s

src/CMakeFiles/chipsum.dir/backend/backend.cpp.o: src/CMakeFiles/chipsum.dir/flags.make
src/CMakeFiles/chipsum.dir/backend/backend.cpp.o: ../src/backend/backend.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/lky/code/git/ChipSum/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_4) "Building CXX object src/CMakeFiles/chipsum.dir/backend/backend.cpp.o"
	cd /home/lky/code/git/ChipSum/build/src && /usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/chipsum.dir/backend/backend.cpp.o -c /home/lky/code/git/ChipSum/src/backend/backend.cpp

src/CMakeFiles/chipsum.dir/backend/backend.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/chipsum.dir/backend/backend.cpp.i"
	cd /home/lky/code/git/ChipSum/build/src && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/lky/code/git/ChipSum/src/backend/backend.cpp > CMakeFiles/chipsum.dir/backend/backend.cpp.i

src/CMakeFiles/chipsum.dir/backend/backend.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/chipsum.dir/backend/backend.cpp.s"
	cd /home/lky/code/git/ChipSum/build/src && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/lky/code/git/ChipSum/src/backend/backend.cpp -o CMakeFiles/chipsum.dir/backend/backend.cpp.s

# Object files for target chipsum
chipsum_OBJECTS = \
"CMakeFiles/chipsum.dir/numeric/operator.cpp.o" \
"CMakeFiles/chipsum.dir/numeric/vector.cpp.o" \
"CMakeFiles/chipsum.dir/numeric/impl/vector_impl.cpp.o" \
"CMakeFiles/chipsum.dir/backend/backend.cpp.o"

# External object files for target chipsum
chipsum_EXTERNAL_OBJECTS =

src/libchipsum.a: src/CMakeFiles/chipsum.dir/numeric/operator.cpp.o
src/libchipsum.a: src/CMakeFiles/chipsum.dir/numeric/vector.cpp.o
src/libchipsum.a: src/CMakeFiles/chipsum.dir/numeric/impl/vector_impl.cpp.o
src/libchipsum.a: src/CMakeFiles/chipsum.dir/backend/backend.cpp.o
src/libchipsum.a: src/CMakeFiles/chipsum.dir/build.make
src/libchipsum.a: src/CMakeFiles/chipsum.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/lky/code/git/ChipSum/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_5) "Linking CXX static library libchipsum.a"
	cd /home/lky/code/git/ChipSum/build/src && $(CMAKE_COMMAND) -P CMakeFiles/chipsum.dir/cmake_clean_target.cmake
	cd /home/lky/code/git/ChipSum/build/src && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/chipsum.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
src/CMakeFiles/chipsum.dir/build: src/libchipsum.a

.PHONY : src/CMakeFiles/chipsum.dir/build

src/CMakeFiles/chipsum.dir/clean:
	cd /home/lky/code/git/ChipSum/build/src && $(CMAKE_COMMAND) -P CMakeFiles/chipsum.dir/cmake_clean.cmake
.PHONY : src/CMakeFiles/chipsum.dir/clean

src/CMakeFiles/chipsum.dir/depend:
	cd /home/lky/code/git/ChipSum/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/lky/code/git/ChipSum /home/lky/code/git/ChipSum/src /home/lky/code/git/ChipSum/build /home/lky/code/git/ChipSum/build/src /home/lky/code/git/ChipSum/build/src/CMakeFiles/chipsum.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : src/CMakeFiles/chipsum.dir/depend

