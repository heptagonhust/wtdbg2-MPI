# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.10

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
CMAKE_SOURCE_DIR = /mnt/c/Users/Alisa/Documents/MyFile/code/ASC/wtdbg2-GPU/src

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /mnt/c/Users/Alisa/Documents/MyFile/code/ASC/wtdbg2-GPU/src

# Include any dependencies generated for this target.
include CMakeFiles/wtdbg-cns.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/wtdbg-cns.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/wtdbg-cns.dir/flags.make

CMakeFiles/wtdbg-cns.dir/wtdbg-cns.cpp.o: CMakeFiles/wtdbg-cns.dir/flags.make
CMakeFiles/wtdbg-cns.dir/wtdbg-cns.cpp.o: wtdbg-cns.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/mnt/c/Users/Alisa/Documents/MyFile/code/ASC/wtdbg2-GPU/src/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/wtdbg-cns.dir/wtdbg-cns.cpp.o"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/wtdbg-cns.dir/wtdbg-cns.cpp.o -c /mnt/c/Users/Alisa/Documents/MyFile/code/ASC/wtdbg2-GPU/src/wtdbg-cns.cpp

CMakeFiles/wtdbg-cns.dir/wtdbg-cns.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/wtdbg-cns.dir/wtdbg-cns.cpp.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /mnt/c/Users/Alisa/Documents/MyFile/code/ASC/wtdbg2-GPU/src/wtdbg-cns.cpp > CMakeFiles/wtdbg-cns.dir/wtdbg-cns.cpp.i

CMakeFiles/wtdbg-cns.dir/wtdbg-cns.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/wtdbg-cns.dir/wtdbg-cns.cpp.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /mnt/c/Users/Alisa/Documents/MyFile/code/ASC/wtdbg2-GPU/src/wtdbg-cns.cpp -o CMakeFiles/wtdbg-cns.dir/wtdbg-cns.cpp.s

CMakeFiles/wtdbg-cns.dir/wtdbg-cns.cpp.o.requires:

.PHONY : CMakeFiles/wtdbg-cns.dir/wtdbg-cns.cpp.o.requires

CMakeFiles/wtdbg-cns.dir/wtdbg-cns.cpp.o.provides: CMakeFiles/wtdbg-cns.dir/wtdbg-cns.cpp.o.requires
	$(MAKE) -f CMakeFiles/wtdbg-cns.dir/build.make CMakeFiles/wtdbg-cns.dir/wtdbg-cns.cpp.o.provides.build
.PHONY : CMakeFiles/wtdbg-cns.dir/wtdbg-cns.cpp.o.provides

CMakeFiles/wtdbg-cns.dir/wtdbg-cns.cpp.o.provides.build: CMakeFiles/wtdbg-cns.dir/wtdbg-cns.cpp.o


# Object files for target wtdbg-cns
wtdbg__cns_OBJECTS = \
"CMakeFiles/wtdbg-cns.dir/wtdbg-cns.cpp.o"

# External object files for target wtdbg-cns
wtdbg__cns_EXTERNAL_OBJECTS =

wtdbg-cns: CMakeFiles/wtdbg-cns.dir/wtdbg-cns.cpp.o
wtdbg-cns: CMakeFiles/wtdbg-cns.dir/build.make
wtdbg-cns: libksw.a
wtdbg-cns: CMakeFiles/wtdbg-cns.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/mnt/c/Users/Alisa/Documents/MyFile/code/ASC/wtdbg2-GPU/src/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable wtdbg-cns"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/wtdbg-cns.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/wtdbg-cns.dir/build: wtdbg-cns

.PHONY : CMakeFiles/wtdbg-cns.dir/build

CMakeFiles/wtdbg-cns.dir/requires: CMakeFiles/wtdbg-cns.dir/wtdbg-cns.cpp.o.requires

.PHONY : CMakeFiles/wtdbg-cns.dir/requires

CMakeFiles/wtdbg-cns.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/wtdbg-cns.dir/cmake_clean.cmake
.PHONY : CMakeFiles/wtdbg-cns.dir/clean

CMakeFiles/wtdbg-cns.dir/depend:
	cd /mnt/c/Users/Alisa/Documents/MyFile/code/ASC/wtdbg2-GPU/src && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /mnt/c/Users/Alisa/Documents/MyFile/code/ASC/wtdbg2-GPU/src /mnt/c/Users/Alisa/Documents/MyFile/code/ASC/wtdbg2-GPU/src /mnt/c/Users/Alisa/Documents/MyFile/code/ASC/wtdbg2-GPU/src /mnt/c/Users/Alisa/Documents/MyFile/code/ASC/wtdbg2-GPU/src /mnt/c/Users/Alisa/Documents/MyFile/code/ASC/wtdbg2-GPU/src/CMakeFiles/wtdbg-cns.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/wtdbg-cns.dir/depend

