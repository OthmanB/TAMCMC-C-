# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.15

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
CMAKE_COMMAND = /usr/local/Cellar/cmake/3.15.5/bin/cmake

# The command to remove a file.
RM = /usr/local/Cellar/cmake/3.15.5/bin/cmake -E remove -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /Users/obenomar/Documents/GitHub/TAMCMC-C

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /Users/obenomar/Documents/GitHub/TAMCMC-C/build

# Include any dependencies generated for this target.
include CMakeFiles/getstats.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/getstats.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/getstats.dir/flags.make

CMakeFiles/getstats.dir/tools/read_stats.cpp.o: CMakeFiles/getstats.dir/flags.make
CMakeFiles/getstats.dir/tools/read_stats.cpp.o: ../tools/read_stats.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/obenomar/Documents/GitHub/TAMCMC-C/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/getstats.dir/tools/read_stats.cpp.o"
	/Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/getstats.dir/tools/read_stats.cpp.o -c /Users/obenomar/Documents/GitHub/TAMCMC-C/tools/read_stats.cpp

CMakeFiles/getstats.dir/tools/read_stats.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/getstats.dir/tools/read_stats.cpp.i"
	/Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /Users/obenomar/Documents/GitHub/TAMCMC-C/tools/read_stats.cpp > CMakeFiles/getstats.dir/tools/read_stats.cpp.i

CMakeFiles/getstats.dir/tools/read_stats.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/getstats.dir/tools/read_stats.cpp.s"
	/Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /Users/obenomar/Documents/GitHub/TAMCMC-C/tools/read_stats.cpp -o CMakeFiles/getstats.dir/tools/read_stats.cpp.s

CMakeFiles/getstats.dir/tamcmc/sources/string_handler.cpp.o: CMakeFiles/getstats.dir/flags.make
CMakeFiles/getstats.dir/tamcmc/sources/string_handler.cpp.o: ../tamcmc/sources/string_handler.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/obenomar/Documents/GitHub/TAMCMC-C/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Building CXX object CMakeFiles/getstats.dir/tamcmc/sources/string_handler.cpp.o"
	/Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/getstats.dir/tamcmc/sources/string_handler.cpp.o -c /Users/obenomar/Documents/GitHub/TAMCMC-C/tamcmc/sources/string_handler.cpp

CMakeFiles/getstats.dir/tamcmc/sources/string_handler.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/getstats.dir/tamcmc/sources/string_handler.cpp.i"
	/Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /Users/obenomar/Documents/GitHub/TAMCMC-C/tamcmc/sources/string_handler.cpp > CMakeFiles/getstats.dir/tamcmc/sources/string_handler.cpp.i

CMakeFiles/getstats.dir/tamcmc/sources/string_handler.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/getstats.dir/tamcmc/sources/string_handler.cpp.s"
	/Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /Users/obenomar/Documents/GitHub/TAMCMC-C/tamcmc/sources/string_handler.cpp -o CMakeFiles/getstats.dir/tamcmc/sources/string_handler.cpp.s

# Object files for target getstats
getstats_OBJECTS = \
"CMakeFiles/getstats.dir/tools/read_stats.cpp.o" \
"CMakeFiles/getstats.dir/tamcmc/sources/string_handler.cpp.o"

# External object files for target getstats
getstats_EXTERNAL_OBJECTS =

getstats: CMakeFiles/getstats.dir/tools/read_stats.cpp.o
getstats: CMakeFiles/getstats.dir/tamcmc/sources/string_handler.cpp.o
getstats: CMakeFiles/getstats.dir/build.make
getstats: /usr/local/lib/libboost_system-mt.dylib
getstats: /usr/local/lib/libboost_filesystem-mt.dylib
getstats: /usr/local/lib/libboost_iostreams-mt.dylib
getstats: /usr/local/Cellar/gsl/2.6/lib/libgslcblas.dylib
getstats: /usr/local/Cellar/gsl/2.6/lib/libgsl.dylib
getstats: CMakeFiles/getstats.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/Users/obenomar/Documents/GitHub/TAMCMC-C/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_3) "Linking CXX executable getstats"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/getstats.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/getstats.dir/build: getstats

.PHONY : CMakeFiles/getstats.dir/build

CMakeFiles/getstats.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/getstats.dir/cmake_clean.cmake
.PHONY : CMakeFiles/getstats.dir/clean

CMakeFiles/getstats.dir/depend:
	cd /Users/obenomar/Documents/GitHub/TAMCMC-C/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /Users/obenomar/Documents/GitHub/TAMCMC-C /Users/obenomar/Documents/GitHub/TAMCMC-C /Users/obenomar/Documents/GitHub/TAMCMC-C/build /Users/obenomar/Documents/GitHub/TAMCMC-C/build /Users/obenomar/Documents/GitHub/TAMCMC-C/build/CMakeFiles/getstats.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/getstats.dir/depend

