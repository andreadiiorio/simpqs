# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.14

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
CMAKE_COMMAND = /home/andysnake/Documenti/IDE/clion-2019.2/bin/cmake/linux/bin/cmake

# The command to remove a file.
RM = /home/andysnake/Documenti/IDE/clion-2019.2/bin/cmake/linux/bin/cmake -E remove -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /home/andysnake/Desktop/tenPrj

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/andysnake/Desktop/tenPrj/cmake-build-debug

# Include any dependencies generated for this target.
include CMakeFiles/SIMQS_Worker.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/SIMQS_Worker.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/SIMQS_Worker.dir/flags.make

CMakeFiles/SIMQS_Worker.dir/utils/utils.c.o: CMakeFiles/SIMQS_Worker.dir/flags.make
CMakeFiles/SIMQS_Worker.dir/utils/utils.c.o: ../utils/utils.c
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/andysnake/Desktop/tenPrj/cmake-build-debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building C object CMakeFiles/SIMQS_Worker.dir/utils/utils.c.o"
	/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -o CMakeFiles/SIMQS_Worker.dir/utils/utils.c.o   -c /home/andysnake/Desktop/tenPrj/utils/utils.c

CMakeFiles/SIMQS_Worker.dir/utils/utils.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/SIMQS_Worker.dir/utils/utils.c.i"
	/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -E /home/andysnake/Desktop/tenPrj/utils/utils.c > CMakeFiles/SIMQS_Worker.dir/utils/utils.c.i

CMakeFiles/SIMQS_Worker.dir/utils/utils.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/SIMQS_Worker.dir/utils/utils.c.s"
	/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -S /home/andysnake/Desktop/tenPrj/utils/utils.c -o CMakeFiles/SIMQS_Worker.dir/utils/utils.c.s

CMakeFiles/SIMQS_Worker.dir/utils/gmp_patch.c.o: CMakeFiles/SIMQS_Worker.dir/flags.make
CMakeFiles/SIMQS_Worker.dir/utils/gmp_patch.c.o: ../utils/gmp_patch.c
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/andysnake/Desktop/tenPrj/cmake-build-debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Building C object CMakeFiles/SIMQS_Worker.dir/utils/gmp_patch.c.o"
	/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -o CMakeFiles/SIMQS_Worker.dir/utils/gmp_patch.c.o   -c /home/andysnake/Desktop/tenPrj/utils/gmp_patch.c

CMakeFiles/SIMQS_Worker.dir/utils/gmp_patch.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/SIMQS_Worker.dir/utils/gmp_patch.c.i"
	/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -E /home/andysnake/Desktop/tenPrj/utils/gmp_patch.c > CMakeFiles/SIMQS_Worker.dir/utils/gmp_patch.c.i

CMakeFiles/SIMQS_Worker.dir/utils/gmp_patch.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/SIMQS_Worker.dir/utils/gmp_patch.c.s"
	/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -S /home/andysnake/Desktop/tenPrj/utils/gmp_patch.c -o CMakeFiles/SIMQS_Worker.dir/utils/gmp_patch.c.s

CMakeFiles/SIMQS_Worker.dir/worker/sievingSIMPQS.c.o: CMakeFiles/SIMQS_Worker.dir/flags.make
CMakeFiles/SIMQS_Worker.dir/worker/sievingSIMPQS.c.o: ../worker/sievingSIMPQS.c
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/andysnake/Desktop/tenPrj/cmake-build-debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_3) "Building C object CMakeFiles/SIMQS_Worker.dir/worker/sievingSIMPQS.c.o"
	/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -o CMakeFiles/SIMQS_Worker.dir/worker/sievingSIMPQS.c.o   -c /home/andysnake/Desktop/tenPrj/worker/sievingSIMPQS.c

CMakeFiles/SIMQS_Worker.dir/worker/sievingSIMPQS.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/SIMQS_Worker.dir/worker/sievingSIMPQS.c.i"
	/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -E /home/andysnake/Desktop/tenPrj/worker/sievingSIMPQS.c > CMakeFiles/SIMQS_Worker.dir/worker/sievingSIMPQS.c.i

CMakeFiles/SIMQS_Worker.dir/worker/sievingSIMPQS.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/SIMQS_Worker.dir/worker/sievingSIMPQS.c.s"
	/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -S /home/andysnake/Desktop/tenPrj/worker/sievingSIMPQS.c -o CMakeFiles/SIMQS_Worker.dir/worker/sievingSIMPQS.c.s

CMakeFiles/SIMQS_Worker.dir/factorization/factorizerQuick.c.o: CMakeFiles/SIMQS_Worker.dir/flags.make
CMakeFiles/SIMQS_Worker.dir/factorization/factorizerQuick.c.o: ../factorization/factorizerQuick.c
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/andysnake/Desktop/tenPrj/cmake-build-debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_4) "Building C object CMakeFiles/SIMQS_Worker.dir/factorization/factorizerQuick.c.o"
	/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -o CMakeFiles/SIMQS_Worker.dir/factorization/factorizerQuick.c.o   -c /home/andysnake/Desktop/tenPrj/factorization/factorizerQuick.c

CMakeFiles/SIMQS_Worker.dir/factorization/factorizerQuick.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/SIMQS_Worker.dir/factorization/factorizerQuick.c.i"
	/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -E /home/andysnake/Desktop/tenPrj/factorization/factorizerQuick.c > CMakeFiles/SIMQS_Worker.dir/factorization/factorizerQuick.c.i

CMakeFiles/SIMQS_Worker.dir/factorization/factorizerQuick.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/SIMQS_Worker.dir/factorization/factorizerQuick.c.s"
	/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -S /home/andysnake/Desktop/tenPrj/factorization/factorizerQuick.c -o CMakeFiles/SIMQS_Worker.dir/factorization/factorizerQuick.c.s

# Object files for target SIMQS_Worker
SIMQS_Worker_OBJECTS = \
"CMakeFiles/SIMQS_Worker.dir/utils/utils.c.o" \
"CMakeFiles/SIMQS_Worker.dir/utils/gmp_patch.c.o" \
"CMakeFiles/SIMQS_Worker.dir/worker/sievingSIMPQS.c.o" \
"CMakeFiles/SIMQS_Worker.dir/factorization/factorizerQuick.c.o"

# External object files for target SIMQS_Worker
SIMQS_Worker_EXTERNAL_OBJECTS =

SIMQS_Worker: CMakeFiles/SIMQS_Worker.dir/utils/utils.c.o
SIMQS_Worker: CMakeFiles/SIMQS_Worker.dir/utils/gmp_patch.c.o
SIMQS_Worker: CMakeFiles/SIMQS_Worker.dir/worker/sievingSIMPQS.c.o
SIMQS_Worker: CMakeFiles/SIMQS_Worker.dir/factorization/factorizerQuick.c.o
SIMQS_Worker: CMakeFiles/SIMQS_Worker.dir/build.make
SIMQS_Worker: CMakeFiles/SIMQS_Worker.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/andysnake/Desktop/tenPrj/cmake-build-debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_5) "Linking C executable SIMQS_Worker"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/SIMQS_Worker.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/SIMQS_Worker.dir/build: SIMQS_Worker

.PHONY : CMakeFiles/SIMQS_Worker.dir/build

CMakeFiles/SIMQS_Worker.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/SIMQS_Worker.dir/cmake_clean.cmake
.PHONY : CMakeFiles/SIMQS_Worker.dir/clean

CMakeFiles/SIMQS_Worker.dir/depend:
	cd /home/andysnake/Desktop/tenPrj/cmake-build-debug && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/andysnake/Desktop/tenPrj /home/andysnake/Desktop/tenPrj /home/andysnake/Desktop/tenPrj/cmake-build-debug /home/andysnake/Desktop/tenPrj/cmake-build-debug /home/andysnake/Desktop/tenPrj/cmake-build-debug/CMakeFiles/SIMQS_Worker.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/SIMQS_Worker.dir/depend

