# minimum version of CMake is 3.5
cmake_minimum_required(VERSION 3.5 FATAL_ERROR)

# define project
enable_language(Fortran)
project( knapsack
	 VERSION 1.0 )

# set testing vars
set(ENABLE_TESTS ON CACHE BOOL "")
set(ENABLE_GTEST OFF CACHE BOOL "")
set(ENABLE_FRUIT ON CACHE BOOL "")

# use BLT
include(blt/SetupBLT.cmake)

# required modules
include(GNUInstallDirs)
include(CheckIncludeFile)
include(CheckIncludeFiles)

# set custom compiler flags

blt_append_custom_compiler_flag(
	FLAGS_VAR BLT_PROD_FLAGS
	DEFAULT  "-ffree-form -fbackslash -fimplicit-none -fdefault-integer-8 -cpp -Wno-tabs"
	INTEL    "-free -i8 -implicitnone -fpp -O2 -assume bscc"
	PGI      "-Mfreeform -Mdclchk -Mstandard -Mallocatable=03 -Mbackslash")
list(APPEND CMAKE_Fortran_FLAGS ${BLT_PROD_FLAGS})

# add subdirectories
add_subdirectory(src)
add_subdirectory(tests)

