# minimum version of CMake is 3.5
cmake_minimum_required(VERSION 3.5 FATAL_ERROR)

# define project
enable_language(Fortran)
project( knapsack
	 VERSION 1.0 )
	 
# Tell it where to find libs
set (LIBPATH "/usr/local")
set (INCPATH "/usr/local")

# required for unit testing
include(CTest)
enable_testing()

# required modules
include(GNUInstallDirs)
include(CheckIncludeFile)
include(CheckIncludeFiles)

find_package(MPI REQUIRED)
include_directories(${MPI_INCLUDE_PATH})

# we need at least F03
if(CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
    set(dialect "-ffree-form -fbackslash -fimplicit-none -fdefault-integer-8 -Wno-tabs")
    set(bounds "-fbounds-check")
endif()
if(CMAKE_Fortran_COMPILER_ID MATCHES "Intel")
    set(dialect "-stand f03 -free -i8 -implicitnone -backslash")
    set(bounds "-check bounds")
endif()
if(CMAKE_Fortran_COMPILER_ID MATCHES "PGI")
    set(dialect "-Mfreeform -Mdclchk -Mstandard -Mallocatable=03 -Mbackslash")
    set(bounds "-C")
endif()

list(APPEND CMAKE_Fortran_FLAGS_DEBUG ${bounds})
list(APPEND CMAKE_Fortran_FLAGS ${dialect})

# add subdirectories
add_subdirectory(src)
add_subdirectory(tests)
