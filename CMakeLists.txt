############################ BASE ######################################

cmake_minimum_required (VERSION 3.14 FATAL_ERROR)
project(ERF-AMRWind CXX C)
list(APPEND CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/CMake")
include(CMakePackageConfigHelpers)

########################## OPTIONS #####################################

#General options for all executables in the project
set(ERF_HOME "${CMAKE_SOURCE_DIR}/../ERF" CACHE STRING "Path to ERF root directory")
set(AMRWIND_HOME "${CMAKE_SOURCE_DIR}/../amr-wind" CACHE STRING "Path to AMR-Wind root directory")
option(ERF_ENABLE_ALL_WARNINGS "Enable all compiler warnings" OFF)

#Options for performance
option(ERF_ENABLE_MPI "Enable MPI" OFF)
option(ERF_ENABLE_OPENMP "Enable OpenMP" OFF)
option(ERF_ENABLE_CUDA "Enable CUDA" OFF)
option(ERF_ENABLE_HIP  "Enable HIP" OFF)
option(ERF_ENABLE_SYCL "Enable SYCL" OFF)

#Options for C++
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

if(ERF_ENABLE_CUDA)
  enable_language(CUDA)
  if(CMAKE_CUDA_COMPILER_VERSION VERSION_LESS "11.0")
    message(FATAL_ERROR "Your nvcc version is ${CMAKE_CUDA_COMPILER_VERSION} which is unsupported."
      "Please use CUDA toolkit version 11.0 or newer.")
  endif()
endif()

########################### AMReX #####################################

set(CMAKE_PREFIX_PATH ${AMREX_DIR} ${CMAKE_PREFIX_PATH})
list(APPEND AMREX_COMPONENTS
  "3D" "PIC" "PARTICLES" "PDOUBLE" "DOUBLE" "LSOLVERS")
if (ERF_ENABLE_MPI)
  list(APPEND AMREX_COMPONENTS "MPI")
endif()
if (ERF_ENABLE_OPENMP)
  list(APPEND AMREX_COMPONENTS "OMP")
endif()
if (ERF_ENABLE_CUDA)
  list(APPEND AMREX_COMPONENTS "CUDA")
endif()
if (ERF_ENABLE_SYCL)
  list(APPEND AMREX_COMPONENTS "SYCL")
endif()
if (ERF_ENABLE_ROCM)
  list(APPEND AMREX_COMPONENTS "HIP")
endif()
if (ERF_ENABLE_HYPRE)
  list(APPEND AMREX_COMPONENTS "HYPRE")
endif()
if (ERF_ENABLE_TINY_PROFILE)
  list(APPEND AMREX_COMPONENTS "TINY_PROFILE")
endif()
separate_arguments(AMREX_COMPONENTS)
find_package(AMReX CONFIG REQUIRED
  COMPONENTS ${AMREX_COMPONENTS})
message(STATUS "Found AMReX = ${AMReX_DIR}")
set(FCOMPARE_EXE ${AMReX_DIR}/../../../bin/amrex_fcompare
  CACHE INTERNAL "Path to fcompare executable for regression tests")

########################### AMReX-Hydro #####################################

find_package(AMReX-Hydro CONFIG REQUIRED)
message(STATUS "Found AMReX-Hydro = ${AMReX-Hydro_DIR}")

########################### ERF #####################################

if(ERF_ENABLE_MPI)
  find_package(MPI REQUIRED)
endif()

# General information about machine, compiler, and build type
message(STATUS "ERF / AMR-Wind Driver Build Information:")
message(STATUS "CMAKE_SYSTEM_NAME = ${CMAKE_SYSTEM_NAME}")
message(STATUS "CMAKE_CXX_COMPILER_ID = ${CMAKE_CXX_COMPILER_ID}")
message(STATUS "CMAKE_CXX_COMPILER_VERSION = ${CMAKE_CXX_COMPILER_VERSION}")
message(STATUS "CMAKE_BUILD_TYPE = ${CMAKE_BUILD_TYPE}")

# Turn on rpath stuff
include(${CMAKE_SOURCE_DIR}/CMake/SetRpath.cmake)

#Build erf executables and intermediate object library and link to amrex library
add_subdirectory(Exec)
