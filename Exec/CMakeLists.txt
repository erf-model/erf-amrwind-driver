set(erf_lib_name erf-amr-wind-driver_srclib)
add_library(${erf_lib_name} OBJECT)
include(${CMAKE_SOURCE_DIR}/CMake/BuildExe.cmake)
build_erf_lib_wrapper(${erf_lib_name})

add_subdirectory(ABL)
add_subdirectory(CouetteFlow)
add_subdirectory(PoiseuilleFlow)
