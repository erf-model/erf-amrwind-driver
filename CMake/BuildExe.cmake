function(target_link_libraries_system target visibility)
  set(libs ${ARGN})
  foreach(lib ${libs})
    get_target_property(lib_include_dirs ${lib} INTERFACE_INCLUDE_DIRECTORIES)
    target_include_directories(${target} SYSTEM ${visibility} ${lib_include_dirs})
    target_link_libraries(${target} ${visibility} ${lib})
  endforeach(lib)
endfunction(target_link_libraries_system)


function(build_erf_lib_amrw erf_lib_name)
  find_package(AMR-Wind REQUIRED)
  message(STATUS "Found AMR-Wind = ${AMR_WIND_INCLUDE_DIR}")
  message(STATUS "Found AMR-Wind = ${AMR_WIND_LIBRARY_DIR}")
  target_link_libraries_system(${erf_lib_name} PUBLIC
    AMR-Wind::amrwind_api)
    target_link_libraries_system(${erf_lib_name} PUBLIC
    AMR-Wind::buildInfoamrwind_obj)
endfunction(build_erf_lib_amrw)


function(build_erf_lib_erf erf_lib_name)
  find_package(ERF REQUIRED)
  message(STATUS "Found ERF = ${ERF_INCLUDE_DIR}")
  message(STATUS "Found ERF = ${ERF_LIBRARY_DIR}")
  target_link_libraries_system(${erf_lib_name} PUBLIC
    ERF::erf_api)
endfunction(build_erf_lib_erf)


function(build_erf_lib_amrex erf_lib_name)
  #Link to amrex library
  target_link_libraries_system(${erf_lib_name} PUBLIC AMReX::amrex)
  target_link_libraries_system(${erf_lib_name} PUBLIC AMReX-Hydro::amrex_hydro_api)
  if(ERF_ENABLE_CUDA)
    set(pctargets "${erf_lib_name}")
    foreach(tgt IN LISTS pctargets)
      get_target_property(ERF_SOURCES ${tgt} SOURCES)
      list(FILTER ERF_SOURCES INCLUDE REGEX "\\.cpp")
      set_source_files_properties(${ERF_SOURCES} PROPERTIES LANGUAGE CUDA)
      message(STATUS "setting cuda for ${ERF_SOURCES}")
    endforeach()
    set_target_properties(
    ${erf_lib_name} PROPERTIES
    LANGUAGE CUDA
    CUDA_SEPARABLE_COMPILATION ON
    CUDA_RESOLVE_DEVICE_SYMBOLS ON)
  endif()
endfunction(build_erf_lib_amrex)


function(build_erf_lib_wrapper erf_lib_name)
  set(SRC_DIR ${CMAKE_SOURCE_DIR}/Source)
  set(BIN_DIR ${CMAKE_BINARY_DIR}/Source/${erf_lib_name})
  set(ERF_SRC_DIR  ${ERF_HOME}/Source)
  set(AMRW_SRC_DIR  ${AMRWIND_HOME}/amr-wind)

  include(${CMAKE_SOURCE_DIR}/CMake/SetERFCompileFlags.cmake)
  set_erf_compile_flags(${erf_lib_name})

  target_compile_definitions(${erf_lib_name} PUBLIC ERF_USE_MOISTURE)

  target_compile_definitions(${erf_lib_name} PUBLIC ERF_MB_EXTERN)
  if(ERF_ENABLE_MULTIBLOCK)
    target_sources(${erf_lib_name} PRIVATE
                   ${SRC_DIR}/main.cpp
                   ${SRC_DIR}/amrwind_Evolve_MB.cpp
                   ${SRC_DIR}/ERF_Evolve_MB.cpp
                   ${SRC_DIR}/MultiBlock/MultiBlockContainer.cpp
                   ${SRC_DIR}/wind_energy/ABLReadERF.cpp)
    target_compile_definitions(${erf_lib_name} PUBLIC ERF_USE_MULTIBLOCK)
    target_include_directories(${erf_lib_name} PRIVATE
                                ${SRC_DIR}/MultiBlock
                                ${SRC_DIR}/wind_energy)
  endif()

  build_erf_lib_amrw(${erf_lib_name})
  build_erf_lib_erf(${erf_lib_name})
  build_erf_lib_amrex(${erf_lib_name})

  #Define what we want to be installed during a make install
  install(TARGETS ${erf_lib_name}
          RUNTIME DESTINATION bin
          ARCHIVE DESTINATION lib
          LIBRARY DESTINATION lib)
endfunction(build_erf_lib_wrapper)


function(build_erf_exe erf_exe_name)
  set(SRC_DIR ${CMAKE_SOURCE_DIR}/Source)
  set(ERF_SRC_DIR  ${ERF_HOME}/Source)

  target_link_libraries(${erf_exe_name} PRIVATE ${erf_lib_name} AMReX-Hydro::amrex_hydro_api)
  target_link_libraries(${erf_exe_name}  PUBLIC ${erf_lib_name})
  include(${CMAKE_SOURCE_DIR}/CMake/SetERFCompileFlags.cmake)
  set_erf_compile_flags(${erf_exe_name})

  if(ERF_ENABLE_CUDA)
    set(pctargets "${erf_exe_name}")
    foreach(tgt IN LISTS pctargets)
      get_target_property(ERF_SOURCES ${tgt} SOURCES)
      list(FILTER ERF_SOURCES INCLUDE REGEX "\\.cpp")
      set_source_files_properties(${ERF_SOURCES} PROPERTIES LANGUAGE CUDA)
      message(STATUS "setting cuda for ${ERF_SOURCES}")
    endforeach()
    set_target_properties(
    ${erf_exe_name} PROPERTIES
    LANGUAGE CUDA
    CUDA_SEPARABLE_COMPILATION ON
    CUDA_RESOLVE_DEVICE_SYMBOLS ON)
  endif()

  install(TARGETS ${erf_exe_name}
          RUNTIME DESTINATION bin
          ARCHIVE DESTINATION lib
          LIBRARY DESTINATION lib)
endfunction()
