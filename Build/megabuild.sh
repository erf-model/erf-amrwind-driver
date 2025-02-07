#!/bin/bash

#set -o errexit

export TOP=$(pwd)
export TOP_MOD=$(pwd)/erf-amrwind-driver/Submodules

export CXXFLAGS=-fPIC

### Init and update submodules
cd ${TOP}/erf-amrwind-driver
git submodule update --init --depth=1

### Build AMReX
cd ${TOP_MOD}/amrex
cmake -DBUILD_SHARED_LIBS=OFF \
      -DAMReX_EB=ON -DAMReX_PIC=YES \
      -B ${TOP}/amrex-build -DCMAKE_INSTALL_PREFIX=${TOP}/amrex-install -S ${TOP_MOD}/amrex \
      -DCMAKE_BUILD_TYPE:STRING=RELEASE
cd ${TOP}/amrex-build
make -j16
make install

### Build AMReX-Hydro
cd ${TOP_MOD}/AMReX-Hydro
cmake -DBUILD_SHARED_LIBS=OFF \
      -DAMReX_EB=OFF -DHYDRO_EB=OFF \
      -B ${TOP}/AMReX-Hydro-build -DCMAKE_INSTALL_PREFIX=${TOP}/AMReX-Hydro-install -S ${TOP_MOD}/AMReX-Hydro \
      -DAMReX_ROOT=${TOP}/amrex-install/lib/cmake/AMReX/ \
      -DCMAKE_PREFIX_PATH==${TOP}/amrex-install/lib/cmake/AMReX \
      -DCMAKE_BUILD_TYPE:STRING=RELEASE
cd ${TOP}/AMReX-Hydro-build
make -j16
make install

### Build ERF
cmake -DBUILD_SHARED_LIBS=OFF \
      -DCMAKE_INSTALL_PREFIX:PATH=${TOP}/ERF-install \
      -DCMAKE_CXX_COMPILER:STRING=mpicxx \
      -DCMAKE_C_COMPILER:STRING=mpicc \
      -DCMAKE_Fortran_COMPILER:STRING=mpifort \
      -DCMAKE_BUILD_TYPE:STRING=Release \
      -DERF_ENABLE_MPI:BOOL=ON \
      -DERF_ENABLE_MULTIBLOCK:BOOL=ON \
      -DERF_USE_INTERNAL_AMREX:BOOL=OFF \
      -DERF_ENABLE_TESTS:BOOL=OFF \
      -DCMAKE_PREFIX_PATH=${TOP}/amrex-install/lib/cmake/AMReX \
      -DERF_ENABLE_FCOMPARE:BOOL=ON \
      -DERF_ENABLE_DOCUMENTATION:BOOL=OFF \
      -DCMAKE_EXPORT_COMPILE_COMMANDS:BOOL=ON \
      -S ${TOP_MOD}/ERF -B ${TOP}/ERF-build
cd ${TOP}/ERF-build
make -j16
make install

### Build AMR-Wind
cd ${TOP_MOD}/amr-wind
cmake -DBUILD_SHARED_LIBS=ON \
      -B ${TOP}/amr-wind-build \
      -DAMR_WIND_USE_INTERNAL_AMREX=OFF -DAMR_WIND_USE_INTERNAL_AMREX_HYDRO=OFF \
      -DCMAKE_INSTALL_PREFIX=${TOP}/amr-wind-install \
      -DAMR_WIND_ENABLE_TESTS:BOOL=OFF \
      -DAMR_WIND_ENABLE_UNIT_TESTS:BOOL=OFF \
      -S ${TOP_MOD}/amr-wind \
      -DCMAKE_PREFIX_PATH="${TOP}/amrex-install/lib/cmake/AMReX;${TOP}/AMReX-Hydro-install/lib/cmake/AMReX-Hydro" \
      -DCMAKE_BUILD_TYPE:STRING=RELEASE \
      -DERF_AMR_WIND_MULTIBLOCK:BOOL=ON
cd ${TOP}/amr-wind-build
make -j16
make install

### Build the coupling driver
cd ${TOP}
cmake -DCMAKE_INSTALL_PREFIX:PATH=${TOP}/erf-amrwind-driver-install \
      -DCMAKE_CXX_COMPILER:STRING=mpicxx \
      -DCMAKE_C_COMPILER:STRING=mpicc \
      -DCMAKE_Fortran_COMPILER:STRING=mpifort \
      -DCMAKE_BUILD_TYPE:STRING=Release \
      -DAMRWIND_HOME:STRING=${TOP_MOD}/amr-wind \
      -DERF_HOME:STRING=${TOP_MOD}/ERF \
      -DERF_ENABLE_MPI:BOOL=ON \
      -DCMAKE_EXPORT_COMPILE_COMMANDS:BOOL=ON \
      -DCMAKE_PREFIX_PATH="${TOP}/amrex-install/lib/cmake/AMReX;${TOP}/AMReX-Hydro-install/lib/cmake/AMReX-Hydro;${TOP}/amr-wind-install/lib/cmake/AMR-Wind;${TOP}/ERF-install/lib/cmake/ERF" \
      -S ${TOP}/erf-amrwind-driver -B ${TOP}/erf-amrwind-driver-build
cd ${TOP}/erf-amrwind-driver-build
make -j16
make install

cd ${TOP}
