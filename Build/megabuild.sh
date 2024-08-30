#!/bin/bash

#set -o errexit

mkdir test
cd test
export TOP=$(pwd)

git clone git@github.com:mukul1992/ERF.git --branch preserve-amrw-coupling --single-branch # <folder>
git clone git@github.com:mukul1992/amr-wind.git --branch preserve-erf-coupling --single-branch # <folder> 
git clone git@github.com:erf-model/erf-amrwind-driver.git

cd ${TOP}/ERF
git submodule update --init

cd ${TOP}/amr-wind
git submodule update --init

cd ${TOP}
git clone amr-wind/submods/amrex
git clone amr-wind/submods/AMReX-Hydro

cd ${TOP}/amrex
cmake -DAMReX_EB=OFF -DAMReX_PIC=YES -B ${TOP}/amrex/build -DCMAKE_INSTALL_PREFIX=${TOP}/amrex/install -S . -DCMAKE_BUILD_TYPE:STRING=RELEASE
cd ${TOP}/amrex/build
make -j16
make install

cd ${TOP}/AMReX-Hydro
cmake -DAMReX_EB=OFF -DHYDRO_EB=OFF -B ${TOP}/AMReX-Hydro/build -DCMAKE_INSTALL_PREFIX=${TOP}/AMReX-Hydro/install -S . -DAMReX_ROOT=${TOP}/amrex/install/lib/cmake/AMReX/ -DCMAKE_PREFIX_PATH==${TOP}/amrex/install/lib/cmake/AMReX -DCMAKE_BUILD_TYPE:STRING=RELEASE
cd ${TOP}/AMReX-Hydro/build
make -j16
make install

cd ${TOP}/ERF
git apply ${TOP}/../external_amrex_erf_fixes.patch
head -n13 Build/cmake_multiblock.sh > ../erf_cmake_multiblock.sh
echo "       -DERF_USE_INTERNAL_AMREX:BOOL=OFF \\" >> ../erf_cmake_multiblock.sh 
echo "       -DERF_ENABLE_TESTS:BOOL=OFF \\" >> ../erf_cmake_multiblock.sh 
echo "       -DCMAKE_PREFIX_PATH=${TOP}/amrex/install/lib/cmake/AMReX \\" >> ../erf_cmake_multiblock.sh
tail -n 4 Build/cmake_multiblock.sh >> ../erf_cmake_multiblock.sh
chmod +x ../erf_cmake_multiblock.sh
#../erf_cmake_multiblock.sh
#cd build
#make install

cmake -DCMAKE_INSTALL_PREFIX:PATH=${TOP}/ERF/install \
      -DCMAKE_CXX_COMPILER:STRING=mpicxx \
      -DCMAKE_C_COMPILER:STRING=mpicc \
      -DCMAKE_Fortran_COMPILER:STRING=mpifort \
      -DMPIEXEC_PREFLAGS:STRING=--oversubscribe \
      -DCMAKE_BUILD_TYPE:STRING=Release \
      -DERF_DIM:STRING=3 \
      -DERF_ENABLE_MPI:BOOL=ON \
      -DERF_ENABLE_MULTIBLOCK:BOOL=ON \
       -DERF_USE_INTERNAL_AMREX:BOOL=OFF \
       -DERF_ENABLE_TESTS:BOOL=OFF \
       -DCMAKE_PREFIX_PATH=${TOP}/amrex/install/lib/cmake/AMReX \
      -DERF_ENABLE_FCOMPARE:BOOL=ON \
      -DERF_ENABLE_DOCUMENTATION:BOOL=OFF \
      -DCMAKE_EXPORT_COMPILE_COMMANDS:BOOL=ON \
      -S ${TOP}/ERF -B ${TOP}/ERF/build
cd build
make -j16
make install

cd ${TOP}
cmake -B ${TOP}/amr-wind-build -DAMR_WIND_USE_INTERNAL_AMREX=OFF -DAMR_WIND_USE_INTERNAL_AMREX_HYDRO=OFF -DCMAKE_INSTALL_PREFIX=${TOP}/amr-wind-install -S ${TOP}/amr-wind -DCMAKE_PREFIX_PATH="${TOP}/amrex/install/lib/cmake/AMReX;${TOP}/AMReX-Hydro/install/lib/cmake/AMReX-Hydro" -DCMAKE_BUILD_TYPE:STRING=RELEASE -DERF_ENABLE_MULTIBLOCK:BOOL=ON
cd amr-wind-build
make -j16
make install

cd ${TOP}
cmake -B ${TOP}/amr-wind-build-internal -DAMR_WIND_USE_INTERNAL_AMREX=ON -DAMR_WIND_USE_INTERNAL_AMREX_HYDRO=ON -DCMAKE_INSTALL_PREFIX=${TOP}/amr-wind-install-internal -S ${TOP}/amr-wind -DCMAKE_BUILD_TYPE:STRING=RELEASE -DERF_ENABLE_MULTIBLOCK:BOOL=ON
cd amr-wind-build-internal
make -j16
make install

cd ${TOP}
cmake -DCMAKE_INSTALL_PREFIX:PATH=${TOP}/erf-amrwind-driver/install \
      -DCMAKE_CXX_COMPILER:STRING=mpicxx \
      -DCMAKE_C_COMPILER:STRING=mpicc \
      -DCMAKE_Fortran_COMPILER:STRING=mpifort \
      -DMPIEXEC_PREFLAGS:STRING=--oversubscribe \
      -DCMAKE_BUILD_TYPE:STRING=Release \
      -DAMRWIND_HOME:STRING=${TOP}/amr-wind \
      -DERF_HOME:STRING=${TOP}/ERF \
      -DERF_DIM:STRING=3 \
      -DERF_ENABLE_MPI:BOOL=ON \
      -DERF_ENABLE_TESTS:BOOL=OFF \
      -DERF_ENABLE_FCOMPARE:BOOL=ON \
      -DERF_ENABLE_DOCUMENTATION:BOOL=OFF \
      -DERF_ENABLE_MULTIBLOCK:BOOL=ON \
      -DCMAKE_EXPORT_COMPILE_COMMANDS:BOOL=ON \
      -DCMAKE_PREFIX_PATH="${TOP}/amrex/install/lib/cmake/AMReX;${TOP}/AMReX-Hydro/install/lib/cmake/AMReX-Hydro;${TOP}/amr-wind-install/lib/cmake/AMR-Wind" \
      -DERF_USE_INTERNAL_AMREX:BOOL=OFF \
      -DAMR_WIND_USE_INTERNAL_AMREX:BOOL=OFF \
      -DAMR_WIND_USE_INTERNAL_AMREX_HYDRO:BOOL=OFF \
      -DDRIVER_USE_INTERNAL_AMRWIND:BOOL=ON \
      -DAMR_WIND_USE_INTERNAL_AMREX=OFF -DAMR_WIND_USE_INTERNAL_AMREX_HYDRO=OFF \
      -S ${TOP}/erf-amrwind-driver -B ${TOP}/erf-amrwind-driver/build
cd erf-amrwind-driver/build
make -j16
make install

cd ${TOP}

cd ${TOP}
cmake -DCMAKE_INSTALL_PREFIX:PATH=${TOP}/erf-amrwind-driver/install-internal \
      -DCMAKE_CXX_COMPILER:STRING=mpicxx \
      -DCMAKE_C_COMPILER:STRING=mpicc \
      -DCMAKE_Fortran_COMPILER:STRING=mpifort \
      -DMPIEXEC_PREFLAGS:STRING=--oversubscribe \
      -DCMAKE_BUILD_TYPE:STRING=Release \
      -DAMRWIND_HOME:STRING=${TOP}/amr-wind \
      -DERF_HOME:STRING=${TOP}/ERF \
      -DERF_DIM:STRING=3 \
      -DERF_ENABLE_MPI:BOOL=ON \
      -DERF_ENABLE_TESTS:BOOL=OFF \
      -DERF_ENABLE_FCOMPARE:BOOL=ON \
      -DERF_ENABLE_DOCUMENTATION:BOOL=OFF \
      -DERF_ENABLE_MULTIBLOCK:BOOL=ON \
      -DCMAKE_EXPORT_COMPILE_COMMANDS:BOOL=ON \
      -DCMAKE_PREFIX_PATH="${TOP}/amrex/install/lib/cmake/AMReX;${TOP}/AMReX-Hydro/install/lib/cmake/AMReX-Hydro;${TOP}/amr-wind-install/lib/cmake/AMR-Wind" \
      -DDRIVER_USE_INTERNAL_AMRWIND:BOOL=OFF \
      -DDRIVER_USE_INTERNAL_ERF:BOOL=ON \
      -DERF_USE_INTERNAL_AMREX:BOOL=OFF \
      -DAMR_WIND_USE_INTERNAL_AMREX=OFF -DAMR_WIND_USE_INTERNAL_AMREX_HYDRO=OFF \
      -S ${TOP}/erf-amrwind-driver -B ${TOP}/erf-amrwind-driver/build-internal
cd erf-amrwind-driver/build-internal
make -j16
make install

cd ${TOP}
