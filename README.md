# ERF / AMR-Wind Driver

This repository contains the C++ API for coupling the [ERF](https://github.com/erf-model/ERF) and [AMR-Wind](https://github.com/erf-model/amr-wind) codes.

## Compile instructions

We recommend cloning this repository under a separate parent directory as the build and install directories for all dependencies and the driver itself will be created at the same level as the repository.
The `megabuild.sh` script updates the submodule dependencies, compiles each of them, and then compiles the driver code.
```
mkdir driver-build
cd driver-build
git clone git@github.com:erf-model/erf-amrwind-driver.git
./erf-amrwind-driver/Build/megabuild.sh
```

For testing any of the programs, for instance, the CouetteFlow problem:
```
cd erf-amrwind-driver-install/bin
cp ../../erf-amrwind-driver/Exec/CouetteFlow/input* .
./erf_mb_couette inputs_amrex
```

[![Linux GCC Build and Run](https://github.com/mukul1992/erf-amrwind-driver/actions/workflows/gcc.yml/badge.svg)](https://github.com/mukul1992/erf-amrwind-driver/actions/workflows/gcc.yml)
