# ERF / AMR-Wind Driver

This repository contains the C++ API for coupling the [ERF](https://github.com/erf-model/ERF) and [AMR-Wind](https://github.com/erf-model/amr-wind) codes.

## Compile instructions

Create a directory for the binary tree, configure using CMake, and compile:
```
mkdir mybuild
cmake <options> ..
make -j8
```

The ERF and AMR-Wind root directory paths need to provided as the `ERF_HOME` and `AMRWIND_HOME` CMake variables. Template CMake scripts are provided in the `Build` directory for reference.

[![Linux GCC Build and Run](https://github.com/mukul1992/erf-amrwind-driver/actions/workflows/gcc.yml/badge.svg)](https://github.com/mukul1992/erf-amrwind-driver/actions/workflows/gcc.yml)
