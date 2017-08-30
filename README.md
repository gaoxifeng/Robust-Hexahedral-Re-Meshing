Robust Hexahedral Re-Meshing
===================

This repository contains the meshing software developed as part of the publication

Robust Structure Simplification for Hex Re-meshing

Xifeng Gao, Daniele Panozzo, Wenping Wang, Zhigang Deng, Guoning Chen

In ACM Transactions on Graphics (Proceedings of SIGGRAPH ASIA 2017)

Compiling
-------------
Compiling from scratch requires CMake and Visual Studio 2015 on Windows 10.

git clone --recursive https://github.com/gaoxifeng/Robust-Hexahedral-Re-Meshing.git

On Windows, open the generated file complex_simplification.sln after CMake compilation and proceed building as usual from within Visual Studio.

Usage
-------------
**complex_simplification.exe c r b s f i**

**c**--processing type. SIM represents simplification of the structure in the input mesh, OPT represents purely geometric optimization of the input mesh,

**r**-- a ratio of the resolution of the output divided by the resolution of the input, which controls the hex-element resolution of the output,

**b**--a base number for the optimization iteration, default value is 2,

**s**--the ratio of the to-be-removed blocks of the structure in the output divided by the number of blocks of the structure in the input,

**f**--denote that the input mesh contains sharp feature or not. 1 means yes, 0 means no,

**i**--the input (.vtk format only).

**An example command for simplification**: 
complex_simplification_SIM.exe SIM 1 2 1 0 ../../Db_data_movies/Octree/airplane1_input_tri_hexa

**An example command for optimization**: 
complex_simplification_SIM.exe OPT 1 2 1 0 ../../Db_data_movies/Octree/airplane1_input_tri_hexa
