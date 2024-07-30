#!/bin/bash

mkdir -p fields
rm fields/*

cmake -S . -B build
(cd build; make)

mpiexec -np 4 build/laplaceSolver2D mesh -ksp_type cg -ksp_monitor
(cd fields; ${SFEM_DIR}/bin/sfemToVTK ../mesh 1 phi U)