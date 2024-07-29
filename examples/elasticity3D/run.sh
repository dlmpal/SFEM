#!/bin/bash

mkdir -p fields
rm fields/*

cmake -S . -B build
(cd build; make)

mpiexec -np 1 build/elasticitySolver3D mesh/mesh_tet2 -ksp_type cg -ksp_monitor
(cd fields; ${SFEM_DIR}/bin/sfemToVTK ../mesh/mesh_tet2 1 U stress) 