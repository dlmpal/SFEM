#!/bin/bash

mkdir -p fields
rm fields/*

cmake -S . -B build
(cd build; make)

mpiexec -np 2 build/elasticitySolver2D mesh/mesh_quad2 -ksp_type cg 1e-8
(cd fields; ${SFEM_DIR}/bin/sfemToVTK ../mesh/mesh_quad2 1 U stress)