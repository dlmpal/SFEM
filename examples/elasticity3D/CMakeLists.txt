cmake_minimum_required(VERSION 3.16)
project(elasticity_solver_3d)

set(sfem_DIR $ENV{SFEM_DIR}/lib/cmake/sfem)
find_package(sfem REQUIRED)

add_executable(elasticitySolver3D ../main.cc)
target_link_libraries(elasticitySolver3D sfem::sfem)