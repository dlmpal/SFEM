cmake_minimum_required(VERSION 3.16)
project(laplace_solver_2d)

set(sfem_DIR $ENV{SFEM_DIR}/lib/cmake/sfem)
find_package(sfem REQUIRED)

add_executable(laplaceSolver2D ../main.cc)
target_link_libraries(laplaceSolver2D sfem::sfem)