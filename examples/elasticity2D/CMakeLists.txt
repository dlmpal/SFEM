cmake_minimum_required(VERSION 3.16)
project(elasticity_solver_2d)

set(CMAKE_CXX_STANDARD 17)
set(sfem_DIR $ENV{SFEM_DIR}/lib/cmake/sfem)
find_package(sfem REQUIRED)

add_executable(elasticitySolver2D ../main.cc)
target_link_libraries(elasticitySolver2D sfem::sfem)