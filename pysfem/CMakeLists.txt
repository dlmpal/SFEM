#==============================================================================
if (CMAKE_VERSION VERSION_LESS 3.18)
  set(DEV_MODULE Development)
else()
  set(DEV_MODULE Development.Module)
endif()
#==============================================================================
find_package(Python 3.8 COMPONENTS Interpreter ${DEV_MODULE} REQUIRED)
#==============================================================================
# Detect the installed nanobind package and import it into CMake
set(tsl-robin-map_DIR ~/Projects/nanobind/robin-map)
execute_process(
  COMMAND "${Python_EXECUTABLE}" -m nanobind --cmake_dir
  OUTPUT_STRIP_TRAILING_WHITESPACE OUTPUT_VARIABLE NB_DIR)
list(APPEND CMAKE_PREFIX_PATH "${NB_DIR}")
find_package(nanobind CONFIG REQUIRED)
#==============================================================================
# Create the module 
nanobind_add_module(pysfem 
    wrappers/common.cc
    wrappers/geo.cc
    wrappers/mesh.cc
    wrappers/la.cc
    wrappers/io.cc
    wrappers/fe.cc
    wrappers/sfem.cc)
#==============================================================================
# Link with the C++ library
target_link_libraries(pysfem PRIVATE sfem)
#==============================================================================
# Installation
include(GNUInstallDirs)
install( 
  TARGETS pysfem
  EXPORT sfemTargets
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}/pysfem)