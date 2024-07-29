#=============================================================================
# gmshToSfem
add_executable(gmshToSfem gmsh_to_sfem.cc)
target_link_libraries(gmshToSfem PRIVATE sfem)
#==============================================================================
# sfemToVTK
add_executable(sfemToVTK sfem_to_vtk.cc)
target_link_libraries(sfemToVTK PRIVATE sfem)
#==============================================================================
# Installation
include(GNUInstallDirs)
install( 
  TARGETS gmshToSfem sfemToVTK
  EXPORT sfemTargets
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})