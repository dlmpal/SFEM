#==============================================================================
target_sources(sfem PRIVATE 
${CMAKE_CURRENT_SOURCE_DIR}/petsc_vec.cc
${CMAKE_CURRENT_SOURCE_DIR}/petsc_mat.cc
${CMAKE_CURRENT_SOURCE_DIR}/petsc_ksp.cc)