#==============================================================================
target_sources(sfem PRIVATE 
${CMAKE_CURRENT_SOURCE_DIR}/dense_matrix.cc
${CMAKE_CURRENT_SOURCE_DIR}/sparsity_pattern.cc)
#==============================================================================
add_subdirectory(petsc)
#==============================================================================
add_subdirectory(slepc)