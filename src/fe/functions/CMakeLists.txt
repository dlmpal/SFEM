#==============================================================================
target_sources(sfem PRIVATE
${CMAKE_CURRENT_SOURCE_DIR}/field_gradient.cc
${CMAKE_CURRENT_SOURCE_DIR}/structural_mass.cc
${CMAKE_CURRENT_SOURCE_DIR}/stress.cc)