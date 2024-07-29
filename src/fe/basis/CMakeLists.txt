#==============================================================================
target_sources(sfem PRIVATE
${CMAKE_CURRENT_SOURCE_DIR}/basis.cc
${CMAKE_CURRENT_SOURCE_DIR}/point_basis.cc
${CMAKE_CURRENT_SOURCE_DIR}/line_basis.cc
${CMAKE_CURRENT_SOURCE_DIR}/triangle_basis.cc
${CMAKE_CURRENT_SOURCE_DIR}/quad_basis.cc
${CMAKE_CURRENT_SOURCE_DIR}/tetrahedral_basis.cc
${CMAKE_CURRENT_SOURCE_DIR}/hexahedral_basis.cc)