#==============================================================================
target_sources(sfem PRIVATE
${CMAKE_CURRENT_SOURCE_DIR}/logger.cc
${CMAKE_CURRENT_SOURCE_DIR}/error.cc
${CMAKE_CURRENT_SOURCE_DIR}/timer.cc
${CMAKE_CURRENT_SOURCE_DIR}/init.cc
${CMAKE_CURRENT_SOURCE_DIR}/index_map.cc
${CMAKE_CURRENT_SOURCE_DIR}/mpi_utils.cc)