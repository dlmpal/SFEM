#==============================================================================
target_sources(sfem PRIVATE 
${CMAKE_CURRENT_SOURCE_DIR}/heat_conduction.cc
${CMAKE_CURRENT_SOURCE_DIR}/heat_convection.cc
${CMAKE_CURRENT_SOURCE_DIR}/heat_flux.cc)