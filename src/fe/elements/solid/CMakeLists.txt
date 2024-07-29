#==============================================================================
target_sources(sfem PRIVATE
${CMAKE_CURRENT_SOURCE_DIR}/linear_elasticity.cc
${CMAKE_CURRENT_SOURCE_DIR}/pressure_load.cc
${CMAKE_CURRENT_SOURCE_DIR}/thermal_load.cc
${CMAKE_CURRENT_SOURCE_DIR}/inertial_load.cc)