#==============================================================================
# Create the library
add_library(sfem SHARED)

set(SFEM_DIRS 
        common
        mesh
        geo
        io
        la
        fe)

foreach(DIR ${SFEM_DIRS})
        add_subdirectory(${DIR})
endforeach()

target_include_directories(sfem 
                           PUBLIC 
                           "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>"
                           "$<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>")
#==============================================================================
# Third-party dependencies

## PETSc
set(SFEM_HAS_PETSC TRUE)
if(DEFINED PETSC_DIR AND DEFINED PETSC_ARCH)
    target_include_directories(sfem PUBLIC ${PETSC_DIR}/include ${PETSC_DIR}/${PETSC_ARCH}/include)
    target_link_directories(sfem PUBLIC ${PETSC_DIR}/${PETSC_ARCH}/lib)
    target_link_libraries(sfem PUBLIC petsc)
    target_compile_definitions(sfem PUBLIC SFEM_HAS_PETSC)
else()
    message(FATAL_ERROR "SFEM currently only supports PETSc as the linear algebra backend!")
    set(SFEM_HAS_PETSC FALSE)
endif()

## SLEPc
if(DEFINED SLEPC_DIR)
    if(${SFEM_HAS_PETSC} EQUAL FALSE)
        message(FATAL_ERROR "SLEPc requires PETSc backend!")
    endif()
    target_include_directories(sfem PUBLIC ${SLEPC_DIR}/include ${SLEPC_DIR}/${PETSC_ARCH}/include)
    target_link_directories(sfem PUBLIC ${SLEPC_DIR}/${PETSC_ARCH}/lib)
    target_link_libraries(sfem PUBLIC slepc)
    target_compile_definitions(sfem PUBLIC SFEM_HAS_SLEPC)
endif()

## METIS
if(DEFINED METIS_DIR)
    target_include_directories(sfem PUBLIC ${METIS_DIR}/include)
    target_link_directories(sfem PUBLIC ${METIS_DIR}/lib)
    target_link_libraries(sfem PUBLIC metis GKlib)
    target_compile_definitions(sfem PUBLIC SFEM_HAS_METIS)
else()
    message(FATAL_ERROR "SFEM currently only supports METIS for mesh partitioning!")
endif()

## MPI
if(DEFINED MPI_DIR)
    target_include_directories(sfem PUBLIC ${MPI_DIR}/include)
    target_link_directories(sfem PUBLIC ${MPI_DIR}/lib)
    target_link_libraries(sfem PUBLIC mpi)
    target_compile_definitions(sfem PUBLIC SFEM_USE_MPI)
endif()
#==============================================================================
# Installation
include(GNUInstallDirs)
install( 
  TARGETS sfem
  EXPORT sfemTargets
  INCLUDES DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR})

## Install headers
install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/
        DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
        FILES_MATCHING PATTERN "*.h")