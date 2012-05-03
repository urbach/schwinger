cmake_minimum_required(VERSION 2.4)

# The name of our project is "HELLO".  CMakeLists files in this project can
# refer to the root source directory of the project as ${HELLO_SOURCE_DIR} and
# to the root binary directory of the project as ${HELLO_BINARY_DIR}.
project (maxmize_tr)


# intel compiling
# SET ( CMAKE_C_COMPILER "icc" )

SET ( CMAKE_C_FLAGS "-std=c99 -pg -fopenmp -DOMP" )
## SET ( CMAKE_C_FLAGS "-std=c99 -pg" )

INCLUDE_DIRECTORIES(./)

add_subdirectory( rand )





SET (BASE_SRC rand/gauss.c  rand/ranlxd.c  rand/ranlxs.c 2MN_integrator.c  fields.c   leapfrog.c   measurements.c  statistics.c  hmc.c   leapfrog2.c  dirac.c  lattice.c  linalg.c  rec_lf_integrator.c )



ADD_EXECUTABLE(qed qed.c ${BASE_SRC})


set( LIBS ${LIBS} "-lm" )

target_link_libraries( qed ${LIBS} )




