#author : Serhat Istin

# We can compile and link pythia examples with ROOT by this CMakeLists 
#top-level CMakeLists file 


cmake_minimum_required (VERSION  3.5)
project(pythia8X3872)



set(INSTALLATION ${PROJECT_BINARY_DIR}/INSTALL)
set(BINDIR ${INSTALLATION}/bin)


find_program( ROOT-CONFIG root-config) 
if( NOT ROOT-CONFIG )
  message( FATAL_ERROR "root-config not found !" )
endif()


execute_process(
  COMMAND ${ROOT-CONFIG} --incdir
  OUTPUT_VARIABLE ROOT_INCDIR
  OUTPUT_STRIP_TRAILING_WHITESPACE
)

execute_process(
  COMMAND ${ROOT-CONFIG} --cflags
  OUTPUT_VARIABLE ROOT_CXX_FLAGS
  OUTPUT_STRIP_TRAILING_WHITESPACE
)

execute_process(
  COMMAND ${ROOT-CONFIG} --libs
  OUTPUT_VARIABLE ROOT_LIBRARIES
  OUTPUT_STRIP_TRAILING_WHITESPACE
)


#Pick up PYTHIA8 from the environment and link against the precompiled library

find_program( PYTHIA8-CONFIG pythia8-config) 
if( NOT PYTHIA8-CONFIG )
  message( FATAL_ERROR "pythia8-config not found !" )
endif()


execute_process(
  COMMAND ${PYTHIA8-CONFIG} --includedir
  OUTPUT_VARIABLE PYTHIA8_INCDIR
  OUTPUT_STRIP_TRAILING_WHITESPACE
)

execute_process(
  COMMAND ${PYTHIA8-CONFIG} --cxxflags
  OUTPUT_VARIABLE PYTHIA8_CXX_FLAGS
  OUTPUT_STRIP_TRAILING_WHITESPACE
)

execute_process(
  COMMAND ${PYTHIA8-CONFIG} --libs
  OUTPUT_VARIABLE PYTHIA8_LIBRARIES
  OUTPUT_STRIP_TRAILING_WHITESPACE
)

message(STATUS "ROOT includes: ${ROOT_INCDIR}")
message(STATUS "ROOT CXX flags: ${ROOT_CXX_FLAGS}")
message(STATUS "ROOT libraries: ${ROOT_LIBRARIES}")


message(STATUS "Pythia8 includes: ${PYTHIA8_INCDIR}")
message(STATUS "Pythia8 CXX flags: ${PYTHIA8_CXX_FLAGS}")
message(STATUS "Pythia8 libraries: ${PYTHIA8_LIBRARIES}")


add_subdirectory(generator)


