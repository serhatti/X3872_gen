# author : Serhat Istin

include(FetchContent)
FetchContent_Declare(fmt
  GIT_REPOSITORY https://github.com/fmtlib/fmt.git
  GIT_TAG master
)
FetchContent_MakeAvailable(fmt)

find_package(TBB REQUIRED)
include_directories(${PYTHIA8_INCDIR} ${ROOT_INCDIR})

#any source file here with .cc extension is considered to be a main program
FILE(GLOB CPPFiles *.cc *.cpp)

FILE(GLOB UserHeaders *.h)

#compile and link all the mains

foreach(EXECSRC ${CPPFiles})
    get_filename_component(appname ${EXECSRC} NAME_WE  )
    add_executable(${appname} ${EXECSRC} )
    target_link_libraries( ${appname} ${PYTHIA8_LIBRARIES} ${ROOT_LIBRARIES} tbb fmt::fmt)
    install (TARGETS ${appname} DESTINATION ${BINDIR} )
endforeach()