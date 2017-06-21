#
# Attempt at a dumb config file for QDP++
# 
#

# Find the qmp-config program to populate teh various
# variables

find_program(QDPXX_CONFIG NAMES qdp++-config)

# If not found Error and exit
if( NOT QDPXX_CONFIG ) 
  message(FATAL_ERROR "qdpxx-config not found")  
endif( NOT QDPXX_CONFIG)


message(STATUS "Found qdp++-config: ${QDPXX_CONFIG}")

# Get Compiler
execute_process(COMMAND ${QDPXX_CONFIG} --cxx 
			OUTPUT_VARIABLE QDPXX_CXX_COMPILER OUTPUT_STRIP_TRAILING_WHITESPACE )
message(STATUS "QDP++ C++ Compiler is: ${QDPXX_CXX_COMPILER}")

# Get CFLAGS
execute_process(COMMAND ${QDPXX_CONFIG} --cxxflags
			OUTPUT_VARIABLE QDPXX_CXX_COMPILE_FLAGS OUTPUT_STRIP_TRAILING_WHITESPACE )
message(STATUS "QDP++ CXXFLAGS are: ${QDPXX_CXX_COMPILE_FLAGS}")

# Get LDFLAGS
execute_process(COMMAND ${QDPXX_CONFIG} --ldflags 
			OUTPUT_VARIABLE QDPXX_CXX_LINK_FLAGS OUTPUT_STRIP_TRAILING_WHITESPACE )
message(STATUS "QDPXX LDFLAGS are: ${QDPXX_CXX_LINK_FLAGS}")

# Get LIBS
execute_process(COMMAND ${QDPXX_CONFIG} --libs
			OUTPUT_VARIABLE QDPXX_CXX_LIBRARIES OUTPUT_STRIP_TRAILING_WHITESPACE )
message(STATUS "QDP++ LIBS are: ${QDPXX_CXX_LIBRARIES}")


# Get Parallel Arch
execute_process(COMMAND ${QDPXX_CONFIG} --parallel-arch
	OUTPUT_VARIABLE QDPXX_PARALLEL_ARCH OUTPUT_STRIP_TRAILING_WHITESPACE )
message(STATUS "QDP++ Parallel Architecture is ${QDPXX_PARALLEL_ARCH}") 
 
#Turn the various strings into lists for appending
separate_arguments(QDPXX_CXX_COMPILE_FLAGS)
separate_arguments(QDPXX_CXX_LINK_FLAGS)
separate_arguments(QDPXX_CXX_LIBRARIES)

# Find the include path by looking for the path containing qmp.h
find_path(QDPXX_CXX_INCLUDE_PATH qdp.h)
message(STATUS "QDP++ Include path is: ${QDPXX_CXX_INCLUDE_PATH}")

# OK, we are done
set ( QDPXX_CXX_FOUND TRUE )
