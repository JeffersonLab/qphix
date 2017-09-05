#
# Attempt at a dumb config file for QMP
# 
#

# Find the qmp-config program to populate teh various
# variables
find_program(QMP_CONFIG NAMES qmp-config)

# If not found Error and exit
if( NOT QMP_CONFIG ) 
  message(FATAL_ERROR "qmp-config not found")  
endif( NOT QMP_CONFIG)


message(STATUS "Found qmp-config: ${QMP_CONFIG}")

# Get Compiler
execute_process(COMMAND ${QMP_CONFIG} --cc 
			OUTPUT_VARIABLE QMP_C_COMPILER OUTPUT_STRIP_TRAILING_WHITESPACE )
message(STATUS "QMP C Compiler is: ${QMP_C_COMPILER}")

# Get CFLAGS
execute_process(COMMAND ${QMP_CONFIG} --copts 
			OUTPUT_VARIABLE QMP_C_COMPILE_FLAGS OUTPUT_STRIP_TRAILING_WHITESPACE )
message(STATUS "QMP CFLAGS are: ${QMP_C_COMPILE_FLAGS}")

# Get LDFLAGS
execute_process(COMMAND ${QMP_CONFIG} --ldflags 
			OUTPUT_VARIABLE QMP_C_LINK_FLAGS OUTPUT_STRIP_TRAILING_WHITESPACE )
message(STATUS "QMP LDFLAGS are: ${QMP_C_LINK_FLAGS}")

# Get LIBS
execute_process(COMMAND ${QMP_CONFIG} --libs
			OUTPUT_VARIABLE QMP_C_LIBRARIES OUTPUT_STRIP_TRAILING_WHITESPACE )
message(STATUS "QMP LIBS are: ${QMP_C_LIBRARIES}")
 
#Turn the various strings into lists for appending
separate_arguments(QMP_C_COMPILE_FLAGS)
separate_arguments(QMP_C_LINK_FLAGS)
separate_arguments(QMP_C_LIBRARIES)

# Find the include path by looking for the path containing qmp.h
find_path(QMP_C_INCLUDE_PATH qmp.h)
message(STATUS "QMP Include path is: ${QMP_C_INCLUDE_PATH}")

# OK, we are done
set ( QMP_C_FOUND TRUE )
