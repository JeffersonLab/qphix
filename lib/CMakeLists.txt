
list(APPEND ${source_files} print_util.cc
 	memory_usage.cc 
	qphix_cli_args.cc)

if( ${isa} STREQUAL qpx )
  list( APPEND ${source_files} bgq_threadbind.cc)
else()
  list( APPEND ${source_files} generic_threadbind.cc)
endif()
 
add_library(qphix_solver STATIC $(source_files})

