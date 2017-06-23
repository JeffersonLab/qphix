
list(APPEND source_files print_utils.cc
 	memory_usage.cc 
	qphix_cli_args.cc)

if( ${isa} STREQUAL "qpx" )
  list( APPEND source_files bgq_threadbind.cc)
else()
  list( APPEND source_files generic_threadbind.cc)
endif()
 
add_library(qphix_solver STATIC ${source_files})
target_link_libraries(qphix_solver qphix_codegen)

target_compile_options(qphix_solver PUBLIC ${EXT_CXXFLAGS})
install(TARGETS qphix_solver EXPORT QPhiXSolverConfig
	 ARCHIVE  DESTINATION lib
         LIBRARY  DESTINATION lib
         RUNTIME  DESTINATION bin)

install(EXPORT QPhiXSolverConfig DESTINATION cmake)

