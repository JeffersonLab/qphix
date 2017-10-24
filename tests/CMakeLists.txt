add_compile_options(${EXT_CXXFLAGS})

list(APPEND framework_sources cli_args.cc)

if( ${testing} )
  if( ${QDPXX_CXX_FOUND} )
     list(APPEND framework_headers unittest.h dslashm_w.h reunit.h invbicgstab_test.h mesfield.h)

     list(APPEND framework_sources dslashm_w.cc
                                   reunit.cc
                                   mesfield.cc
                                   tolerance.cpp)

     if( ${clover} OR ${tm_clover} )
       list(APPEND framework_headers clover_fermact_params_w.h )
   
       if( ${qdpjit} )
        list(APPEND framework_headers clover_term_llvm_w.h)
       else()
        list(APPEND framework_headers clover_term_qdp_w.h)
       endif()
     endif()
    
     add_library(testlib ${framework_sources} ${framework_headers} ${header_files})
    target_link_libraries(testlib qphix_codegen)
    
    add_executable(time_dslash_noqdp time_dslash_noqdp.cc timeDslashNoQDP.h timeDslashNoQDP.cc)
    target_link_libraries(time_dslash_noqdp testlib qphix_solver ${EXT_LIBS})
    list( APPEND built_targets time_dslash_noqdp)
    
    if( ${clover} )
      add_executable(time_clov_noqdp time_clov_noqdp.cc timeClovNoQDP.h timeClovNoQDP.cc)
      target_link_libraries(time_clov_noqdp testlib qphix_codegen qphix_solver ${EXT_LIBS})
      list( APPEND built_targets time_clov_noqdp)
      add_test(NAME clover COMMAND mpirun -n 2 ./t_clov_dslash ${test_runner_args})
    endif()
    
    if( ${twisted_mass} )
     add_executable(time_twm_noqdp time_twm_noqdp.cc timeTWMNoQDP.h timeTWMNoQDP.cc)
     target_link_libraries(time_twm_noqdp testlib qphix_codegen qphix_solver ${EXT_LIBS})
     list( APPEND built_targets time_twm_noqdp)
      add_test(NAME tm COMMAND mpirun -n 2 ./t_twm_dslash ${test_runner_args})
    endif()
    
    if( ${tm_clover} )
      add_executable(time_twm_clover time_twm_clover.cc timeTWMClover.h timeTWMClover.cc)
      target_link_libraries(time_twm_clover testlib qphix_codegen qphix_solver ${EXT_LIBS})
      add_test(NAME tm_clover COMMAND mpirun -n 2 ./t_twm_clover ${test_runner_args})
    endif()
 
    add_executable(t_dslash t_dslash.cc testDslashFull.h testDslashFull.cc ${framework_headers})
    target_link_libraries(t_dslash testlib qphix_codegen qphix_solver ${EXT_LIBS})
    list(APPEND built_targets t_dslash)
      add_test(NAME wilson COMMAND mpirun -n 2 ./t_dslash ${test_runner_args})

    add_executable(test_blas test_blas.cc testBlas.h testBlas.cc ${framework_headers})
    target_link_libraries(test_blas testlib qphix_codegen qphix_solver ${EXT_LIBS})
    list(APPEND built_targets test_blas)

    if( ${clover} )
	add_executable(t_clov_dslash  t_clov_dslash.cc testClovDslashFull.h testClovDslashFull.cc ${framework_headers})
	target_link_libraries(t_clov_dslash testlib qphix_codegen qphix_solver ${EXT_LIBS})

	add_executable(t_clov_invert_from_file t_clov_invert_from_file.cc
					       testClovInvertFromFile.h 
					       testClovInvertFromFile.cc
						${framework_headers})
	target_link_libraries(t_clov_invert_from_file testlib qphix_codegen qphix_solver ${EXT_LIBS})
        list(APPEND built_targets t_clov_dslash t_clov_invert_from_file)
    endif()

    if( ${twisted_mass} )
	add_executable(t_twm_dslash t_twm_dslash.cc testTWMDslashFull.h testTWMDslashFull.cc ${framework_headers})
	target_link_libraries(t_twm_dslash testlib qphix_codegen qphix_solver ${EXT_LIBS})
	list(APPEND built_targets t_twm_dslash)
    endif()

    if( ${tm_clover}) 
        add_executable(t_twm_clover t_twm_clover.cc testTWMCloverFull.h testTWMCloverFull.cc ${framework_headers})
	target_link_libraries(t_twm_clover testlib qphix_codegen qphix_solver ${EXT_LIBS})
	list(APPEND buil_targets t_twm_clover)
    endif()
  endif()
endif()

INSTALL(TARGETS ${built_targets} DESTINATION bin)

