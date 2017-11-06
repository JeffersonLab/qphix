

# Include GTest Headers
include_directories(${PROJECT_SOURCE_DIR}/external)

include_directories(${PROJECT_SOURCE_DIR})

if( ${testing})
  if ( ${QDPXX_CXX_FOUND} )
     add_library(gtestlib STATIC process_args.cc random_gauge.cc qdpxx_test_env.cc)
     add_dependencies(gtestlib qphix_solver qphix_codegen)
     target_compile_options(gtestlib PUBLIC ${EXT_CXXFLAGS})

     list(APPEND HEADERS compare_spinors.h tolerance.h veclen.h process_args.h qdpxx_test_env.h)
     if( !${qdpjit} ) 
       if( ${tm_clover} )
         add_executable(tests_g main.cc clover_product.cc ndtm_herm.cc a_imugamma5.cc richardson_cg.cc)	 
         target_link_libraries(tests_g testlib qphix_codegen qphix_solver gtest_all gtestlib ${EXT_LIBS})
         add_test(NAME gtest-runner COMMAND mpirun -n 1 ./tests_g ${test_runner_args} -x 8 -y 4 -z 4 -t 4 -geom 1 1 1 1)
       endif()
     endif()

     add_executable(solver_test solver_test.cc)
     target_link_libraries(solver_test testlib qphix_solver gtest_all gtestlib ${EXT_LIBS})
     add_test(NAME solver-test COMMAND mpirun -n 1 ./solver_test ${test_runner_args} -x 8 -y 4 -z 4 -t 4 -geom 1 1 1 1)

     add_executable(dslash_dir_test dslash_dir_test.cc)
     target_link_libraries(dslash_dir_test testlib qphix_solver gtest_all gtestlib ${EXT_LIBS})
     add_test(NAME dslash-dir-test COMMAND mpirun -n 1 ./dslash_dir_test ${test_runner_args} -x 8 -y 4 -z 4 -t 4 -geom 1 1 1 1)
     add_test(NAME dslash-dir-test-T COMMAND mpirun -n 2 ./dslash_dir_test ${test_runner_args} -x 8 -y 4 -z 4 -t 4 -geom 1 1 1 2)
     add_test(NAME dslash-dir-test-Z COMMAND mpirun -n 2 ./dslash_dir_test ${test_runner_args} -x 8 -y 4 -z 8 -t 4 -geom 1 1 2 1)
     add_test(NAME dslash-dir-test-Y COMMAND mpirun -n 2 ./dslash_dir_test ${test_runner_args} -x 8 -y 8 -z 4 -t 4 -geom 1 2 1 1)
     if(x_comms_testing)
       add_test(NAME dslash-dir-test-X COMMAND mpirun -n 2 ./dslash_dir_test ${test_runner_args} -x 16 -y 4 -z 4 -t 4 -geom 2 1 1 1)
     endif()
         
     add_executable(full_blas_test full_blas_test.cc)
     target_link_libraries(full_blas_test testlib qphix_solver gtest_all gtestlib ${EXT_LIBS})
     add_test(NAME full-blas-test COMMAND mpirun -n 1 ./full_blas_test ${test_runner_args} -x 8 -y 4 -z 4 -t 4 -geom 1 1 1 1)

 endif()
endif()
