

# Include GTest Headers
include_directories(${PROJECT_SOURCE_DIR}/external)

if( ${testing})
  if ( ${QDPXX_CXX_FOUND} )
     add_executable(tests_g
         main.cc
         ndtm_herm.cc
         clover_product.cc
         a_imugamma5.cc
         richardson_cg.cc
         )

     add_dependencies(tests_g qphix_codegen qphix_solver)

    #testlip is defined in ../tests
     target_link_libraries(tests_g testlib qphix_solver gtest_all ${EXT_LIBS})

     add_test(NAME gtest-runner COMMAND mpirun -n 1 ./tests_g ${test_runner_args} -x 8 -y 4 -z 4 -t 4 -geom 1 1 1 1)

 endif()
endif()
