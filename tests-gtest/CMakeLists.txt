

# Include GTest Headers
include_directories(${PROJECT_SOURCE_DIR}/external)

if( ${testing})
  if ( ${QDPXX_CXX_FOUND} )
     list(APPEND HEADERS compare_spinors.h tolerance.h veclen.h)
     add_executable(tests_g main.cc clover_product.cc)

    #testlip is defined in ../tests
     target_link_libraries(tests_g testlib qphix_solver gtest_all ${EXT_LIBS})

     add_test(NAME gtest-runner COMMAND mpirun -n 2 ./tests_g -x ${test_runner_args})

 endif()
endif()
