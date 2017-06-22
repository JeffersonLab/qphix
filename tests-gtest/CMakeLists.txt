# Pick up all the yucky compile options if not yet done

# Include GTest Headers
include_directories(${PROJECT_SOURCE_DIR}/external)

if( ${testing})
  if ( ${QDPXX_CXX_FOUND} )
     list(APPEND HEADERS compare_spinors.h tolerance.h veclen.h)
     add_executable(tests_g main.cc clover_product.cc)

    #testlip is defined in ../tests
     target_link_libraries(tests_g testlib qphix_solver gtest_all ${EXT_LIBS})

     add_test(NAME TestCloverProduct COMMAND ./tests_g -x 8 -y 4 -z 4 -t 4 -by 4 -bz 4 -pxy 0 -pxyz 0 -minct 1 -c 1 -sy 1 -sz 2 -geom 1 1 1 1)

 endif()
endif()
