
add_executable(time_dslash_noqdp )
target_link_libraries(time_dslash_noqdp    )

if( QDPXX_CXX_FOUND ) 
  add_executable( t_dslash  )
  target_link_libraries( t_dslash )

  add_executable( t_minvcg  )
  target_link_libraries( t_minvcg )

  add_executable( test_blas )
  target_link_libraries( test_blas   )
endif()

if( ${clover} ) 
  add_executable(time_clov_noqdp  )
  target_link_libraries( time_clov_noqdp ) 

  if( QDPXX_CXX_FOUND ) 
     add_executable( t_clov_dslash  )
     target_link_libraries( t_clov_dslash  
