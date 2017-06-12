add_compile_options(${EXT_CXXFLAGS})

add_executable(time_dslash_noqdp time_dslash_noqdp.cc timeDslashNoQDP.h timeDslashNoQDP.cc cli_args.cc )
target_link_libraries(time_dslash_noqdp qphix_codegen qphix_solver ${EXT_LIBS})

