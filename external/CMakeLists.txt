add_library(gtest_all gtest/gtest-all.cc)
include_directories(${PROJECT_SOURCE_DIR}/external)
target_compile_options(gtest_all PUBLIC ${EXT_CXXFLAGS})
