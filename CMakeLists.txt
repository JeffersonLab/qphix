cmake_minimum_required(VERSION 3.1)
project(QPhiX)
include(CTest)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_EXTENSIONS OFF)

#add_definitions(-Wall -pedantic
#  -Wno-unused-variable
#  -Wno-sign-compare
#  )

###############################################################################
#                                User Options                                 #
###############################################################################

set(parallel_arch scalar CACHE STRING  "Build for parallel arch: scalar/parscalar")
option(testing "Enable building of tests. This can take a long time" OFF)
option(cean "Enable CILK array Notation" OFF)
option(mm_malloc "Use mm_malloc to allocate memory (unless qdpalloc is enabled, in which case QDP++ allocator will be used" ON)
option(clover "Enable building of clover term. Default is ON" ON)
option(twisted_mass "Enable building of twisted mass term. Default is OFF" OFF)
option(tm_clover "Enable building of twisted mass + clover term. Default if OFF" OFF)
option(fake_comms "Enable fake comms. Does face processing but no actual comms" OFF)
option(extra_messages "Emit extram messages during compilation time" OFF)
set(isa avx CACHE STRING "Enable use of code for particular isa. Allowed are: avx avx2 avx512 sse scalar (defaults to avx)")
option(qdpjit "When QDP++ is used it is QDP-JIT so use QDP-JIT interfaces" OFF)
option(qdpalloc "When QDP++ is used, use its memory allocator instead of mm_malloc or posix_memalign" OFF)
option(mpi_comms "Use plain MPI comms instead of QMP in the kernels" ON)

# Codegen stuff
set(host_cxx ${CMAKE_CXX_COMPILER} CACHE STRING "select target CXX Compiler for building libqphix-codegen.a" )
set(host_cxxflags ${CMAKE_CXX_FLAGS} CACHE STRING  "select target CXXFLAGS for building libqphix-codege.a")
set(recusive_jN '' CACHE STRING "select -j value for recursive make (defaults to 1)")
option(skip_build "Developer only: Reduce building of certain kernels for faster build" OFF)

###############################################################################

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/CMakeModules/")
set(QDPXX_CXX_FOUND FALSE)
set(QMP_C_FOUND FALSE)

# Include our own files first.
include_directories(${PROJECT_SOURCE_DIR}/include)
include_directories(${PROJECT_BINARY_DIR}/include)
if ( ${testing} ) 
  include_directories(${PROJECT_SOURCE_DIR}/external )
endif()

###############################################################################
#                                   OpenMP                                    #
###############################################################################

find_package(OpenMP)
if(OPENMP_FOUND)
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
else()
  add_definitions("-Wno-unknown-pragmas")
  message(ERROR "OpenMP not found")
endif()

###############################################################################
#                                    QDP++                                    #
###############################################################################

if(QDPXX_DIR)
  message(STATUS "Looking for QDP++")
  set(CMAKE_PREFIX_PATH "${CMAKE_PREFIX_PATH}" "${QDPXX_DIR}")
  find_package(QDPXX)
  if(${QDPXX_CXX_FOUND})
    include_directories(${QDPXX_CXX_INCLUDE_PATH})
    list(APPEND EXT_CXXFLAGS "${QDPXX_CXX_COMPILE_FLAGS}")
    list(REMOVE_DUPLICATES EXT_CXXFLAGS)
    list(APPEND EXT_LIBS ${QDPXX_CXX_LINK_FLAGS})
    list(APPEND EXT_LIBS ${QDPXX_CXX_LIBRARIES})
    set(QPHIX_USE_QDPXX TRUE)
    set(parallel_arch ${QDPXX_PARALLEL_ARCH})
    if(${qdpalloc}) 
      message(INFO "Using QDPAlloc for memory allocation")
      set(QPHIX_USE_QDPXX_ALLOC 1)
    endif()
  endif()
else()
  message(STATUS "Not Looking for QDP++")
  set(QPHIX_USE_QDPXX FALSE)
endif()

###############################################################################
#                                     QMP                                     #
###############################################################################

if(QMP_DIR)
  message(STATUS "Looking for QMP")
  set(CMAKE_PREFIX_PATH "${CMAKE_PREFIX_PATH}" "${QMP_DIR}")

  # Find package QMP
  find_package(QMP)

  # Check for QMP
  if(${QMP_C_FOUND})
    include_directories(${QMP_C_INCLUDE_PATH})
    list(APPEND ${EXT_LIBS} ${QMP_C_LINK_FLAGS})
    list(APPEND ${EXT_LIBS} ${QMP_C_LIBRARIES})
  endif()
else()
  message(STATUS "Not Looking for QMP")
endif()

###############################################################################
#                            Parallel Architecture                            #
###############################################################################

if(${parallel_arch} STREQUAL "scalar")
  # set in configure file
  set(QPHIX_SCALAR 1)
  message(STATUS "Scalar Comms Mode (single node)")
  if(${fake_comms}) 
    message(STATUS "Fake Comms will be used") 
    set(QPHIX_DO_COMMS 1)
    set(QPHIX_FAKE_COMMS 1)
  endif()
elseif(${parallel_arch} STREQUAL "parscalar")
  message(STATUS "Parscalar comms mode (QMP + OpenMP)")
  set(QPHIX_PARSCALAR 1)

  if(${fake_comms})  
    message(FATAL_ERROR "Cannot do fake comms in parscalar mode")
  endif()

  if(NOT (${QDPXX_CXX_FOUND} OR ${QMP_C_FOUND}))
    message(FATAL_ERROR "Parscalar arch requested but neither QMP nor QDPXX specified")
  else()
    message(STATUS "Regular QMP Comms")
    set(QPHIX_DO_COMMS 1)
    set(QPHIX_QMP_COMMS 1)
  endif()
else()
  message(FATAL_ERROR "Only scalar or parscalar modes supported")
endif()

###############################################################################

if(${cean})
  message(STATUS "Setting CEAN")
  set(QPHIX_USE_CEAN 1)
endif()

###############################################################################

if (mpi_comms)
  set(QPHIX_MPI_COMMS_CALLS TRUE)
endif()

###############################################################################

if(NOT ${QPHIX_USE_QDPXX_ALLOC})
  if(${mm_malloc})
    message(STATUS "Setting use mm_malloc")
    set(QPHIX_USE_MM_MALLOC 1)
  else()
    message(STATUS "Posix Memalign allocator will be used")
  endif()
endif()

###############################################################################
#             Types of Action (Clover, Wilson, Twisted Mass, ...)             #
###############################################################################

if(${clover}) 
  message(STATUS "Enabling Clover")
  set(QPHIX_BUILD_CLOVER 1)
endif()

if(${twisted_mass}) 
  message(STATUS "Enabling Twisted-Mass")
  set(QPHIX_BUILD_TWISTED_MASS 1)
endif()

if(${tm_clover})
  message(STATUS "Enabling Twisted-Mass with Clover")
  set(QPHIX_BUILD_TWISTED_MASS_WITH_CLOVER 1)
endif()

###############################################################################
#                                     ISA                                     #
###############################################################################
set(x_comms_testing TRUE)
set(ISA ${isa})
set(QPHIX_LLC_CACHE_ALIGN 64)
if(${isa} STREQUAL "avx")
  message(STATUS "Building for AVX")
  set(QPHIX_AVX_SOURCE 1)
  set(test_soalen 4)
elseif(${isa} STREQUAL "avx2")
  message(STATUS "Building for AVX2")
  set(QPHIX_AVX2_SOURCE 1)
  set(test_soalen 4)
elseif(${isa} STREQUAL "avx512")
  message(STATUS "Building for AVX512")
  set(QPHIX_AVX512_SOURCE 1)
  set(test_soalen 8)
elseif(${isa} STREQUAL "sse")
  message(STATUS "Building for SSE")
  set(QPHIX_SSE_SOURCE 1)
elseif(${isa} STREQUAL "mic")
  message(STATUS "Building for MIC (KNC)")
  set(QPHIX_MIC_SOURCE 1)
elseif(${isa} STREQUAL "scalar")
  message(STATUS "Building for straight C (nonvectorized")
  set(QPHIX_SCALAR_SOURCE 1)
  set(test_soalen 1)
  set(x_comms_testing FALSE)
  message(WARNING "Disabling X-directional communicating tests")
  message(WARNING "This is due to a known issue (issue #10) for isa=scalar")
elseif(${isa} STREQUAL "qpx")
  message(STATUS "Building for QPX (BG/Q)")
  set(QPHIX_QPX_SOURCE 1)
  set(QPHIX_LLC_CACHE_ALIGN 128)
endif()

if( x_comms_testing )
   message(STATUS "Will do X-comms tests")
else()
   message(STATUS "Disabled X-comms tests")
endif()

###############################################################################

set(${build_qdpjit_clover} FALSE)
if(${qdpjit}) 
  message(STATUS "Building with QDPJIT")
  set(QPHIX_BUILD_QDPJIT 1)
  if(${clover})
    message(STATUS "Enabling LLVM Clover Term")
    set(${build_qdpjit_clover} TRUE)
  endif()
endif()

###############################################################################

if(${testing}) 
  message(STATUS "Building Tests")
else()
  message(STATUS "Not Building Tests")
endif()

###############################################################################

if(${extra_messages}) 
  message(STATUS "Setting Extra Messages")
  set(QPHIX_EMIT_MESSAGES 1)
endif()

###############################################################################
#                            QPhiX Code Generator                             #
###############################################################################

if(QPHIX_CODEGEN) 
  message(STATUS "Importing QPhix Codegen from ${QPHIX_CODEGEN}")
  include_directories(${QPHIX_CODEGEN}/include)
  add_library(qphix_codegen STATIC IMPORTED)
  set_property(TARGET qphix_codegen PROPERTY IMPORTED_LOCATION
    ${QPHIX_CODEGEN}/lib/libqphix_codegen.a)
else()

  MESSAGE(INFO, "Looking for Python")
  find_package(PythonInterp 3.0 REQUIRED)
  MESSAGE(INFO, "Python interpreter is ${PYTHON_EXECUTABLE}")
  # A build directory for building qphix-codegen
  file(MAKE_DIRECTORY ${CMAKE_BINARY_DIR}/codegen)

  if ("${recursive_jN}" STREQUAL "")
    set(make_program_arguments "")
  else()
    set(make_program_arguments "-j${recursive_jN}")
  endif()
  message(STATUS "make_program_arguments: ${make_program_arguments}")

  # The custom command will make at least this `libqphix_codegen.a`. It will
  # use the build directory (specified by the `WORKING_DIRECTORY`) When this
  # command is executed, it will also install the library to
  # `CMAKE_INSTALL_PREFIX`.
  add_custom_command(OUTPUT ${CMAKE_INSTALL_PREFIX}/lib/libqphix_codegen.a
    COMMAND CXX=${host_cxx} CXXFLAGS=${host_cxxflags}
      cmake
      -Disa=${isa}
      -Dtarget_cxx=${CMAKE_CXX_COMPILER}
      -Dtarget_cxxflags=${CMAKE_CXX_FLAGS}
      -Dtarget_jN=${recursive_jN}
      -DCMAKE_INSTALL_PREFIX=${CMAKE_INSTALL_PREFIX}
      -DPYTHON_EXECUTABLE=${PYTHON_EXECUTABLE}
      -DPYTHON_LIBRARY=${PYTHON_LIBRARY}
      -DPYTHON_INCLUDE_DIR=${PYTHON_INCLUDE_DIR}
      -Dskip_build=${skip_build}
      -G ${CMAKE_GENERATOR}
      ${CMAKE_SOURCE_DIR}/codegen
    COMMAND ${CMAKE_MAKE_PROGRAM} ${make_program_arguments} VERBOSE=1
                   COMMENT "Generating Libcodegen"
                   WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/codegen )

  #
  # now add a target which depends on this file existing. 
  # This is kind of liek an alias. I can then make later targets depend on this
  #
  add_custom_target(build_codegen DEPENDS ${CMAKE_INSTALL_PREFIX}/lib/libqphix_codegen.a )

  #
  # I want to be able to use the libqphix_codegen.a  library as proper CMake library
  # so this defines a library as IMPORTED
  add_library(qphix_codegen STATIC IMPORTED)

  # This says where the library file is. So now adding
  # `target_link_libraries(foo libcodegen)` will link exeutable `foo` against
  # `libqphix_codegen.a`.
  set_property(TARGET qphix_codegen PROPERTY IMPORTED_LOCATION
    ${CMAKE_INSTALL_PREFIX}/lib/libqphix_codegen.a)

  # I will also make `libcodegen` depend on the `build_codegen` target so if
  # something depends on `libcodegen`, that will trigger the `build_codegen`
  # target which will cause the `libqphix_codegen.a` library to be rebuilt if
  # it does not yet exist.
  add_dependencies(qphix_codegen build_codegen)
  include_directories(${CMAKE_INSTALL_PREFIX}/include)
endif()

###############################################################################

set(test_runner_args -by 4 -bz 4 -c 2 -sy 1 -sz 1 -pxy 1 -pxyz 0 -minct 1 -x 64 -y 16 -z 8 -t 8 -prec f -geom 1 1 1 2 -soalen ${test_soalen})

###############################################################################
configure_file(
  ${PROJECT_SOURCE_DIR}/include/qphix/qphix_config_internal.h.in.cmake
  ${PROJECT_BINARY_DIR}/include/qphix/qphix_config_internal.h 
)

add_subdirectory(include/qphix)
add_subdirectory(include/qphix_codegen)
add_subdirectory(lib)
if( ${testing} )
  # GTest lives here
  add_subdirectory(external)
  add_subdirectory(tests-gtest)
endif()
add_subdirectory(tests)

install(FILES ${PROJECT_SOURCE_DIR}/cmake/QPhiXConfig.cmake DESTINATION cmake)

# vim: sts=2 sw=2
