cmake_minimum_required(VERSION 3.1)
project(QPhiX-Codegen)

function(setup_defs target precision soalen lowprec)

    if( ${SERIAL_SPIN} )  
    	list(APPEND yesnolist SERIAL_SPIN ) 
    endif()
    if( ${PREF_L1_SPINOR_IN} ) 
    	list( APPEND yesnolist PREF_L1_SPINOR_IN)  
    endif()
    if( ${PREF_L2_SPINOR_IN} )
    	 list( APPEND yesnolist PREF_L2_SPINOR_IN)  
    endif()
    if( ${PREF_L1_SPINOR_OUT} ) 
    	list( APPEND yesnolist PREF_L1_SPINOR_OUT)  
   endif()
    if( ${PREF_L2_SPINOR_OUT} ) 
    	list( APPEND yesnolist PREF_L2_SPINOR_OUT)  
    endif()
    if( ${PREF_L1_GAUGE} )
    	 list( APPEND yesnolist PREF_L1_GAUGE)  
    endif()
    if( ${PREF_L2_GAUGE} )
    	 list( APPEND yesnolist PREF_L2_GAUGE) 
    endif()
    if( ${PREF_L1_CLOVER} ) 
    	list( APPEND yesnolist PREF_L1_CLOVER)  
    endif()
    if( ${PREF_L2_CLOVER} )
    	 list( APPEND yesnolist PREF_L2_CLOVER)  
    endif()
    if( ${USE_LDUNPK} )
    	 list( APPEND yesnolist USE_LDUNPK)  
    endif()
    if( ${USE_PKST} ) 
    	list( APPEND yesnolist USE_PKST)  
    endif()
    if( ${USE_SHUFFLES} ) 
    	list( APPEND yesnolist USE_SHUFFLES) 
    endif()
    if( ${USE_NO_GPREF_L1} )
    	 list( APPEND yesnolist USE_NO_GPREF_L2)  
    endif()
    if( ${USE_NO_GPREF_L2} ) 
    	list( APPEND yesnolist USE_NO_GPREF_L2) 
    endif()
    if( ${ENABLE_STREAMING_STORES} ) 
    	list( APPEND yesnolist ENABLE_STREAMING_STORES)  
    endif()
    if( ${USE_PACKED_GAUGES} ) 
    	list( APPEND yesnolist USE_PACKED_GAUGES)  
    endif()
    if( ${USE_PACKED_CLOVER} ) 
    	list( APPEND yesnolist USE_PACKED_CLOVER)  
    endif()
    if( ${USE_NO_HW_MASKING} ) 
    	list( APPEND yesnolist NO_HW_MASKING ) 
    endif()
    if( ${USE_NO_MASKS} ) 
    	list( APPEND yesnolist NO_MASKS) 
    endif()

   #Enable AVX2 or AVX512
    if( ${isa} STREQUAL "avx2" )
      LIST(APPEND deflist AVX2 )
    elseif ( ${isa} STREQUAL "avx512" )
      LIST(APPEND deflist AVX512)
    endif()

    # Set Low precision related defines
    if( ${lowprec} )
       LIST(APPEND deflist ENABLE_LOW_PRECISION)
       LIST(APPEND deflist USE_LP_SPINOR)
       LIST(APPEND deflist USE_LP_GAUGE)
       LIST(APPEND deflist USE_LP_CLOVER)
    endif()

    LIST(APPEND deflist SOALEN=${soalen})

   LIST(APPEND deflist PRECISION=${precision})
   if( ${precision} EQUAL 1 )

   # Single precision vector lengths
      if (  ( ${isa} STREQUAL "avx" ) OR ( ${isa} STREQUAL "avx2" ) )
          LIST(APPEND deflist VECLEN=8)
      elseif ( ${isa} STREQUAL "avx512" )
     	  LIST(APPEND deflist VECLEN=16)
      elseif ( ${isa} STREQUAL "sse" )
    	  LIST(APPEND deflist VECLEN=4)
      elseif( ${isa} STREQUAL "scalar" )
          LIST(APPEND deflist VECLEN=1)
      endif()

   else() 

      #Double precision vector lengths
      if (  ( ${isa} STREQUAL "avx" ) OR ( ${isa} STREQUAL "avx2" ) )
            LIST(APPEND deflist VECLEN=4)
      elseif ( ${isa} STREQUAL "avx512" )
      	    LIST(APPEND deflist VECLEN=8)
      elseif ( ${isa} STREQUAL "sse" )
       	    LIST(APPEND deflist VECLEN=2)
      elseif( ${isa} STREQUAL "scalar" )
            LIST(APPEND deflist VECLEN=1)
      endif()

endif()

MESSAGE(INFO, "Yesnolist  is ${yesnolist}")
MESSAGE(INFO, "Deflist2 is ${deflist}")
target_compile_definitions(${target} PRIVATE ${yesnolist} PRIVATE ${deflist} )
endfunction()

LIST(APPEND codegen_sources codegen.cc
	    data_types.cc
	    dslash.cc
	    dslash_common.cc
	    inst_dp_vec8.cc 
	    inst_sp_vec16.cc
 	    inst_dp_vec4.cc 
	    inst_sp_vec8.cc 
	    inst_sp_vec4.cc 
	    inst_dp_vec2.cc 
	    inst_scalar.cc
	    address_types.h 
	    data_types.h 
	    dslash.h 
	    instructions.h)

OPTION( isa  "select ISA to generate for" )

if( ${isa} STREQUAL "avx" OR ${isa} STREQUAL "avx2" )

    set( SERIAL_SPIN TRUE)
    set( PREF_L1_SPINOR_IN FALSE )
    set( PREF_L2_SPINOR_IN FALSE )
    set( PREF_L1_SPINOR_OUT FALSE )
    set( PREF_L2_SPINOR_OUT FALSE )
    set( PREF_L1_GAUGE FALSE )
    set( PREF_L2_GAUGE FALSE )
    set( PREF_L1_CLOVER FALSE )
    set( PREF_L2_CLOVER FALSE )
    set( USE_LDUNPK  FALSE)
    set( USE_PKST, TRUE)
    set( USE_SUFFLES  FALSE)
    set( NO_GPREF_L1  TRUE )
    set( NO_GPREF_L2  TRUE )
    set( ENABLE_STREAMING_STORES  TRUE)
    set( USE_PACKED_GAUGES TRUE)
    set( USE_PACKED_CLOVER TRUE)
    set( NO_HW_MASKING TRUE )
    set( NO_MASKS FALSE) 

elseif( ${isa} STREQUAL "avx512" )

    set( SERIAL_SPIN  TRUE)
    set( PREF_L1_SPINOR_IN  TRUE)
    set( PREF_L2_SPINOR_IN TRUE )
    set( PREF_L1_SPINOR_OUT TRUE )
    set( PREF_L2_SPINOR_OUT TRUE )
    set( PREF_L1_GAUGE  TRUE )
    set( PREF_L2_GAUGE TRUE )
    set( PREF_L1_CLOVER TRUE )
    set( PREF_L2_CLOVER TRUE )
    set( USE_LDUNPK TRUE)
    set( USE_PKST TRUE) 
    set( USE_SHUFFLES  FALSE)
    set( NO_GPREF_L1 TRUE )
    set( NO_GPREF_L2 TRUE )
    set( ENABLE_STREAMING_STORES TRUE)
    set( USE_PACKED_GAUGES TRUE)
    set( USE_PACKED_CLOVER TRUE)
    set( NO_HW_MASKING FALSE )
    set( NO_MASKS FALSE) 

elseif( ${isa} STREQUAL "sse" ) 

    set( SERIAL_SPIN  FALSE)
    set( PREF_L1_SPINOR_IN  FALSE)
    set( PREF_L2_SPINOR_IN  FALSE )
    set( PREF_L1_SPINOR_OUT FALSE )
    set( PREF_L2_SPINOR_OUT FALSE )
    set( PREF_L1_GAUGE  FALSE )
    set( PREF_L2_GAUGE FALSE )
    set( PREF_L1_CLOVER FALSE )
    set( PREF_L2_CLOVER FALSE )
    set( USE_LDUNPK TRUE)
    set( USE_PKST TRUE) 
    set( USE_SHUFFLES  FALSE)
    set( NO_GPREF_L1 TRUE )
    set( NO_GPREF_L2 TRUE )
    set( ENABLE_STREAMING_STORES TRUE)
    set( USE_PACKED_GAUGES TRUE)
    set( USE_PACKED_CLOVER TRUE)
    set( NO_HW_MASKING TRUE )
    set( NO_MASKS TRUE) 

elseif( ${isa} STREQUAL "scalar" )

    set( SERIAL_SPIN  FALSE)
    set( PREF_L1_SPINOR_IN  FALSE)
    set( PREF_L2_SPINOR_IN  FALSE )
    set( PREF_L1_SPINOR_OUT FALSE )
    set( PREF_L2_SPINOR_OUT FALSE )
    set( PREF_L1_GAUGE  FALSE )
    set( PREF_L2_GAUGE FALSE )
    set( PREF_L1_CLOVER FALSE )
    set( PREF_L2_CLOVER FALSE )
    set( USE_LDUNPK TRUE)
    set( USE_PKST TRUE) 
    set( USE_SHUFFLES  FALSE)
    set( NO_GPREF_L1 TRUE )
    set( NO_GPREF_L2 TRUE )
    set( ENABLE_STREAMING_STORES TRUE)
    set( USE_PACKED_GAUGES TRUE)
    set( USE_PACKED_CLOVER TRUE)
    set( NO_HW_MASKING TRUE )
    set( NO_MASKS TRUE) 

endif()


file(MAKE_DIRECTORY  ${CMAKE_CURRENT_BINARY_DIR}/generated/${isa}/generated)

if( ${isa} STREQUAL "avx")

    

    set(target avx_PRECISION-2_SOALEN-2)
    add_executable(${target} ${codegen_sources})
    setup_defs(${target} 2 2 FALSE)
    list(APPEND target_list ${target})

    set(target avx_PRECISION-2_SOALEN-4)
    add_executable(${target} ${codegen_sources})
    setup_defs(${target} 2 4 FALSE)
    list(APPEND target_list ${target})

    set(target avx_PRECISION-1_SOALEN-4)
    add_executable(${target} ${codegen_sources})
    setup_defs(${target} 1 4 FALSE)
    list(APPEND target_list ${target})

    set(target avx_PRECISION-1_SOALEN-8)
    add_executable(${target} ${codegen_sources})
    setup_defs(${target} 1 8 FALSE)
    list(APPEND target_list ${target})

   elseif( ${isa} STREQUAL "avx2" ) 

	set(target avx_PRECISION-2_SOALEN_2_AVX2-1)
	add_executable(${target} ${codegen_sources})
	setup_defs(${target} 2 2 FALSE)
        list(APPEND target_list ${target})
 
	set(target avx_PRECISION-2_SOALEN_4_AVX2-1)
	add_executable(${target} ${codegen_sources})
	setup_defs(${target} 2 4 FALSE)
        list(APPEND target_list ${target})

	set(target avx_PRECISION-1_SOALEN_4_AVX2-1)
	add_executable(${target} ${codegen_sources})
	setup_defs(${target} 1 4 FALSE)
        list(APPEND target_list ${target})

	set(target avx_PRECISION-1_SOALEN_8_AVX2-1)
	add_executable(${target} ${codegen_sources})
	setup_defs(${target} 1 8 FALSE)
        list(APPEND target_list ${target})

	set(target avx_PRECISION-1_SOALEN_4_AVX2-1_ENABLE_LOW_PRECISION-1)
	add_executable(${target} ${codegen_sources})
	setup_defs(${target} 1 4 TRUE)
        list(APPEND target_list ${target})

	set(target avx_PRECISION-1_SOALEN_8_AVX2-1_ENABLE_LOW_PRECISION-1)
	add_executable(${target} ${codegen_sources})
	setup_defs(${target} 1 8 TRUE)
        list(APPEND target_list ${target})

elseif( ${isa} STREQUAL "avx512" )

	set(target mic_PRECISION-2_SOALEN-8)
	add_executable(${target} ${codegen_sources})
	setup_defs(${target} 2 8 FALSE)
        list(APPEND target_list ${target})

	set(target mic_PRECISION-2_SOALEN-4)
	add_executable(${target} ${codegen_sources})
	setup_defs(${target} 2 4 FALSE)
        list(APPEND target_list ${target})
	
	set(target mic_PRECISION-1_SOALEN-16)
	add_executable(${target} ${codegen_sources})
	setup_defs(${target} 1 16  FALSE)
        list(APPEND target_list ${target})

	set(target mic_PRECISION-1_SOALEN-8)
	add_executable(${target} ${codegen_sources})
	setup_defs(${target} 1 8  FALSE)
        list(APPEND target_list ${target})

	set(target mic_PRECISION-1_SOALEN-4)
	add_executable(${target} ${codegen_sources})
	setup_defs(${target} 1 4  FALSE)
        list(APPEND target_list ${target})

	set(target mic_PRECISION-1_SOALEN-16_ENABLE_LOW_PRECISION-1)
	add_executable(${target} ${codegen_sources})
	setup_defs(${target} 1 16 TRUE)
        list(APPEND target_list ${target})

	set(target mic_PRECISION-1_SOALEN-8_ENABLE_LOW_PRECISION-1)
	add_executable(${target} ${codegen_sources})
	setup_defs(${target} 1 8 TRUE)
        list(APPEND target_list ${target})

	set(target mic_PRECISION-1_SOALEN-4_ENABLE_LOW_PRECISION-1)
	add_executable(${target} ${codegen_sources})
	setup_defs(${target} 1 4 TRUE)
        list(APPEND target_list ${target})

elseif( ${isa} STREQUAL "sse" ) 
	set(target, sse_PRECISION-2)
	add_executable(${target} ${codegen_sources})
	setup_defs(${target} 2 2 TRUE)
        list(APPEND target_list ${target})


	set(target, sse_PRECISION-1)
	add_executable(${target} ${codegen_sources})
	setup_defs(${target} 1 2 TRUE)
        list(APPEND target_list ${target})

elseif( ${isa} STREQUAL "scalar" )


	set(target scalar_PRECISION-2)
	add_executable(${target} ${codegen_sources})
	setup_defs(${target} 2 1 FALSE)
        list(APPEND target_list ${target})

	set(target scalar_PRECISION-1)
	add_executable(${target} ${codegen_sources})
	setup_defs(${target} 1 1 FALSE)
        list(APPEND target_list ${target})

else()
	MESSAGE( ERROR, "Unknown isa: ${isa}")
endif()

export(TARGETS ${target_list} FILE qphix-codegen-exports.cmake )

file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/jinja DESTINATION ${CMAKE_CURRENT_BINARY_DIR})




