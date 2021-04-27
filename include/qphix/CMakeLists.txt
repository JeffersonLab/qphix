list(APPEND header_files qphix_config.h dslash_utils.h
        arith_type.h
        abs_solver.h
        unprec_solver_wrapper.h
        full_spinor.h
        geometry.h
        comm.h
        dslash_def.h
        dslash_body.h
        face.h
        Barrier_mic.h
        Barrier_stubs.h
        wilson.h
        invcg.h
        invbicgstab.h
				invmr.h
        inv_richardson_multiprec.h
        inv_dummy_hermtest.h
        minvcg.h
        linearOp.h
        memmap.h
        tsc.h
        blas_utils.h
        real_functors.h
        complex_functors.h
        kernel_selector.h
        site_loops.h
        blas.h
        blas_c.h
        blas_new_c.h
        blas_full_spinor.h
        print_utils.h
        threadbind.h
        thread_limits.h
        qphix_cli_args.h
        diagnostics.h)

if( ${QDPXX_CXX_FOUND} ) 
	list(APPEND header_files qdp_packer.h qdp_packer_parscalar.h qdp_packer_qdpjit.h)
endif()

if( ${clover} )
	list(APPEND header_files clover.h clover_dslash_def.h clover_dslash_body.h clov_face.h clover_product.h)
endif()

if( ${twisted_mass} ) 
	list(APPEND header_files ndtm_reuse_operator.h twisted_mass.h tm_dslash_def.h tm_dslash_body.h tm_dslash_face.h)
endif()

if( ${tm_clover} )
	list(APPEND header_files ndtm_reuse_operator_clover.h twisted_mass_clover.h tm_clov_dslash_def.h tm_clov_dslash_body.h tm_clov_dslash_face.h)
endif()

list(APPEND header_files ${PROJECT_BINARY_DIR}/include/qphix/qphix_config_internal.h)
install( FILES ${header_files} DESTINATION include/qphix )
