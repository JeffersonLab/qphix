SUBDIRS = ./avx ./avx2 ./mic ./scalar ./qpx

qphix_includedir=$(includedir)/qphix

qphix_include_HEADERS= \
	./qphix_config.h \
	./dslash_utils.h

qphix_include_HEADERS += \
	./abs_solver.h \
	./geometry.h \
	./comm.h \
	./dslash_def.h \
	./dslash_body.h \
	./dslash_generated.h \
	./face.h \
	./Barrier_mic.h \
	./Barrier_stubs.h \
	./wilson.h \
	./invcg.h \
	./invbicgstab.h \
	./inv_richardson_multiprec.h \
	./minvcg.h \
	./linearOp.h \
	./memmap.h \
	./tsc.h \
	./blas_utils.h \
	./real_functors.h \
	./complex_functors.h \
	./site_loops.h \
	./blas_new_c.h \
	./print_utils.h \
	./threadbind.h


#
# not used currently: blas.h blas_c.h blas_mic.h
#

if QPHIX_BUILD_WITH_QDP
qphix_include_HEADERS += ./qdp_packer.h
endif


if QPHIX_BUILD_CLOVER
qphix_include_HEADERS += ./clover.h \
	./clover_dslash_def.h\
        ./clover_dslash_body.h \
        ./clov_dslash_generated.h \
	./clov_face.h
endif

nodist_qphix_include_HEADERS=./qphix_config_internal.h

