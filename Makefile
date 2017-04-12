# User facing Makefile.

SHELL = /bin/bash

make_codegen = -f codegen.mak

.PHONY: mic avx avx2 avx512 sse scalar

all: mic avx avx2 avx512 sse scalar

mic: mic_dir \
	 	mic_PRECISION-2_SOALEN-8 mic_PRECISION-2_SOALEN-4 mic_PRECISION-1_SOALEN-16 mic_PRECISION-1_SOALEN-8 mic_PRECISION-1_SOALEN-4 mic_PRECISION-1_SOALEN-16_ENABLE_LOW_PRECISION-1 mic_PRECISION-1_SOALEN-8_ENABLE_LOW_PRECISION-1 mic_PRECISION-1_SOALEN-4_ENABLE_LOW_PRECISION-1

mic_dir:
	mkdir -p generated/mic/generated

mic_PRECISION-2_SOALEN-8:
	$(MAKE) $(make_codegen) mode=mic PRECISION=2 SOALEN=8

mic_PRECISION-2_SOALEN-4:
	$(MAKE) $(make_codegen) mode=mic PRECISION=2 SOALEN=4

mic_PRECISION-1_SOALEN-16:
	$(MAKE) $(make_codegen) mode=mic PRECISION=1 SOALEN=16

mic_PRECISION-1_SOALEN-8:
	$(MAKE) $(make_codegen) mode=mic PRECISION=1 SOALEN=8

mic_PRECISION-1_SOALEN-4:
	$(MAKE) $(make_codegen) mode=mic PRECISION=1 SOALEN=4

mic_PRECISION-1_SOALEN-16_ENABLE_LOW_PRECISION-1:
	$(MAKE) $(make_codegen) mode=mic PRECISION=1 SOALEN=16 ENABLE_LOW_PRECISION=1

mic_PRECISION-1_SOALEN-8_ENABLE_LOW_PRECISION-1:
	$(MAKE) $(make_codegen) mode=mic PRECISION=1 SOALEN=8 ENABLE_LOW_PRECISION=1

mic_PRECISION-1_SOALEN-4_ENABLE_LOW_PRECISION-1:
	$(MAKE) $(make_codegen) mode=mic PRECISION=1 SOALEN=4 ENABLE_LOW_PRECISION=1

avx512: avx512_dir \
		mic_AVX512-1_PRECISION-2_SOALEN-8 mic_AVX512-1_PRECISION-2_SOALEN-4 mic_AVX512-1_PRECISION-1_SOALEN-16 mic_AVX512-1_PRECISION-1_SOALEN-8 mic_AVX512-1_PRECISION-1_SOALEN-4 mic_AVX512-1_PRECISION-1_SOALEN-16_ENABLE_LOW_PRECISION-1 mic_AVX512-1_PRECISION-1_SOALEN-8_ENABLE_LOW_PRECISION-1 mic_AVX512-1_PRECISION-1_SOALEN-4_ENABLE_LOW_PRECISION-1

avx512_dir:
	mkdir -p generated/avx512/generated

mic_AVX512-1_PRECISION-2_SOALEN-8:
	$(MAKE) $(make_codegen) mode=mic AVX512=1 PRECISION=2 SOALEN=8

mic_AVX512-1_PRECISION-2_SOALEN-4:
	$(MAKE) $(make_codegen) mode=mic AVX512=1 PRECISION=2 SOALEN=4

mic_AVX512-1_PRECISION-1_SOALEN-16:
	$(MAKE) $(make_codegen) mode=mic AVX512=1 PRECISION=1 SOALEN=16

mic_AVX512-1_PRECISION-1_SOALEN-8:
	$(MAKE) $(make_codegen) mode=mic AVX512=1 PRECISION=1 SOALEN=8

mic_AVX512-1_PRECISION-1_SOALEN-4:
	$(MAKE) $(make_codegen) mode=mic AVX512=1 PRECISION=1 SOALEN=4

mic_AVX512-1_PRECISION-1_SOALEN-16_ENABLE_LOW_PRECISION-1:
	$(MAKE) $(make_codegen) mode=mic AVX512=1 PRECISION=1 SOALEN=16 ENABLE_LOW_PRECISION=1

mic_AVX512-1_PRECISION-1_SOALEN-8_ENABLE_LOW_PRECISION-1:
	$(MAKE) $(make_codegen) mode=mic AVX512=1 PRECISION=1 SOALEN=8 ENABLE_LOW_PRECISION=1

mic_AVX512-1_PRECISION-1_SOALEN-4_ENABLE_LOW_PRECISION-1:
	$(MAKE) $(make_codegen) mode=mic AVX512=1 PRECISION=1 SOALEN=4 ENABLE_LOW_PRECISION=1

avx: avx_dir \
		avx_PRECISION-2_SOALEN-2 avx_PRECISION-2_SOALEN-4 avx_PRECISION-1_SOALEN-8 avx_PRECISION-1_SOALEN-4

avx_dir:
	mkdir -p generated/avx/generated

avx_PRECISION-2_SOALEN-2:
	$(MAKE) $(make_codegen) mode=avx PRECISION=2 SOALEN=2

avx_PRECISION-2_SOALEN-4:
	$(MAKE) $(make_codegen) mode=avx PRECISION=2 SOALEN=4

avx_PRECISION-1_SOALEN-8:
	$(MAKE) $(make_codegen) mode=avx PRECISION=1 SOALEN=8

avx_PRECISION-1_SOALEN-4:
	$(MAKE) $(make_codegen) mode=avx PRECISION=1 SOALEN=4

avx2: avx2_dir \
		avx_PRECISION-2_SOALEN-2_AVX2-1 avx_PRECISION-2_SOALEN-4_AVX2-1 avx_PRECISION-1_SOALEN-8_AVX2-1 avx_PRECISION-1_SOALEN-4_AVX2-1 avx_PRECISION-1_SOALEN-8_AVX2-1_ENABLE_LOW_PRECISION-1 avx_PRECISION-1_SOALEN-4_AVX2-1_ENABLE_LOW_PRECISION-1

avx2_dir:
	mkdir -p generated/avx2/generated

avx_PRECISION-2_SOALEN-2_AVX2-1:
	$(MAKE) $(make_codegen) mode=avx PRECISION=2 SOALEN=2 AVX2=1

avx_PRECISION-2_SOALEN-4_AVX2-1:
	$(MAKE) $(make_codegen) mode=avx PRECISION=2 SOALEN=4 AVX2=1

avx_PRECISION-1_SOALEN-8_AVX2-1:
	$(MAKE) $(make_codegen) mode=avx PRECISION=1 SOALEN=8 AVX2=1

avx_PRECISION-1_SOALEN-4_AVX2-1:
	$(MAKE) $(make_codegen) mode=avx PRECISION=1 SOALEN=4 AVX2=1

avx_PRECISION-1_SOALEN-8_AVX2-1_ENABLE_LOW_PRECISION-1:
	$(MAKE) $(make_codegen) mode=avx PRECISION=1 SOALEN=8 AVX2=1 ENABLE_LOW_PRECISION=1

avx_PRECISION-1_SOALEN-4_AVX2-1_ENABLE_LOW_PRECISION-1:
	$(MAKE) $(make_codegen) mode=avx PRECISION=1 SOALEN=4 AVX2=1 ENABLE_LOW_PRECISION=1

sse: sse_dir sse_PRECISION-2 sse_PRECISION-1

sse_dir:
	mkdir -p generated/sse/generated

sse_PRECISION-2:
	$(MAKE) $(make_codegen) mode=sse PRECISION=2

sse_PRECISION-1:
	$(MAKE) $(make_codegen) mode=sse PRECISION=1

scalar: scalar_dir scalar_PRECISION-2 scalar_PRECISION-1

scalar_dir:
	mkdir -p generated/scalar/generated

scalar_PRECISION-2:
	$(MAKE) $(make_codegen) mode=scalar PRECISION=2

scalar_PRECISION-1:
	$(MAKE) $(make_codegen) mode=scalar PRECISION=1

clean: 
	rm -rf build

cleanall: clean
	rm -rf generated
