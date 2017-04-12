# User facing Makefile.

SHELL = /bin/bash

make_codegen := $(MAKE) -f codegen.mak

.PHONY: mic avx avx2 avx512 sse scalar

all: mic avx avx2 avx512 sse scalar

mic:
	mkdir -p generated/mic/generated
	$(make_codegen) mode=mic PRECISION=2 SOALEN=8
	$(make_codegen) mode=mic PRECISION=2 SOALEN=4
	$(make_codegen) mode=mic PRECISION=1 SOALEN=16
	$(make_codegen) mode=mic PRECISION=1 SOALEN=8
	$(make_codegen) mode=mic PRECISION=1 SOALEN=4
	$(make_codegen) mode=mic PRECISION=1 SOALEN=16 ENABLE_LOW_PRECISION=1
	$(make_codegen) mode=mic PRECISION=1 SOALEN=8 ENABLE_LOW_PRECISION=1
	$(make_codegen) mode=mic PRECISION=1 SOALEN=4 ENABLE_LOW_PRECISION=1

avx512:
	mkdir -p generated/avx512/generated
	$(make_codegen) mode=mic AVX512=1 PRECISION=2 SOALEN=8
	$(make_codegen) mode=mic AVX512=1 PRECISION=2 SOALEN=4
	$(make_codegen) mode=mic AVX512=1 PRECISION=1 SOALEN=16
	$(make_codegen) mode=mic AVX512=1 PRECISION=1 SOALEN=8
	$(make_codegen) mode=mic AVX512=1 PRECISION=1 SOALEN=4
	$(make_codegen) mode=mic AVX512=1 PRECISION=1 SOALEN=16 ENABLE_LOW_PRECISION=1
	$(make_codegen) mode=mic AVX512=1 PRECISION=1 SOALEN=8 ENABLE_LOW_PRECISION=1
	$(make_codegen) mode=mic AVX512=1 PRECISION=1 SOALEN=4 ENABLE_LOW_PRECISION=1

avx:
	mkdir -p generated/avx/generated
	$(make_codegen) mode=avx PRECISION=2 SOALEN=2
	$(make_codegen) mode=avx PRECISION=2 SOALEN=4
	$(make_codegen) mode=avx PRECISION=1 SOALEN=8
	$(make_codegen) mode=avx PRECISION=1 SOALEN=4

avx2:
	mkdir -p generated/avx2/generated
	$(make_codegen) mode=avx PRECISION=2 SOALEN=2 AVX2=1
	$(make_codegen) mode=avx PRECISION=2 SOALEN=4 AVX2=1
	$(make_codegen) mode=avx PRECISION=1 SOALEN=8 AVX2=1
	$(make_codegen) mode=avx PRECISION=1 SOALEN=4 AVX2=1
	$(make_codegen) mode=avx PRECISION=1 SOALEN=8 AVX2=1 ENABLE_LOW_PRECISION=1
	$(make_codegen) mode=avx PRECISION=1 SOALEN=4 AVX2=1 ENABLE_LOW_PRECISION=1

sse:
	mkdir -p generated/sse/generated
	$(make_codegen) mode=sse PRECISION=2
	$(make_codegen) mode=sse PRECISION=1

scalar:
	mkdir -p generated/scalar/generated
	$(make_codegen) mode=scalar PRECISION=2
	$(make_codegen) mode=scalar PRECISION=1

clean: 
	rm -rf *.o ./codegen 

cleanall: 
	rm -rf *.o ./codegen
	rm -rf generated
