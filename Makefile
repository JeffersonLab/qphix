# User facing Makefile.

SHELL = /bin/bash

make_codegen := $(MAKE) -f codegen.mak

.PHONY: mic avx avx2 avx512 sse scalar

all: mic avx avx2 avx512 sse scalar

mic:
	mkdir -p generated/mic/generated
	make clean && $(make_codegen) mode=mic PRECISION=2 SOALEN=8 && ./codegen
	make clean && $(make_codegen) mode=mic PRECISION=2 SOALEN=4 && ./codegen
	make clean && $(make_codegen) mode=mic PRECISION=1 SOALEN=16 && ./codegen
	make clean && $(make_codegen) mode=mic PRECISION=1 SOALEN=8 && ./codegen
	make clean && $(make_codegen) mode=mic PRECISION=1 SOALEN=4 && ./codegen
	make clean && $(make_codegen) mode=mic PRECISION=1 SOALEN=16 ENABLE_LOW_PRECISION=1 && ./codegen
	make clean && $(make_codegen) mode=mic PRECISION=1 SOALEN=8 ENABLE_LOW_PRECISION=1 && ./codegen
	make clean && $(make_codegen) mode=mic PRECISION=1 SOALEN=4 ENABLE_LOW_PRECISION=1 && ./codegen

avx512:
	mkdir -p generated/avx512/generated
	make clean && $(make_codegen) mode=mic AVX512=1 PRECISION=2 SOALEN=8 && ./codegen
	make clean && $(make_codegen) mode=mic AVX512=1 PRECISION=2 SOALEN=4 && ./codegen
	make clean && $(make_codegen) mode=mic AVX512=1 PRECISION=1 SOALEN=16 && ./codegen
	make clean && $(make_codegen) mode=mic AVX512=1 PRECISION=1 SOALEN=8 && ./codegen
	make clean && $(make_codegen) mode=mic AVX512=1 PRECISION=1 SOALEN=4 && ./codegen
	make clean && $(make_codegen) mode=mic AVX512=1 PRECISION=1 SOALEN=16 ENABLE_LOW_PRECISION=1 && ./codegen
	make clean && $(make_codegen) mode=mic AVX512=1 PRECISION=1 SOALEN=8 ENABLE_LOW_PRECISION=1 && ./codegen
	make clean && $(make_codegen) mode=mic AVX512=1 PRECISION=1 SOALEN=4 ENABLE_LOW_PRECISION=1 && ./codegen

avx:
	mkdir -p generated/avx/generated
	make clean && $(make_codegen) mode=avx PRECISION=2 SOALEN=2 && ./codegen
	make clean && $(make_codegen) mode=avx PRECISION=2 SOALEN=4 && ./codegen
	make clean && $(make_codegen) mode=avx PRECISION=1 SOALEN=8 && ./codegen
	make clean && $(make_codegen) mode=avx PRECISION=1 SOALEN=4 && ./codegen

avx2:
	mkdir -p generated/avx2/generated
	make clean && $(make_codegen) mode=avx PRECISION=2 SOALEN=2 AVX2=1 && ./codegen
	make clean && $(make_codegen) mode=avx PRECISION=2 SOALEN=4 AVX2=1 && ./codegen
	make clean && $(make_codegen) mode=avx PRECISION=1 SOALEN=8 AVX2=1 && ./codegen
	make clean && $(make_codegen) mode=avx PRECISION=1 SOALEN=4 AVX2=1 && ./codegen
	make clean && $(make_codegen) mode=avx PRECISION=1 SOALEN=8 AVX2=1 ENABLE_LOW_PRECISION=1 && ./codegen
	make clean && $(make_codegen) mode=avx PRECISION=1 SOALEN=4 AVX2=1 ENABLE_LOW_PRECISION=1 && ./codegen

sse:
	mkdir -p generated/sse/generated
	make clean && $(make_codegen) mode=sse PRECISION=2 && ./codegen
	make clean && $(make_codegen) mode=sse PRECISION=1 && ./codegen

scalar:
	mkdir -p generated/scalar/generated
	make clean && $(make_codegen) mode=scalar PRECISION=2 && ./codegen
	make clean && $(make_codegen) mode=scalar PRECISION=1 && ./codegen

clean: 
	rm -rf *.o ./codegen 

cleanall: 
	rm -rf *.o ./codegen
	rm -rf generated
