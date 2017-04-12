# Makefile for the code generator. The user does not call this directly.

SHELL = /bin/bash

mode ?= sentinel

mode:=$(strip $(mode))

CONFFILE=options-$(mode).mak
include $(CONFFILE)

# Check if the Intel C++ compiler is available. If not, use the default CXX of
# the system.
has_icpc = $(shell command -v icpc 2> /dev/null)
ifneq (,$(has_icpc))
    CXX = icpc
endif

CXXHOST = $(CXX) -O3 -g --std=c++11

ifeq ($(mode),mic)
ifeq ($(PRECISION),1)
override VECLEN=16
else
override VECLEN=8
endif
yesnolist += AVX512
endif

ifeq ($(mode),avx)
ifeq ($(PRECISION),1)
override VECLEN=8
else
override VECLEN=4
endif
DEFS += -DNO_HW_MASKING
yesnolist += AVX2
endif

ifeq ($(mode),sse)
ifeq ($(PRECISION),1)
override VECLEN=4
override SOALEN=4
else
override VECLEN=2
override SOALEN=2
endif
DEFS += -DNO_HW_MASKING
yesnolist += NO_MASKS
endif

ifeq ($(mode),scalar)
override VECLEN=1
override SOALEN=1
DEFS += -DNO_HW_MASKING
DEFS += -DNO_MASKS
endif


# If streaming stores are enabled we
# should definitely disable prefetching 
# of output spinors
ifeq ($(ENABLE_STREAMING_STORES),1)
PREF_L1_SPINOR_OUT=0
PREF_L2_SPINOR_OUT=0
endif

ifeq ($(ENABLE_LOW_PRECISION),1)
USE_LP_SPINOR=1
USE_LP_GAUGE=1
USE_LP_CLOVER=1
endif


yesnolist += PREF_L1_SPINOR_IN 
yesnolist += PREF_L2_SPINOR_IN 
yesnolist += PREF_L1_SPINOR_OUT 
yesnolist += PREF_L2_SPINOR_OUT 
yesnolist += PREF_L1_GAUGE
yesnolist += PREF_L2_GAUGE
yesnolist += PREF_L1_CLOVER
yesnolist += PREF_L2_CLOVER
yesnolist += USE_LDUNPK
yesnolist += USE_PKST
yesnolist += USE_PACKED_GAUGES
yesnolist += USE_PACKED_CLOVER
yesnolist += USE_SHUFFLES
yesnolist += NO_GPREF_L1
yesnolist += NO_GPREF_L2
yesnolist += ENABLE_STREAMING_STORES

yesnolist += USE_LP_SPINOR
yesnolist += USE_LP_GAUGE
yesnolist += USE_LP_CLOVER
yesnolist += SERIAL_SPIN
yesnolist += TESTFLAG

deflist += SOALEN
deflist += VECLEN
deflist += PRECISION

DEFS += $(strip $(foreach var, $(yesnolist), $(if $(filter 1, $($(var))), -D$(var))))
DEFS += $(strip $(foreach var, $(deflist), $(if $($(var)), -D$(var)=$($(var)))))

sources = codegen.cc data_types.cc dslash.cc dslash_common.cc inst_dp_vec8.cc inst_sp_vec16.cc inst_dp_vec4.cc inst_sp_vec8.cc inst_sp_vec4.cc inst_dp_vec2.cc inst_scalar.cc
headers = address_types.h data_types.h dslash.h instructions.h Makefile $(CONFFILE)

buildroot := build/codegen-precision$(PRECISION)-soalen$(SOALEN)-low$(ENABLE_LOW_PRECISION)

objects := $(sources:%.cc=$(buildroot)/%.o)

$(buildroot)/codegen: $(sources) $(headers) $(objects) | $(buildroot)
	@if [ "$(mode)" = sentinel ]; then echo 'Fatal error: You have neither specified a mode nor a target. Omitting the target compile the code generator, but that needs a `mode`. If you just want to generate code for some architecture (say avx2), then call `make avx2`.'; exit 1; fi
	$(CXXHOST) $(DEFS) $(objects) -o $@
	# Generate the kernels.
	./$@

$(buildroot)/%.o: %.cc $(headers) | $(buildroot)
	$(CXXHOST) -c $(DEFS) $< -o $@

$(buildroot):
	mkdir -p $@
