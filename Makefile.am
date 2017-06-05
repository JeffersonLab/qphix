generate_code:
	cd @top_srcdir@/codegen && ./generate-and-compile @ISA@ "$(CXX)" "$(CXXFLAGS)"
	touch generate_code

all: generate_code
all-am: generate_code
all-hook: generate_code
all-recursive: generate_code

SUBDIRS = include lib
if QPHIX_BUILD_TESTS
SUBDIRS += tests #tests-gtest
endif
