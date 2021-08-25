# names of the basic directories
BUILD_DIR    = build
DOC_DIR      = documentation
EXAMPLE_DIR  = $(DOC_DIR)/examples
TUTORIAL_DIR = $(DOC_DIR)/tutorials
DEVELOPER_DIR = $(DOC_DIR)/developer

#list all source files in DOC_DIR
TUTORIAL_FILES  = $(wildcard $(TUTORIAL_DIR)/*.cpp)
EXAMPLE_FILES  = $(wildcard $(EXAMPLE_DIR)/*.cpp) $(wildcard $(DEVELOPER_DIR)/*.cpp)

#list all dependency files in BUILD_DIR
DEP_FILES   = $(wildcard $(BUILD_DIR)/*.d)
DEP_FILES  += $(wildcard $(BUILD_DIR)/**/*.d)
DEP_FILES  += $(wildcard $(BUILD_DIR)/**/**/*.d)

MAJOR_VERSION = $(shell grep -oP 'define CODI_MAJOR_VERSION \K\d+' include/codi.hpp)
MINOR_VERSION = $(shell grep -oP 'define CODI_MINOR_VERSION \K\d+' include/codi.hpp)
BUILD_VERSION = $(shell grep -oP 'define CODI_BUILD_VERSION \K\d+' include/codi.hpp)
CODI_VERSION = $(MAJOR_VERSION).$(MINOR_VERSION).$(BUILD_VERSION)

CODI_DIR := .

FLAGS = -Wall -pedantic -DCODI_OptIgnoreInvalidJacobians=true -DCODI_EnableAssert=true -I$(CODI_DIR)/include -fopenmp

ifndef CLANG_FORMAT
  CLANG_FORMAT := clang-format
else
  CLANG_FORMAT := $(CLANG_FORMAT)
endif

ifeq ($(CPP14), yes)
  FLAGS += -std=c++14
else
  FLAGS += -std=c++11
endif
ifeq ($(OPT), yes)
  FLAGS += -O3
else
  FLAGS += -O0 -g
endif
ifeq ($(MPI), yes)
  FLAGS += -DCODI_EnableMPI
  ifdef MEDI_DIR
    FLAGS += -I$(MEDI_DIR)/include -I$(MEDI_DIR)/src
  else
    $(error Error: 'MEDI_DIR' is not defined for the MPI build. You can get it at 'https://www.scicomp.uni-kl.de/software/medi/' or with 'git clone https://github.com/SciCompKL/MeDiPack.git'.)
  endif
endif

ifndef CXX
  ifeq ($(MPI), yes)
    CXX := mpic++
  else
    CXX := g++
  endif
else
  CXX := $(CXX)
endif

TUTORIALS = $(patsubst %.cpp,$(BUILD_DIR)/%.exe,$(TUTORIAL_FILES))
EXAMPLES = $(patsubst %.cpp,$(BUILD_DIR)/%.exe,$(EXAMPLE_FILES))

# set default rule
all: tutorials examples

$(BUILD_DIR)/%.exe : %.cpp $(BUILD_DIR)/compiler_flags
	@mkdir -p $(@D)
	$(CXX) $(FLAGS) $< -o $@
	@$(CXX) $(FLAGS) $< -MM -MP -MT $@ -MF $@.d

tutorials: $(TUTORIALS)

examples: $(EXAMPLES)

doc:
	@mkdir -p $(BUILD_DIR)
	@mkdir -p $(BUILD_DIR)/documentation
	CODI_VERSION=$(CODI_VERSION) doxygen

.PHONY: format
format:
	find include tests/include tests/src -type f -exec $(CLANG_FORMAT) -i {} \;

.PHONY: clean
clean:
	rm -fr $(BUILD_DIR)

.PHONY: force
$(BUILD_DIR)/compiler_flags: force
	@mkdir -p $(@D)
	@echo '$(FLAGS)' | cmp -s - $@ || echo '$(FLAGS)' > $@

-include $(DEP_FILES)
