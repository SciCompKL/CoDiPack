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

CODI_DIR := .

FLAGS = -Wall -pedantic -DCODI_OptIgnoreInvalidJacobies=true -DCODI_EnableAssert=true -I$(CODI_DIR)/include -fopenmp

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

ifdef CXX
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
	doxygen

.PHONY: clean
clean:
	rm -fr $(BUILD_DIR)

.PHONY: force
$(BUILD_DIR)/compiler_flags: force
	@mkdir -p $(@D)
	@echo '$(FLAGS)' | cmp -s - $@ || echo '$(FLAGS)' > $@

-include $(DEP_FILES)
