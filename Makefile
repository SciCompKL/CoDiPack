#
# CoDiPack, a Code Differentiation Package
#
# Copyright (C) 2015-2022 Chair for Scientific Computing (SciComp), TU Kaiserslautern
# Homepage: http://www.scicomp.uni-kl.de
# Contact:  Prof. Nicolas R. Gauger (codi@scicomp.uni-kl.de)
#
# Lead developers: Max Sagebaum, Johannes Blühdorn (SciComp, TU Kaiserslautern)
#
# This file is part of CoDiPack (http://www.scicomp.uni-kl.de/software/codi).
#
# CoDiPack is free software: you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation, either version 3 of the
# License, or (at your option) any later version.
#
# CoDiPack is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty
# of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
#
# See the GNU General Public License for more details.
# You should have received a copy of the GNU
# General Public License along with CoDiPack.
# If not, see <http://www.gnu.org/licenses/>.
#
# For other licensing options please contact us.
#
# Authors:
#  - SciComp, TU Kaiserslautern:
#    - Max Sagebaum
#    - Johannes Blühdorn
#    - Former members:
#      - Tim Albring
#

# names of the basic directories
BUILD_DIR    = build
DOC_DIR      = documentation
EXAMPLE_DIR  = $(DOC_DIR)/examples
TUTORIAL_DIR = $(DOC_DIR)/tutorials
DEVELOPER_DIR = $(DOC_DIR)/developer

EIGEN_DEFINE=
ifdef EIGEN_DIR
  EIGEN_DEFINE=-I$(EIGEN_DIR) -DCODI_EnableEigen=true
endif

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

FLAGS = -Wall -Werror=return-type -pedantic -DCODI_OptIgnoreInvalidJacobians=true -DCODI_EnableAssert=true -I$(CODI_DIR)/include -fopenmp $(EIGEN_DEFINE)

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

# disable the deletion of secondary targets
.SECONDARY:

# set default rule
all: tutorials examples algorithmTests

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

# Tests and runners for algorithm tests
TUTORIALS = $(patsubst %.cpp,$(BUILD_DIR)/%.exe,$(TUTORIAL_FILES))
EXAMPLES = $(patsubst %.cpp,$(BUILD_DIR)/%.exe,$(EXAMPLE_FILES))

ALGORITHM_TEST_DIR = tests/include/algorithms
ALGORITHM_TEST_FILES  = $(wildcard $(ALGORITHM_TEST_DIR)/*.cpp)
ALGORITHM_TESTS = $(patsubst %.cpp,$(BUILD_DIR)/%.run,$(ALGORITHM_TEST_FILES))

%.run : TEST_NAME = $(basename $(<F))
%.run : %.exe $(BUILD_DIR)/compiler_flags
	@echo "Running $(TEST_NAME)"
	cd $(<D); ./$(<F)
	@echo "Comparing"
	@diff -ruw tests/resultsAlgorithms/$(TEST_NAME) $(@D)/$(TEST_NAME)

algorithmTests: $(ALGORITHM_TESTS)

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
