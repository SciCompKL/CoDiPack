#
# CoDiPack, a Code Differentiation Package
#
# Copyright (C) 2015-2023 Chair for Scientific Computing (SciComp), University of Kaiserslautern-Landau
# Homepage: http://www.scicomp.uni-kl.de
# Contact:  Prof. Nicolas R. Gauger (codi@scicomp.uni-kl.de)
#
# Lead developers: Max Sagebaum, Johannes Blühdorn (SciComp, University of Kaiserslautern-Landau)
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
#  - SciComp, University of Kaiserslautern-Landau:
#    - Max Sagebaum
#    - Johannes Blühdorn
#    - Former members:
#      - Tim Albring
#

# names of the basic directories
BUILD_DIR      = build
DEFINITION_DIR = definitions
DOC_DIR        = documentation
EXAMPLE_DIR    = $(DOC_DIR)/examples
TUTORIAL_DIR   = $(DOC_DIR)/tutorials
DEVELOPER_DIR  = $(DOC_DIR)/developer

EIGEN_DEFINE=
ifdef EIGEN_DIR
  EIGEN_DEFINE=-I$(EIGEN_DIR) -DCODI_EnableEigen=true
endif

ENZYME_DEFINE=
ifdef ENZYME_DIR
  CLANG_VERSION := 14
  ENZYME_DEFINE = -DCODI_EnableEnzyme=true -flegacy-pass-manager -Xclang -load -Xclang $(ENZYME_DIR)/lib/ClangEnzyme-$(CLANG_VERSION).so -Xclang -plugin-arg-enzyme -Xclang -enzyme-globals-default-inactive=1
endif

#list all source files in DOC_DIR
TUTORIAL_FILES  = $(wildcard $(TUTORIAL_DIR)/*.cpp)
EXAMPLE_FILES  = $(wildcard $(EXAMPLE_DIR)/*.cpp) $(wildcard $(DEVELOPER_DIR)/*.cpp)
DEFINITION_FILES  = $(wildcard $(DEFINITION_DIR)/*/*/*.hpp)

#list all dependency files in BUILD_DIR
DEP_FILES   = $(wildcard $(BUILD_DIR)/*.d)
DEP_FILES  += $(wildcard $(BUILD_DIR)/**/*.d)
DEP_FILES  += $(wildcard $(BUILD_DIR)/**/**/*.d)

MAJOR_VERSION = $(shell grep -oP 'define CODI_MAJOR_VERSION \K\d+' include/codi.hpp)
MINOR_VERSION = $(shell grep -oP 'define CODI_MINOR_VERSION \K\d+' include/codi.hpp)
BUILD_VERSION = $(shell grep -oP 'define CODI_BUILD_VERSION \K\d+' include/codi.hpp)
CODI_VERSION = $(MAJOR_VERSION).$(MINOR_VERSION).$(BUILD_VERSION)

CODI_DIR := .

FLAGS = -Wall -Werror=return-type -pedantic -DCODI_OptIgnoreInvalidJacobians=true -DCODI_EnableAssert=true -I$(CODI_DIR)/include -fopenmp $(EIGEN_DEFINE) $(ENZYME_DEFINE) -DCODI_StatementEvents

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
    $(error Error: 'MEDI_DIR' is not defined for the MPI build. You can get it at 'https://www.scicomp.uni-kl.de/software/medi/' or with 'git clone https://github.com/SciCompKL/MeDiPack.git')
  endif
endif
ifeq ($(OPENMP), yes)
  FLAGS += -DCODI_EnableOpenMP -fopenmp
  ifeq ($(OPDILIB), yes)
    FLAGS += -DCODI_EnableOpDiLib
    ifdef OPDI_DIR
      FLAGS += -I$(OPDI_DIR)/include
    else
      $(error Error: 'OPDI_DIR' is not defined for the OpDiLib build. You can get it at 'https://www.scicomp.uni-kl.de/software/opdi/' or with 'git clone https://github.com/SciCompKL/OpDiLib.git')
    endif
  endif
else
  ifeq ($(OPDILIB), yes)
    $(error Error: OPDILIB=yes requires OPENMP=yes)
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
GENERATED = $(patsubst $(DEFINITION_DIR)/%.hpp,include/codi/tools/%.hpp,$(DEFINITION_FILES))

# set default rule
all: tutorials examples

$(BUILD_DIR)/%.exe : %.cpp $(BUILD_DIR)/compiler_flags
	@mkdir -p $(@D)
	$(CXX) $(FLAGS) $< -o $@
	@$(CXX) $(FLAGS) $< -MM -MP -MT $@ -MF $@.d

tutorials: $(TUTORIALS)

examples: $(EXAMPLES)


include/codi/tools/%.hpp: $(DEFINITION_DIR)/%.hpp
	@mkdir -p $(@D)
	python $(EXT_FUNC_GENERATOR_DIR)/externalfunctionparser/externalfunctionparser.py $< -o $@
	$(CLANG_FORMAT) -i $@

generate: $(GENERATED)

doc:
	@mkdir -p $(BUILD_DIR)
	@mkdir -p $(BUILD_DIR)/documentation
	CODI_VERSION=$(CODI_VERSION) doxygen

.PHONY: format
format:
	find include tests/general/include tests/general/src tests/events/include tests/events/src -type f -exec $(CLANG_FORMAT) -i {} \;

.PHONY: clean
clean:
	rm -fr $(BUILD_DIR)

.PHONY: force
$(BUILD_DIR)/compiler_flags: force
	@mkdir -p $(@D)
	@echo '$(FLAGS)' | cmp -s - $@ || echo '$(FLAGS)' > $@

-include $(DEP_FILES)
