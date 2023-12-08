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
BUILD_DIR = build
CODI_DIR := ../..
TEST_DIR = src

TEST_FILES  = $(wildcard $(TEST_DIR)/*.cpp)
DEP_FILES= $(shell find $(BUILD_DIR) -name '*.d')
ALL_PROBLEMS = $(patsubst %.cpp,$(BUILD_DIR)/%.run,$(TEST_FILES))

ifeq ($(TESTS),)
  SELECTED_TESTS = ALL
  SELECTED_PROBLEMS = $(ALL_PROBLEMS)
else
  SELECTED_TESTS = $(TESTS)
  SELECTED_PROBLEMS = $(patsubst %/, $(BUILD_DIR)/src/%.run, $(dir $(TESTS)))
endif

# disable the deletion of secondary targets
.SECONDARY:

FLAGS = -Wall -Werror=return-type -pedantic -DCODI_OptIgnoreInvalidJacobians=true -DCODI_EnableAssert=true -I$(CODI_DIR)/include -fopenmp $(EIGEN_DEFINE) -DCODI_StatementEvents

ifeq ($(OPT), no)
  FLAGS += -O0 -g
else
  FLAGS += -O3
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

all: functionTests

.PHONY: force
$(BUILD_DIR)/compiler_flags: force
	@mkdir -p $(@D)
	@echo '$(FLAGS)' | cmp -s - $@ || echo '$(FLAGS)' > $@

$(BUILD_DIR)/%.exe : %.cpp $(BUILD_DIR)/compiler_flags
	@mkdir -p $(@D)
	$(CXX) $(FLAGS) $< -o $@
	@$(CXX) $(FLAGS) $< -MM -MP -MT $@ -MF $@.d

%.run : TEST_NAME = $(basename $(<F))
%.run : %.exe $(BUILD_DIR)/compiler_flags
	@echo "Running $(TEST_NAME)"
	@mkdir -p $(<D)/$(TEST_NAME)_run
	cd $(<D)/$(TEST_NAME)_run; ../$(<F) $(SELECTED_TESTS)
	@echo "Comparing"
	@scripts/compare.sh $(TEST_NAME) $(filter ALL $(TEST_NAME)/%, $(SELECTED_TESTS))

functionTests: $(SELECTED_PROBLEMS)

.PHONY: clean
clean:
	rm -fr $(BUILD_DIR)

-include $(DEP_FILES)
