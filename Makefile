#
# CoDiPack, a Code Differentiation Package
#
# Copyright (C) 2015-2018 Chair for Scientific Computing (SciComp), TU Kaiserslautern
# Homepage: http://www.scicomp.uni-kl.de
# Contact:  Prof. Nicolas R. Gauger (codi@scicomp.uni-kl.de)
#
# Lead developers: Max Sagebaum, Tim Albring (SciComp, TU Kaiserslautern)
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
# Authors: Max Sagebaum, Tim Albring, (SciComp, TU Kaiserslautern)
#

# names of the basic deriectories
BUILD_DIR = build
DOC_DIR   = documentation

#list all source files in TEST_DIR
DOC_FILES   = $(wildcard $(DOC_DIR)/*.cpp)

#list all dependency files in BUILD_DIR
DEP_FILES   = $(wildcard $(BUILD_DIR)/*.d)

CODI_DIR := .

FLAGS = -Wall -pedantic -std=c++11 -DCODI_OptIgnoreInvalidJacobies=true -DCODI_EnableAssert=true -I$(CODI_DIR)/include -fopenmp

ifeq ($(OPT), yes)
  CXX_FLAGS := -O3 $(FLAGS)
else
  CXX_FLAGS := -O0 -g $(FLAGS)
endif

ifeq ($(CXX), )
	CXX := g++
else
	CXX := $(CXX)
endif

# Complete list of test files
TUTORIALS = $(patsubst $(DOC_DIR)/%.cpp,$(BUILD_DIR)/%.exe,$(DOC_FILES))

# set default rule
tutorials:

$(BUILD_DIR)/%.exe : $(DOC_DIR)/%.cpp
	@mkdir -p $(@D)
	$(CXX) $(CXX_FLAGS) $< -o $@
	@$(CXX) $(CXX_FLAGS) $< -MM -MP -MT $@ -MF $@.d

tutorials: $(TUTORIALS)
	@mkdir -p $(BUILD_DIR)

.PHONY: clean
clean:
	rm -fr $(BUILD_DIR)

-include $(DEP_FILES)
