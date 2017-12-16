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
TEST_DIR = tests
BUILD_DIR = build
DRIVER_DIR = drivers
RESULT_BASE_DIR = results_base
RESULT_DIR = results

#list all source files in TEST_DIR
TEST_FILES   = $(wildcard $(TEST_DIR)/Test**.cpp)
TEST_FILES  += $(wildcard $(TEST_DIR)/**/Test**.cpp)
TEST_FILES  += $(wildcard $(TEST_DIR)/**/**/Test**.cpp)

#list all dependency files in BUILD_DIR
DEP_FILES   = $(wildcard $(BUILD_DIR)/Test**.d)
DEP_FILES  += $(wildcard $(BUILD_DIR)/**/Test**.d)
DEP_FILES  += $(wildcard $(BUILD_DIR)/**/**/Test**.d)

FLAGS = -Wall -Wextra -pedantic -std=c++11 -DCODI_OptIgnoreInvalidJacobies=true -DCODI_EnableAssert=true

# The default is to run all drives
DRIVERS?=ALL

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

CODI_DIR := ..

# Complete list of test files
TESTS = $(patsubst $(TEST_DIR)/%.cpp,$(RESULT_DIR)/%.test,$(TEST_FILES))

# set default rule
all:

# disable the deletion of secondary targets
.SECONDARY:

### define general sets for tests
# Tests for all tapes
BASIC_TESTS = $(wildcard $(TEST_DIR)/basic/Test**.cpp) $(wildcard $(TEST_DIR)/expressions/Test**.cpp) $(wildcard $(TEST_DIR)/exceptions/Test**.cpp) $(wildcard $(TEST_DIR)/jacobiOptimizations/Test**.cpp) $(wildcard $(TEST_DIR)/helpers/Test**.cpp)
# Tests that run only for reverse mode tapes
REVERSE_TESTS = $(wildcard $(TEST_DIR)/external_functions/Test**.cpp) $(wildcard $(TEST_DIR)/io/Test**.cpp)
# Tests that run for non vector mode tapes
REVERSE_VALUE_TESTS = $(wildcard $(TEST_DIR)/preaccumulation/Test**.cpp)

# The build rules for all drivers.
define DRIVER_RULE
$(BUILD_DIR)/%_$(DRIVER_NAME)_bin : $(TEST_DIR)/%.cpp
	@mkdir -p $(@D)
	$(CXX) $(CXX_FLAGS) $< -o $@ $(DRIVER)
	@$(CXX) $(CXX_FLAGS) $(DRIVER) $< -MM -MP -MT $@ -MF $@.d
DRIVER_TEST_FILES := $(patsubst $(TEST_DIR)/%.cpp,$(RESULT_DIR)/%.test,$(DRIVER_TESTS))
$(DRIVER_TEST_FILES): $(RESULT_DIR)/%.test: $(RESULT_DIR)/%_$(DRIVER_NAME).out
endef

# The general template for each driver. It checks if the driver should be run
define DRIVER_INST
  ifeq (${DRIVER_NAME}, $(filter ${DRIVER_NAME}, ${DRIVERS}))
    $(eval $(value DRIVER_RULE))
  else
    ifeq (ALL, $(filter ALL, ${DRIVERS}))
      $(eval $(value DRIVER_RULE))
    endif
  endif
endef

# The first 3 lines of each driver define the basic parameters for this driver.
# DRIVER_NAME:    The short name for the driver. This used to create the specific files for the driver
# DRIVER_TESTS:   The full list of tests which are run for this driver. See the general test list for details.
# DRIVER:         The special flags for this driver. It needs to be defined inside the build rule.
# The eval line creates the rules and dependencies for the driver

# Driver for RealForward
DRIVER_NAME  := FWD
DRIVER_TESTS := $(BASIC_TESTS)
$(BUILD_DIR)/%_$(DRIVER_NAME)_bin : DRIVER = $(DRIVER_DIR)/forward/forwardDriver.cpp -I$(CODI_DIR)/include -I$(CODI_DIR)/source -I$(DRIVER_DIR)/forward
$(eval $(value DRIVER_INST))

# Driver for RealForwardVector
DRIVER_NAME  := FWD_Vec
DRIVER_TESTS := $(BASIC_TESTS)
$(BUILD_DIR)/%_$(DRIVER_NAME)_bin : DRIVER = $(DRIVER_DIR)/forwardVector/forwardDriver.cpp -I$(CODI_DIR)/include -I$(CODI_DIR)/source -I$(DRIVER_DIR)/forwardVector
$(eval $(value DRIVER_INST))

# Driver for 2nd order type but first derivative evaluation both forward.
DRIVER_NAME  := FWD2nd
DRIVER_TESTS := $(BASIC_TESTS)
$(BUILD_DIR)/%_$(DRIVER_NAME)_bin : DRIVER = $(DRIVER_DIR)/forward2ndOrder/forwardDriver.cpp -I$(CODI_DIR)/include -I$(CODI_DIR)/source -I$(DRIVER_DIR)/forward2ndOrder
$(eval $(value DRIVER_INST))

# Driver for RealReverseUnchecked
DRIVER_NAME  := RWS_Unch
DRIVER_TESTS := $(BASIC_TESTS) $(REVERSE_TESTS) $(REVERSE_VALUE_TESTS)
$(BUILD_DIR)/%_$(DRIVER_NAME)_bin : DRIVER := $(DRIVER_DIR)/reverseSimple/reverseDriver.cpp -I$(CODI_DIR)/include -I$(DRIVER_DIR)/reverseSimple
$(eval $(value DRIVER_INST))

# Driver for RealReverseIndexUncheckedIndex
DRIVER_NAME  := RWS_UnchInd
DRIVER_TESTS := $(BASIC_TESTS) $(REVERSE_TESTS) $(REVERSE_VALUE_TESTS)
$(BUILD_DIR)/%_$(DRIVER_NAME)_bin : DRIVER := $(DRIVER_DIR)/reverseSimpleIndex/reverseDriver.cpp -I$(CODI_DIR)/include -I$(DRIVER_DIR)/reverseSimpleIndex
$(eval $(value DRIVER_INST))

# Driver for RealReversePrimalUnchecked
DRIVER_NAME  := RWS_PrimUnch
DRIVER_TESTS := $(BASIC_TESTS) $(REVERSE_TESTS) $(REVERSE_VALUE_TESTS)
$(BUILD_DIR)/%_$(DRIVER_NAME)_bin : DRIVER := $(DRIVER_DIR)/reversePrimalSimple/reverseDriver.cpp -I$(CODI_DIR)/include -I$(DRIVER_DIR)/reversePrimalSimple
$(eval $(value DRIVER_INST))

# Driver for RealReversePrimal
DRIVER_NAME  := RWS_Prim
DRIVER_TESTS := $(BASIC_TESTS) $(REVERSE_TESTS) $(REVERSE_VALUE_TESTS)
$(BUILD_DIR)/%_$(DRIVER_NAME)_bin : DRIVER := $(DRIVER_DIR)/reversePrimal/reverseDriver.cpp -I$(CODI_DIR)/include -I$(DRIVER_DIR)/reversePrimal
$(eval $(value DRIVER_INST))

# Driver for RealReversePrimalVec
DRIVER_NAME  := RWS_PrimVec
DRIVER_TESTS := $(BASIC_TESTS) $(REVERSE_TESTS)
$(BUILD_DIR)/%_$(DRIVER_NAME)_bin : DRIVER := $(DRIVER_DIR)/reversePrimalVector/reverseDriver.cpp -I$(CODI_DIR)/include -I$(DRIVER_DIR)/reversePrimalVector
$(eval $(value DRIVER_INST))

# Driver for RealReversePrimalVectorAdapter
DRIVER_NAME  := RWS_PrimVecA
DRIVER_TESTS := $(BASIC_TESTS) $(REVERSE_TESTS)
$(BUILD_DIR)/%_$(DRIVER_NAME)_bin : DRIVER := $(DRIVER_DIR)/reversePrimalVectorAdapter/reverseDriver.cpp -I$(CODI_DIR)/include -I$(DRIVER_DIR)/reversePrimalVectorAdapter -DCODI_EnableVariableAdjointInterfaceInPrimalTapes
$(eval $(value DRIVER_INST))

# Driver for RealReversePrimalIndex
DRIVER_NAME  := RWS_PrimIndex
DRIVER_TESTS := $(BASIC_TESTS) $(REVERSE_TESTS) $(REVERSE_VALUE_TESTS)
$(BUILD_DIR)/%_$(DRIVER_NAME)_bin : DRIVER := $(DRIVER_DIR)/reversePrimalIndex/reverseDriver.cpp -I$(CODI_DIR)/include -I$(DRIVER_DIR)/reversePrimalIndex
$(eval $(value DRIVER_INST))

# Driver for RealReverse
DRIVER_NAME  := RWS_Chunk
DRIVER_TESTS := $(BASIC_TESTS) $(REVERSE_TESTS) $(REVERSE_VALUE_TESTS)
$(BUILD_DIR)/%_$(DRIVER_NAME)_bin : DRIVER = $(DRIVER_DIR)/reverseChunk/reverseDriver.cpp -I$(CODI_DIR)/include -I$(DRIVER_DIR)/reverseChunk
$(eval $(value DRIVER_INST))

# Driver for RealReverseVector
DRIVER_NAME  := RWS_ChunkVec
DRIVER_TESTS := $(BASIC_TESTS) $(REVERSE_TESTS)
$(BUILD_DIR)/%_$(DRIVER_NAME)_bin : DRIVER = $(DRIVER_DIR)/reverseChunkVector/reverseDriver.cpp -I$(CODI_DIR)/include -I$(DRIVER_DIR)/reverseChunkVector
$(eval $(value DRIVER_INST))

# Driver for RealReverseVectorAdapter
DRIVER_NAME  := RWS_ChunkVecA
DRIVER_TESTS := $(BASIC_TESTS) $(REVERSE_TESTS)
$(BUILD_DIR)/%_$(DRIVER_NAME)_bin : DRIVER = $(DRIVER_DIR)/reverseChunkVectorAdapter/reverseDriver.cpp -I$(CODI_DIR)/include -I$(DRIVER_DIR)/reverseChunkVectorAdapter
$(eval $(value DRIVER_INST))

# Driver for RealReverseIndex
DRIVER_NAME  := RWS_ChunkInd
DRIVER_TESTS := $(BASIC_TESTS) $(REVERSE_TESTS) $(REVERSE_VALUE_TESTS)
$(BUILD_DIR)/%_$(DRIVER_NAME)_bin : DRIVER = $(DRIVER_DIR)/reverseChunkIndex/reverseDriver.cpp -I$(CODI_DIR)/include -I$(DRIVER_DIR)/reverseChunkIndex
$(eval $(value DRIVER_INST))

# Driver for RealReverseIndex with tape swap
DRIVER_NAME  := RWS_ChunkIndSwap
DRIVER_TESTS := $(BASIC_TESTS) $(REVERSE_TESTS) $(REVERSE_VALUE_TESTS)
$(BUILD_DIR)/%_$(DRIVER_NAME)_bin : DRIVER = $(DRIVER_DIR)/reverseChunkIndexTapeSwap/reverseDriver.cpp -I$(CODI_DIR)/include -I$(DRIVER_DIR)/reverseChunkIndexTapeSwap
$(eval $(value DRIVER_INST))

# Driver for RealReverseIndexVector
DRIVER_NAME  := RWS_ChunkIndVec
DRIVER_TESTS := $(BASIC_TESTS) $(REVERSE_TESTS)
$(BUILD_DIR)/%_$(DRIVER_NAME)_bin : DRIVER = $(DRIVER_DIR)/reverseChunkIndexVector/reverseDriver.cpp -I$(CODI_DIR)/include -I$(DRIVER_DIR)/reverseChunkIndexVector
$(eval $(value DRIVER_INST))

# Driver for reverse over forward
DRIVER_NAME  := RWS2nd
DRIVER_TESTS := $(BASIC_TESTS) $(REVERSE_TESTS)
$(BUILD_DIR)/%_$(DRIVER_NAME)_bin : DRIVER = $(DRIVER_DIR)/reverseOverForward/reverseOverForwardDriver.cpp -I$(CODI_DIR)/include -I$(DRIVER_DIR)/reverseOverForward
$(eval $(value DRIVER_INST))

# rules for generating the test files
$(RESULT_DIR)/%.out : $(BUILD_DIR)/%_bin
	@mkdir -p $(@D)
	$< > $@

# rule for printing the results (dependencies are generated by the drivers)
$(RESULT_DIR)/%.test:
	@./compare.sh -n $* -b $(RESULT_BASE_DIR)/$*.out $^

all: $(TESTS)
	@mkdir -p $(BUILD_DIR)

.PHONY: clean
clean:
	rm -fr $(BUILD_DIR)
	rm -fr $(RESULT_DIR)

-include $(DEP_FILES)
