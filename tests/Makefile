#
# CoDiPack, a Code Differentiation Package
#
# Copyright (C) 2015-2021 Chair for Scientific Computing (SciComp), TU Kaiserslautern
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
NOOP=
SPACE = $(NOOP) $(NOOP)
KOMMA = $(NOOP),$(NOOP)

# names of the basic directories
INCLUDE_DIR = include
SRC_DIR = src
BUILD_DIR = build
TEST_DIR = $(INCLUDE_DIR)/tests
DRIVER_DIR = $(INCLUDE_DIR)/drivers
CODI_DIR := ..
CODI_INCLUDE = -I$(CODI_DIR)/include

FLAGS = -Wall -Wextra -pedantic -std=c++11 -DCODI_IgnoreInvalidJacobians=true -DCODI_EnableAssert=true -I$(INCLUDE_DIR)

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

#list all source files in TEST_DIR
TEST_FILES   = $(wildcard $(TEST_DIR)/test**.hpp)
TEST_FILES  += $(wildcard $(TEST_DIR)/**/test**.hpp)
TEST_FILES  += $(wildcard $(TEST_DIR)/**/**/test**.hpp)

#list all dependency files in BUILD_DIR
DEP_FILES   = $(wildcard $(BUILD_DIR)/**.d)
DEP_FILES  += $(wildcard $(BUILD_DIR)/**/**.d)
DEP_FILES  += $(wildcard $(BUILD_DIR)/**/**/**.d)
DEP_FILES  += $(wildcard $(BUILD_DIR)/**/**/**/**.d)

# Test selection by user
ALL_TESTS = $(patsubst test%.hpp,Test%,$(notdir $(TEST_FILES)))
ifeq ($(TESTS),)
  CONFID = $(BUILD_DIR)/allTests.conf
  SELECTED_TESTS = $(ALL_TESTS)
else
  CONFID = $(BUILD_DIR)/$(subst $(SPACE),_,$(strip $(SELECTED_TESTS))).conf
  SELECTED_TESTS = $(patsubst %,Test%,$(TESTS))
endif

ifeq ($(UPDATE), yes)
  UPDATE_FLAG=-u
endif

# ALL_DRIVERS will contain the complete list of all drivers
ALL_DRIVERS =

# set default rule
all:

# disable the deletion of secondary targets
.SECONDARY:

#include drivers
include Makefile.drivers

# selection of drivers to run
ifeq ($(DRIVERS),)
  SELECTED_DRIVERS = $(ALL_DRIVERS)
else
  SELECTED_DRIVERS = $(DRIVERS)
endif

DRIVER_RUNS = $(patsubst %,$(BUILD_DIR)/%.run, $(SELECTED_DRIVERS))

# regular build rules
$(BUILD_DIR)/%.o : $(SRC_DIR)/%.cpp
	@mkdir -p $(@D)
	$(CXX) $(CXX_FLAGS) $< -c -o $@
	@$(CXX) $(CXX_FLAGS) $< -MM -MP -MT $@ -MF $@.d

$(BUILD_DIR)/%.exe :
	@mkdir -p $(@D)
	$(CXX) $(CXX_FLAGS) $^ -o $@

# Update rule for dependency on TESTS parameter
CREATE_CONFID := $(shell mkdir -p $(BUILD_DIR) && touch -a $(CONFID))

$(BUILD_DIR)/conf.d: $(CONFID)
	mkdir -p $(BUILD_DIR)
	rm -f $(BUILD_DIR)/*.conf
	rm -rf $(BUILD_DIR)/results
	touch $(CONFID)
	touch $@

# build rules for test runners
define createTests
	$(subst $(SPACE),$(KOMMA),$(filter $1,$(SELECTED_TESTS)))
endef

$(BUILD_DIR)/%.runner.o : $(SRC_DIR)/runner.cpp $(BUILD_DIR)/conf.d
	@mkdir -p $(@D)
	echo "#define DRIVER_TESTS $(call createTests,$(DRIVER_TESTS))" > $@.tests.hpp
	$(CXX) $(CXX_FLAGS) $(DRIVER_INC) -I. -DDRIVER_TESTS_INC='"$@.tests.hpp"' $< -c -o $@
	@$(CXX) $(CXX_FLAGS) $(DRIVER_INC) -I. -DDRIVER_TESTS_INC='"$@.tests.hpp"' $< -MM -MP -MT $@ -MF $@.d

$(BUILD_DIR)/%.runner.exe : $(BUILD_DIR)/%.runner.o
	@mkdir -p $(@D)
	$(CXX) $(CXX_FLAGS) $(DRIVER_LIB) $^ -o $@

$(BUILD_DIR)/%.run : $(BUILD_DIR)/%.runner.exe
	@mkdir -p $(@D)
	$< | tee $@

# rule for all tests file generation
$(BUILD_DIR)/allTests.hpp: $(INCLUDE_DIR)/tests/allTests.hpp
	echo "#define ALL_TESTS $(call createTests,$(ALL_TESTS))" > $@

# output dependencies and rule
$(BUILD_DIR)/listAllNames.o : $(BUILD_DIR)/allTests.hpp
$(BUILD_DIR)/compare.exe : $(BUILD_DIR)/compare.o $(BUILD_DIR)/compareFiles.o $(BUILD_DIR)/listAllNames.o


.PHONY: output
output: $(BUILD_DIR)/compare.exe $(DRIVER_RUNS)
	$(BUILD_DIR)/compare.exe -t 1e-14 $(patsubst %, -d %, $(SELECTED_DRIVERS))

.PHONY: all
all: output

.PHONY: clean
clean:
	rm -fr $(BUILD_DIR)

-include $(DEP_FILES)
