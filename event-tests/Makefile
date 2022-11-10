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

CODI_DIR = ..

# regards CXX already set
CXX ?= g++

# set to yes to compile with debugging flags
DEBUG = no

# regards CXXFLAGS already set
FLAGS = $(CXXFLAGS)

ifeq ($(DEBUG),yes)
	FLAGS += -O0 -ggdb
else
	FLAGS += -O3
endif

# default target
all:

# disable deletion of intermediate targets
.SECONDARY:

# scalar types covered by the tests
scalarTypes = RealForward RealReverse RealReverseIndex RealReversePrimal RealReversePrimalIndex

# vector types derived from scalar types
vectorTypes = $(foreach type,$(scalarTypes),$(type)Vec)

allTypes = $(scalarTypes) $(vectorTypes)

allTests = Tape Statement Preacc

# test cases are all combinations of types and tests
allCases = $(foreach type,$(allTypes),$(foreach test,$(allTests),$(type)$(test)))

# ensure that the build directory exists
.PHONY: prepare
prepare:
	@mkdir -p build

# delete executables and output files
.PHONY: clean
clean:
	rm -f build/*

HEADERS = $(shell find include -name '*.hpp')

# given a type, decide whether it needs the forward or the reverse driver
getDriver = $(if $(findstring Reverse,$1),Reverse,Forward)

# given a type, add the vector dimension if needed (e.g., RealReverseVec -> RealReverseVec<4>)
addVectorDimIfNeeded = $(if $(findstring Vec,$1),$1<4>,$1)

# define variables that describe the characteristics of a test case

define setTest
build/%$1.exe: TEST=$1
endef

$(foreach test,$(allTests),$(eval $(call setTest,$(test))))

define setType
build/$1%.exe: TYPE=$1
build/$1%.exe: DRIVER=$(call getDriver,$1)
endef

$(foreach type,$(allTypes),$(eval $(call setType,$(type))))

# compile tests for a specific type, create executable
build/%.exe: src/driver$(call getDriver,$*).cpp $(HEADERS) prepare
	$(CXX) src/driver$(DRIVER).cpp -o $@ $(FLAGS) -DNUMBER='codi::$(call addVectorDimIfNeeded,$(TYPE))' -I $(CODI_DIR)/include -DCODI_ADWorkflowEvents -DCODI_PreaccEvents -DCODI_StatementEvents -DCODI_IndexEvents -DCODI_VariableAdjointInterfaceInPrimalTapes -DTEST=$(TEST)

# run executable, create output file
build/%.out: build/%.exe
	./$< > $@

# compare output file against reference file
build/%.diff: build/%.out
	@if cmp -s $< results/$*.ref; then echo -e $* "\e[0;32mOK\e[0m"; else echo -e $* "\e[0;31mFAILED\e[0m"; diff $< results/$*.ref || true; fi

# delete reference files
.PHONY: cleanReference
cleanReference:
	rm -f results/*

# run executable, create reference file
results/%.ref: build/%.exe
	./$< > $@

# create all reference files
.PHONY: reference
reference: $(foreach case,$(allCases),results/$(case).ref)

# run all tests
.PHONY: all
all: $(foreach case,$(allCases),build/$(case).diff)
