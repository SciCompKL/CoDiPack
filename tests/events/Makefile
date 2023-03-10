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

CODI_DIR = ../..

# regards CXX already set
CXX ?= g++

# set to yes to compile with optimization flags
OPT ?= no

# regards CXXFLAGS already set
FLAGS = $(CXXFLAGS)

ifeq ($(OPT),no)
	FLAGS += -O0 -ggdb
else
	FLAGS += -O3
endif

# vector dimension used for vector mode tests
VECTOR_DIM ?= 2

# whether there should be events for the inner tape in second order tests
USE_INNER_CALLBACKS = no

# default target
all:

# delete executables and output files
.PHONY: clean
clean:
	rm -rf build/*

# disable deletion of intermediate targets
.SECONDARY:

# basic types covered by the tests
basicTypes = RealForward RealReverse RealReverseIndex RealReversePrimal RealReversePrimalIndex

# identify primal types among basic types
primalTypes = $(foreach type,$(basicTypes),$(if $(findstring Primal,$(type)),$(type),))

# primal types with variable adjoint interface
primalTypesVariableAdjointInterface = $(foreach type,$(primalTypes),$(type)VariableAdjointInterface)

# all scalar types
scalarTypes = $(basicTypes) $(primalTypesVariableAdjointInterface)

# vector types derived from scalar types
vectorTypes = $(foreach type,$(scalarTypes),$(type)Vec)

# second order types derived from scalar types
secondOrderTypes = $(foreach type,$(scalarTypes),$(type)Sec)

allTypes = $(scalarTypes) $(vectorTypes) $(secondOrderTypes)

allTests = Tape Statement Preacc

$(info [info] all types: $(allTypes))
$(info [info] all tests: $(allTests))

# user selection of types
TYPES ?= $(allTypes)

# user selection of tests
TESTS ?= $(allTests)

# filter against known types
usedTypes = $(filter $(allTypes),$(TYPES))
$(if $(filter-out $(allTypes),$(TYPES)),$(error Unknown type(s) $(filter-out $(allTypes),$(TYPES))))

# filter against known tests
usedTests = $(filter $(allTests),$(TESTS))
$(if $(filter-out $(allTests),$(TESTS)),$(error Unknown test(s) $(filter-out $(allTests),$(TESTS))))

$(info [info] used types: $(usedTypes))
$(info [info] used tests: $(usedTests))

# test cases are combinations of used types and tests
usedCases = $(foreach type,$(usedTypes),$(foreach test,$(usedTests),$(type)$(test)))

# given a type, decide whether it needs the forward, the reverse, or the reversePrimal driver
getDriver = $(if $(findstring ReversePrimal,$1),ReversePrimal,$(if $(findstring Reverse,$1),Reverse,Forward))

# given a type, decide whether it uses a forward or reverse tape, and produce a corresponding preprocessor flag
getTapeFlag = $(if $(findstring Reverse,$1),-DREVERSE_TAPE,-DFORWARD_TAPE)

# given a type, decide whether it is a second order type, and if configured, request events for the inner tape
ifeq ($(USE_INNER_CALLBACKS),yes)
getInnerCallbacksFlag = $(if $(findstring Sec,$1),-DUSE_INNER_CALLBACKS)
else
getInnerCallbacksFlag = 
endif

# to not clutter the result files of second order tests with primal value vector allocations, reduce the chunk sizes
# this is especially important for primal value tapes, but for consistency, it is applied everywhere
getChunkSizeFlags = -DCODI_ChunkSize=1 -DCODI_SmallChunkSize=1

# if needed, produce a flag that enables variable adjoint interfaces
getVariableAdjointInterfaceFlag = $(if $(findstring VariableAdjointInterface,$1),-DCODI_VariableAdjointInterfaceInPrimalTapes)

# given a type, remove the indicators from the name that are not relevant for deriving the CoDiPack type
removeSecondaryIndicators = $(subst VariableAdjointInterface,,$1)

# given a type, add the vector dimension if needed (e.g., RealReverseVec -> RealReverseVec<4>)
addVectorDimIfNeeded = $(patsubst %Vec,%Vec<$(VECTOR_DIM)>,$1)

# given a type, expand it to a second order CoDiPackType type (e.g., RealReverseSec -> RealReverseGen<RealForward>)
expandSecondOrderTypeIfNeeded = $(patsubst %Sec,%Gen<codi::RealForward>,$1)

# given a type, convert it into the corresponding CoDiPack type
getCoDiPackType = codi::$(call addVectorDimIfNeeded,$(call expandSecondOrderTypeIfNeeded,$(call removeSecondaryIndicators,$1)))

# define variables that describe the characteristics of a test case

# test-specific variables
define setTest
build/%$1.exe: TEST_NAME=$1
build/%$1.exe: TEST_FLAGS=
endef

$(foreach test,$(usedTests),$(eval $(call setTest,$(test))))

# type-specific variables
define setType
build/$1%.exe: TYPE=$1
build/$1%.exe: DRIVER=$(call getDriver,$1)
build/$1%.exe: TYPE_FLAGS=$(call getTapeFlag,$1) $(call getInnerCallbacksFlag,$1) $(call getChunkSizeFlags,$1) $(call getVariableAdjointInterfaceFlag,$1)
endef

$(foreach type,$(usedTypes),$(eval $(call setType,$(type))))

# compile tests for a specific type, create executable
build/%.exe:
	@mkdir -p build;
	$(CXX) src/driver$(DRIVER).cpp -o $@ $(FLAGS) -DNUMBER='$(call getCoDiPackType,$(TYPE))' -I $(CODI_DIR)/include -DCODI_ADWorkflowEvents -DCODI_PreaccEvents -DCODI_StatementEvents -DCODI_IndexEvents -DTEST_NAME=Test$(TEST_NAME) $(TYPE_FLAGS) $(TEST_FLAGS) $(DRIVER_FLAGS);
	@$(CXX) src/driver$(DRIVER).cpp $(FLAGS) -DNUMBER='$(call getCoDiPackType,$(TYPE))' -I $(CODI_DIR)/include -DCODI_ADWorkflowEvents -DCODI_PreaccEvents -DCODI_StatementEvents -DCODI_IndexEvents -DTEST_NAME=Test$(TEST_NAME) $(TYPE_FLAGS) $(TEST_FLAGS) $(DRIVER_FLAGS) -MM -MP -MT $@ -MF $@.d

# run executable, create output file
build/%.out: build/%.exe
	./$< > $@

# compare output file against reference file
build/%.diff: build/%.out
	@if cmp -s $< results/$*.ref; then printf "%s \e[0;32mOK\e[0m\n" "$*"; else printf "%s \e[0;31mFAILED\e[0m\n" "$*"; diff $< results/$*.ref || true; fi

# delete reference files
.PHONY: cleanReference
cleanReference: clean
	rm -f results/*

# run executable, create reference file
results/%.ref: build/%.exe
	./$< > $@

# create all reference files
.PHONY: reference
reference: $(foreach case,$(usedCases),results/$(case).ref)

# run all tests
.PHONY: all
all: $(foreach case,$(usedCases),build/$(case).diff)

DEPENDENCIES = $(shell find build -name '*.d')

-include $(DEPENDENCIES)