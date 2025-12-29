#
# CoDiPack, a Code Differentiation Package
#
# Copyright (C) 2015-2026 Chair for Scientific Computing (SciComp), University of Kaiserslautern-Landau
# Homepage: http://scicomp.rptu.de
# Contact:  Prof. Nicolas R. Gauger (codi@scicomp.uni-kl.de)
#
# Lead developers: Max Sagebaum, Johannes Blühdorn (SciComp, University of Kaiserslautern-Landau)
#
# This file is part of CoDiPack (http://scicomp.rptu.de/software/codi).
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

CODI_DIR = ../..

# regards CXX already set
CXX ?= g++

# set to yes to compile with optimization flags
OPT ?= no

# regards CXXFLAGS already set
FLAGS = $(CXXFLAGS) -std=c++17

ifeq ($(OPT),no)
	FLAGS += -O0 -ggdb
else
	FLAGS += -O3
endif

# vector dimension used for vector mode tests
VECTOR_DIM ?= 2

# whether there should be events for the inner tape in second order tests
USE_INNER_CALLBACKS = no

# general flags for OpenMP tests
OPEN_MP_FLAGS=-DREVERSE_TAPE -DCODI_EnableOpenMP -fopenmp -DCODI_SpecializeStdComplex=0

# default target
all:

# delete executables and output files
.PHONY: clean
clean:
	rm -rf build/*

# disable deletion of intermediate targets
.SECONDARY:

# to not clutter the result files of second order tests with primal value vector allocations, reduce the chunk sizes
# this is especially important for primal value tapes, but for consistency, it is applied everywhere
commonFlags = -DCODI_ChunkSize=1 -DCODI_SmallChunkSize=1

# type-specific variables
define setType
build/$1%.exe: TYPE=$1
build/$1%.exe: CODI_TYPE=$2
build/$1%.exe: DRIVER=$3
build/$1%.exe: TYPE_FLAGS=$4 $(commonFlags)
allTypes := $(allTypes) $1
endef

# scalar types

$(eval $(call setType,RealForward,codi::RealForward,Forward,-DFORWARD_TAPE))
$(eval $(call setType,RealReverse,codi::RealReverse,Reverse,-DREVERSE_TAPE))
$(eval $(call setType,RealReverseIndex,codi::RealReverseIndex,Reverse,-DREVERSE_TAPE))
$(eval $(call setType,RealReverseIndexOpenMP,codi::RealReverseIndexOpenMP,Reverse,$(OPEN_MP_FLAGS)))
$(eval $(call setType,RealReversePrimal,codi::RealReversePrimal,ReversePrimal,-DREVERSE_TAPE))
$(eval $(call setType,RealReversePrimalIndex,codi::RealReversePrimalIndex,ReversePrimal,-DREVERSE_TAPE))

$(eval $(call setType,RealReversePrimalVariableAdjointInterface,codi::RealReversePrimal,ReversePrimal,-DREVERSE_TAPE -DCODI_VariableAdjointInterfaceInPrimalTapes))
$(eval $(call setType,RealReversePrimalIndexVariableAdjointInterface,codi::RealReversePrimalIndex,ReversePrimal,-DREVERSE_TAPE -DCODI_VariableAdjointInterfaceInPrimalTapes))

# vector types

$(eval $(call setType,RealForwardVec,codi::RealForwardVec<$(VECTOR_DIM)>,Forward,-DFORWARD_TAPE))
$(eval $(call setType,RealReverseVec,codi::RealReverseVec<$(VECTOR_DIM)>,Reverse,-DREVERSE_TAPE))
$(eval $(call setType,RealReverseIndexVec,codi::RealReverseIndexVec<$(VECTOR_DIM)>,Reverse,-DREVERSE_TAPE))
$(eval $(call setType,RealReverseIndexVecOpenMP,codi::RealReverseIndexVecOpenMP<$(VECTOR_DIM)>,Reverse,$(OPEN_MP_FLAGS)))
$(eval $(call setType,RealReversePrimalVec,codi::RealReversePrimalVec<$(VECTOR_DIM)>,ReversePrimal,-DREVERSE_TAPE))
$(eval $(call setType,RealReversePrimalIndexVec,codi::RealReversePrimalIndexVec<$(VECTOR_DIM)>,ReversePrimal,-DREVERSE_TAPE))

$(eval $(call setType,RealReversePrimalVariableAdjointInterfaceVec,codi::RealReversePrimalVec<$(VECTOR_DIM)>,ReversePrimal,-DREVERSE_TAPE -DCODI_VariableAdjointInterfaceInPrimalTapes))
$(eval $(call setType,RealReversePrimalIndexVariableAdjointInterfaceVec,codi::RealReversePrimalIndexVec<$(VECTOR_DIM)>,ReversePrimal,-DREVERSE_TAPE -DCODI_VariableAdjointInterfaceInPrimalTapes))

# set flag for inner callbacks, depending on config
ifeq ($(USE_INNER_CALLBACKS),yes)
innerCallbacksFlag = -DUSE_INNER_CALLBACKS
else
innerCallbacksFlag = 
endif

# second order types

$(eval $(call setType,RealForwardSec,codi::RealForwardGen<codi::RealForward>,Forward,-DFORWARD_TAPE $(innerCallbacksFlag)))
$(eval $(call setType,RealReverseSec,codi::RealReverseGen<codi::RealForward>,Reverse,-DREVERSE_TAPE $(innerCallbacksFlag)))
$(eval $(call setType,RealReverseIndexSec,codi::RealReverseIndexGen<codi::RealForward>,Reverse,-DREVERSE_TAPE $(innerCallbacksFlag)))
$(eval $(call setType,RealReverseIndexSecOpenMP,codi::RealReverseIndexOpenMPGen<codi::RealForward>,Reverse,$(OPEN_MP_FLAGS) $(innerCallbacksFlag)))
$(eval $(call setType,RealReversePrimalSec,codi::RealReversePrimalGen<codi::RealForward>,ReversePrimal,-DREVERSE_TAPE $(innerCallbacksFlag)))
$(eval $(call setType,RealReversePrimalIndexSec,codi::RealReversePrimalIndexGen<codi::RealForward>,ReversePrimal,-DREVERSE_TAPE $(innerCallbacksFlag)))

$(eval $(call setType,RealReversePrimalVariableAdjointInterfaceSec,codi::RealReversePrimalGen<codi::RealForward>,ReversePrimal,-DREVERSE_TAPE -DCODI_VariableAdjointInterfaceInPrimalTapes $(innerCallbacksFlag)))
$(eval $(call setType,RealReversePrimalIndexVariableAdjointInterfaceSec,codi::RealReversePrimalIndexGen<codi::RealForward>,ReversePrimal,-DREVERSE_TAPE -DCODI_VariableAdjointInterfaceInPrimalTapes $(innerCallbacksFlag)))

# tests

define setTest
build/%$1.exe: TEST_NAME=$1
build/%$1.exe: TEST_FLAGS=$2
allTests := $(allTests) $1
endef

$(eval $(call setTest,Tape,))
$(eval $(call setTest,Statement,))
$(eval $(call setTest,Preacc,))
$(eval $(call setTest,Aggregates,))

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

# compile tests for a specific type, create executable
build/%.exe:
	@mkdir -p build;
	$(CXX) src/driver$(DRIVER).cpp -o $@ $(FLAGS) -DNUMBER='$(CODI_TYPE)' -I $(CODI_DIR)/include -DCODI_ADWorkflowEvents -DCODI_PreaccEvents -DCODI_StatementEvents -DCODI_IndexEvents -DTEST_NAME=Test$(TEST_NAME) $(TYPE_FLAGS) $(TEST_FLAGS);
	@$(CXX) src/driver$(DRIVER).cpp $(FLAGS) -DNUMBER='$(CODI_TYPE)' -I $(CODI_DIR)/include -DCODI_ADWorkflowEvents -DCODI_PreaccEvents -DCODI_StatementEvents -DCODI_IndexEvents -DTEST_NAME=Test$(TEST_NAME) $(TYPE_FLAGS) $(TEST_FLAGS) -MM -MP -MT $@ -MF $@.d

# run executable, create output file
build/%.out: build/%.exe
	./$< > $@

# compare output file against reference file
build/%.compare: build/%.out
	@if cmp -s $< results/$*.ref; then printf "%s \e[0;32mOK\e[0m\n" "$*"; else printf "%s \e[0;31mFAILED\e[0m\n" "$*"; diff $< results/$*.ref | tee build/$*.diff; fi

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
all: $(foreach case,$(usedCases),build/$(case).compare)
	@[ -z "$$(find build -name '*.diff')" ];

DEPENDENCIES = $(shell find build -name '*.d')

-include $(DEPENDENCIES)
