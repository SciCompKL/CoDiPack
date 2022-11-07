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

CXX ?= g++

allTypes = RealForward RealReverse RealReverseIndex RealReversePrimal RealReversePrimalIndex

.PHONY: clean
clean:
	rm -f build/*

HEADERS = $(shell find include -name '*.hpp')

getDriver = $(if $(findstring Forward,$1),Forward,Reverse)

build/%.exe: src/driver$(call getDriver,$*).cpp $(HEADERS)
	$(CXX) src/driver$(call getDriver,$*).cpp -o $@ $(CXXFLGS) -DNUMBER=codi::$* -I $(CODI_DIR)/include -DCODI_ADWorkflowEvents -DCODI_PreaccEvents -DCODI_StatementEvents -DCODI_IndexEvents
	
build/%.out: build/%.exe
	./$< > $@
	
build/%.diff: build/%.out
	@if cmp -s $< results/$*.ref; then echo -e $* "\e[0;32mOK\e[0m"; else echo -e $* "\e[0;31mFAILED\e[0m"; diff $< results/$*.ref || true; fi

.PHONY: all
all: $(foreach type,$(allTypes),build/$(type).diff build/$(type).out build/$(type).exe)

.PHONY: cleanReference
cleanReference:
	rm -f results/*

results/%.ref: build/%.exe
	./$< > $@

.PHONY: reference
reference: $(foreach type,$(allTypes),results/$(type).ref build/$(type).exe)
