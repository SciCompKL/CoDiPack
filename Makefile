# names of the basic directories
BUILD_DIR = build
DOC_DIR   = documentation

#list all source files in DOC_DIR
DOC_FILES   = $(wildcard $(DOC_DIR)/*.cpp)

#list all dependency files in BUILD_DIR
DEP_FILES   = $(wildcard $(BUILD_DIR)/*.d)

CODI_DIR := .

FLAGS = -Wall -pedantic -DCODI_OptIgnoreInvalidJacobies=true -DCODI_EnableAssert=true -I$(CODI_DIR)/include -fopenmp

ifeq ($(CPP14), yes)
  FLAGS += -std=c++14
else
  FLAGS += -std=c++11
endif
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

doc:
	@mkdir -p $(BUILD_DIR)
	doxygen

.PHONY: clean
clean:
	rm -fr $(BUILD_DIR)

-include $(DEP_FILES)
