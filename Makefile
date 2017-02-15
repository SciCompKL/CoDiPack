# names of the basic deriectories
BUILD_DIR = build
DOC_DIR   = documentation

#list all source files in TEST_DIR
DOC_FILES   = $(wildcard $(DOC_DIR)/*.cpp)

#list all dependency files in BUILD_DIR
DEP_FILES   = $(wildcard $(BUILD_DIR)/*.d)

CODI_DIR := .

FLAGS = -Wall -pedantic -std=c++11 -DCODI_OptIgnoreInvalidJacobies=true -DCODI_EnableAssert=true -I$(CODI_DIR)/include

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
