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

FLAGS = -Wall -pedantic -std=gnu++0x -DCODI_OptIgnoreInvalidJacobies=true -DCODI_EnableAssert=true

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

# define general sets for tests
BASIC_TESTS = $(wildcard $(TEST_DIR)/basic/Test**.cpp) $(wildcard $(TEST_DIR)/expressions/Test**.cpp) $(wildcard $(TEST_DIR)/exceptions/Test**.cpp)
REVERSE_TESTS = $(wildcard $(TEST_DIR)/external_functions/Test**.cpp)

# The first 3 lines of each driver define the basic parameters for this driver.
# DRIVER_NAME:    The short name for the driver. This used to create the specific files for the driver
# DRIVER_TESTS:   The full list of tests which are run for this driver. See the general test list for details.
# DRIVER:         The special flags for this driver. It needs to be defined inside the build rule.
# The other lines create the rules and dependencies to run the tests for this driver.

# Driver for RealForward
DRIVER_NAME  := FWD
DRIVER_TESTS := $(BASIC_TESTS)
$(BUILD_DIR)/%_$(DRIVER_NAME)_bin : DRIVER = $(DRIVER_DIR)/forward/forwardDriver.cpp -I$(CODI_DIR)/include -I$(CODI_DIR)/source -I$(DRIVER_DIR)/forward
$(BUILD_DIR)/%_$(DRIVER_NAME)_bin : $(TEST_DIR)/%.cpp
	@mkdir -p $(@D)
	$(CXX) $(CXX_FLAGS) $< -o $@ $(DRIVER)
	@$(CXX) $(CXX_FLAGS) $(DRIVER) $< -MM -MP -MT $@ -MF $@.d
DRIVER_TEST_FILES := $(patsubst $(TEST_DIR)/%.cpp,$(RESULT_DIR)/%.test,$(DRIVER_TESTS))
$(DRIVER_TEST_FILES): $(RESULT_DIR)/%.test: $(RESULT_DIR)/%_$(DRIVER_NAME).out

# Driver for RealForwardVector
DRIVER_NAME  := FWD_Vec
DRIVER_TESTS := $(BASIC_TESTS)
$(BUILD_DIR)/%_$(DRIVER_NAME)_bin : DRIVER = $(DRIVER_DIR)/forwardVector/forwardDriver.cpp -I$(CODI_DIR)/include -I$(CODI_DIR)/source -I$(DRIVER_DIR)/forwardVector
$(BUILD_DIR)/%_$(DRIVER_NAME)_bin : $(TEST_DIR)/%.cpp
	@mkdir -p $(@D)
	$(CXX) $(CXX_FLAGS) $< -o $@ $(DRIVER)
	@$(CXX) $(CXX_FLAGS) $(DRIVER) $< -MM -MP -MT $@ -MF $@.d
DRIVER_TEST_FILES := $(patsubst $(TEST_DIR)/%.cpp,$(RESULT_DIR)/%.test,$(DRIVER_TESTS))
$(DRIVER_TEST_FILES): $(RESULT_DIR)/%.test: $(RESULT_DIR)/%_$(DRIVER_NAME).out

# Driver for 2nd order type but first derivative evaluation both forward.
DRIVER_NAME  := FWD2nd
DRIVER_TESTS := $(BASIC_TESTS)
$(BUILD_DIR)/%_$(DRIVER_NAME)_bin : DRIVER = $(DRIVER_DIR)/forward2ndOrder/forwardDriver.cpp -I$(CODI_DIR)/include -I$(CODI_DIR)/source -I$(DRIVER_DIR)/forward2ndOrder
$(BUILD_DIR)/%_$(DRIVER_NAME)_bin : $(TEST_DIR)/%.cpp
	@mkdir -p $(@D)
	$(CXX) $(CXX_FLAGS) $< -o $@ $(DRIVER)
	@$(CXX) $(CXX_FLAGS) $(DRIVER) $< -MM -MP -MT $@ -MF $@.d
DRIVER_TEST_FILES := $(patsubst $(TEST_DIR)/%.cpp,$(RESULT_DIR)/%.test,$(DRIVER_TESTS))
$(DRIVER_TEST_FILES): $(RESULT_DIR)/%.test: $(RESULT_DIR)/%_$(DRIVER_NAME).out

# Driver for RealReverseUnchecked
DRIVER_NAME  := RWS_Unch
DRIVER_TESTS := $(BASIC_TESTS) $(REVERSE_TESTS)
$(BUILD_DIR)/%_$(DRIVER_NAME)_bin : DRIVER := $(DRIVER_DIR)/reverseSimple/reverseDriver.cpp -I$(CODI_DIR)/include -I$(DRIVER_DIR)/reverseSimple
$(BUILD_DIR)/%_$(DRIVER_NAME)_bin : $(TEST_DIR)/%.cpp
	@mkdir -p $(@D)
	$(CXX) $(CXX_FLAGS) $< -o $@ $(DRIVER)
	@$(CXX) $(CXX_FLAGS) $(DRIVER) $< -MM -MP -MT $@ -MF $@.d
DRIVER_TEST_FILES := $(patsubst $(TEST_DIR)/%.cpp,$(RESULT_DIR)/%.test,$(DRIVER_TESTS))
$(DRIVER_TEST_FILES): $(RESULT_DIR)/%.test: $(RESULT_DIR)/%_$(DRIVER_NAME).out

# Driver for RealReverseIndexUncheckedIndex
DRIVER_NAME  := RWS_UnchInd
DRIVER_TESTS := $(BASIC_TESTS) $(REVERSE_TESTS)
$(BUILD_DIR)/%_$(DRIVER_NAME)_bin : DRIVER := $(DRIVER_DIR)/reverseSimpleIndex/reverseDriver.cpp -I$(CODI_DIR)/include -I$(DRIVER_DIR)/reverseSimpleIndex
$(BUILD_DIR)/%_$(DRIVER_NAME)_bin : $(TEST_DIR)/%.cpp
	@mkdir -p $(@D)
	$(CXX) $(CXX_FLAGS) $< -o $@ $(DRIVER)
	@$(CXX) $(CXX_FLAGS) $(DRIVER) $< -MM -MP -MT $@ -MF $@.d
DRIVER_TEST_FILES := $(patsubst $(TEST_DIR)/%.cpp,$(RESULT_DIR)/%.test,$(DRIVER_TESTS))
$(DRIVER_TEST_FILES): $(RESULT_DIR)/%.test: $(RESULT_DIR)/%_$(DRIVER_NAME).out

# Driver for RealReversePrimalUnchecked
DRIVER_NAME  := RWS_PrimUnch
DRIVER_TESTS := $(BASIC_TESTS) $(REVERSE_TESTS)
$(BUILD_DIR)/%_$(DRIVER_NAME)_bin : DRIVER := $(DRIVER_DIR)/reversePrimalSimple/reverseDriver.cpp -I$(CODI_DIR)/include -I$(DRIVER_DIR)/reversePrimalSimple
$(BUILD_DIR)/%_$(DRIVER_NAME)_bin : $(TEST_DIR)/%.cpp
	@mkdir -p $(@D)
	$(CXX) $(CXX_FLAGS) $< -o $@ $(DRIVER)
	@$(CXX) $(CXX_FLAGS) $(DRIVER) $< -MM -MP -MT $@ -MF $@.d
DRIVER_TEST_FILES := $(patsubst $(TEST_DIR)/%.cpp,$(RESULT_DIR)/%.test,$(DRIVER_TESTS))
$(DRIVER_TEST_FILES): $(RESULT_DIR)/%.test: $(RESULT_DIR)/%_$(DRIVER_NAME).out

# Driver for RealReversePrimal
DRIVER_NAME  := RWS_Prim
DRIVER_TESTS := $(BASIC_TESTS) $(REVERSE_TESTS)
$(BUILD_DIR)/%_$(DRIVER_NAME)_bin : DRIVER := $(DRIVER_DIR)/reversePrimal/reverseDriver.cpp -I$(CODI_DIR)/include -I$(DRIVER_DIR)/reversePrimal
$(BUILD_DIR)/%_$(DRIVER_NAME)_bin : $(TEST_DIR)/%.cpp
	@mkdir -p $(@D)
	$(CXX) $(CXX_FLAGS) $< -o $@ $(DRIVER)
	@$(CXX) $(CXX_FLAGS) $(DRIVER) $< -MM -MP -MT $@ -MF $@.d
DRIVER_TEST_FILES := $(patsubst $(TEST_DIR)/%.cpp,$(RESULT_DIR)/%.test,$(DRIVER_TESTS))
$(DRIVER_TEST_FILES): $(RESULT_DIR)/%.test: $(RESULT_DIR)/%_$(DRIVER_NAME).out

# Driver for RealReversePrimalIndex
DRIVER_NAME  := RWS_PrimIndex
DRIVER_TESTS := $(BASIC_TESTS) $(REVERSE_TESTS)
$(BUILD_DIR)/%_$(DRIVER_NAME)_bin : DRIVER := $(DRIVER_DIR)/reversePrimalIndex/reverseDriver.cpp -I$(CODI_DIR)/include -I$(DRIVER_DIR)/reversePrimalIndex
$(BUILD_DIR)/%_$(DRIVER_NAME)_bin : $(TEST_DIR)/%.cpp
	@mkdir -p $(@D)
	$(CXX) $(CXX_FLAGS) $< -o $@ $(DRIVER)
	@$(CXX) $(CXX_FLAGS) $(DRIVER) $< -MM -MP -MT $@ -MF $@.d
DRIVER_TEST_FILES := $(patsubst $(TEST_DIR)/%.cpp,$(RESULT_DIR)/%.test,$(DRIVER_TESTS))
$(DRIVER_TEST_FILES): $(RESULT_DIR)/%.test: $(RESULT_DIR)/%_$(DRIVER_NAME).out

# Driver for RealReverse
DRIVER_NAME  := RWS_Chunk
DRIVER_TESTS := $(BASIC_TESTS) $(REVERSE_TESTS)
$(BUILD_DIR)/%_$(DRIVER_NAME)_bin : DRIVER = $(DRIVER_DIR)/reverseChunk/reverseDriver.cpp -I$(CODI_DIR)/include -I$(DRIVER_DIR)/reverseChunk
$(BUILD_DIR)/%_$(DRIVER_NAME)_bin : $(TEST_DIR)/%.cpp
	@mkdir -p $(@D)
	$(CXX) $(CXX_FLAGS) $< -o $@ $(DRIVER)
	@$(CXX) $(CXX_FLAGS) $(DRIVER) $< -MM -MP -MT $@ -MF $@.d
DRIVER_TEST_FILES := $(patsubst $(TEST_DIR)/%.cpp,$(RESULT_DIR)/%.test,$(DRIVER_TESTS))
$(DRIVER_TEST_FILES): $(RESULT_DIR)/%.test: $(RESULT_DIR)/%_$(DRIVER_NAME).out

# Driver for RealReverseVector
DRIVER_NAME  := RWS_ChunkVec
DRIVER_TESTS := $(BASIC_TESTS) $(REVERSE_TESTS)
$(BUILD_DIR)/%_$(DRIVER_NAME)_bin : DRIVER = $(DRIVER_DIR)/reverseChunkVector/reverseDriver.cpp -I$(CODI_DIR)/include -I$(DRIVER_DIR)/reverseChunkVector
$(BUILD_DIR)/%_$(DRIVER_NAME)_bin : $(TEST_DIR)/%.cpp
	@mkdir -p $(@D)
	$(CXX) $(CXX_FLAGS) $< -o $@ $(DRIVER)
	@$(CXX) $(CXX_FLAGS) $(DRIVER) $< -MM -MP -MT $@ -MF $@.d
DRIVER_TEST_FILES := $(patsubst $(TEST_DIR)/%.cpp,$(RESULT_DIR)/%.test,$(DRIVER_TESTS))
$(DRIVER_TEST_FILES): $(RESULT_DIR)/%.test: $(RESULT_DIR)/%_$(DRIVER_NAME).out

# Driver for RealReverseIndex
DRIVER_NAME  := RWS_ChunkInd
DRIVER_TESTS := $(BASIC_TESTS) $(REVERSE_TESTS)
$(BUILD_DIR)/%_$(DRIVER_NAME)_bin : DRIVER = $(DRIVER_DIR)/reverseChunkIndex/reverseDriver.cpp -I$(CODI_DIR)/include -I$(DRIVER_DIR)/reverseChunkIndex
$(BUILD_DIR)/%_$(DRIVER_NAME)_bin : $(TEST_DIR)/%.cpp
	@mkdir -p $(@D)
	$(CXX) $(CXX_FLAGS) $< -o $@ $(DRIVER)
	@$(CXX) $(CXX_FLAGS) $(DRIVER) $< -MM -MP -MT $@ -MF $@.d
DRIVER_TEST_FILES := $(patsubst $(TEST_DIR)/%.cpp,$(RESULT_DIR)/%.test,$(DRIVER_TESTS))
$(DRIVER_TEST_FILES): $(RESULT_DIR)/%.test: $(RESULT_DIR)/%_$(DRIVER_NAME).out

# Driver for RealReverseIndexAssignOpt
DRIVER_NAME  := RWS_ChunkIndAssignOpt
DRIVER_TESTS := $(BASIC_TESTS) $(REVERSE_TESTS)
$(BUILD_DIR)/%_$(DRIVER_NAME)_bin : DRIVER = $(DRIVER_DIR)/reverseChunkIndexAssignOpt/reverseDriver.cpp -I$(CODI_DIR)/include -I$(DRIVER_DIR)/reverseChunkIndex
$(BUILD_DIR)/%_$(DRIVER_NAME)_bin : $(TEST_DIR)/%.cpp
	@mkdir -p $(@D)
	$(CXX) $(CXX_FLAGS) $< -o $@ $(DRIVER)
	@$(CXX) $(CXX_FLAGS) $(DRIVER) $< -MM -MP -MT $@ -MF $@.d
DRIVER_TEST_FILES := $(patsubst $(TEST_DIR)/%.cpp,$(RESULT_DIR)/%.test,$(DRIVER_TESTS))
$(DRIVER_TEST_FILES): $(RESULT_DIR)/%.test: $(RESULT_DIR)/%_$(DRIVER_NAME).out

# Driver for RealReverseIndexVector
DRIVER_NAME  := RWS_ChunkIndVec
DRIVER_TESTS := $(BASIC_TESTS) $(REVERSE_TESTS)
$(BUILD_DIR)/%_$(DRIVER_NAME)_bin : DRIVER = $(DRIVER_DIR)/reverseChunkIndexVector/reverseDriver.cpp -I$(CODI_DIR)/include -I$(DRIVER_DIR)/reverseChunkIndexVector
$(BUILD_DIR)/%_$(DRIVER_NAME)_bin : $(TEST_DIR)/%.cpp
	@mkdir -p $(@D)
	$(CXX) $(CXX_FLAGS) $< -o $@ $(DRIVER)
	@$(CXX) $(CXX_FLAGS) $(DRIVER) $< -MM -MP -MT $@ -MF $@.d
DRIVER_TEST_FILES := $(patsubst $(TEST_DIR)/%.cpp,$(RESULT_DIR)/%.test,$(DRIVER_TESTS))
$(DRIVER_TEST_FILES): $(RESULT_DIR)/%.test: $(RESULT_DIR)/%_$(DRIVER_NAME).out

# Driver for reverse over forward
DRIVER_NAME  := RWS2nd
DRIVER_TESTS := $(BASIC_TESTS) $(REVERSE_TESTS)
$(BUILD_DIR)/%_$(DRIVER_NAME)_bin : DRIVER = $(DRIVER_DIR)/reverseOverForward/reverseOverForwardDriver.cpp -I$(CODI_DIR)/include -I$(DRIVER_DIR)/reverseOverForward
$(BUILD_DIR)/%_$(DRIVER_NAME)_bin : $(TEST_DIR)/%.cpp
	@mkdir -p $(@D)
	$(CXX) $(CXX_FLAGS) $< -o $@ $(DRIVER)
	@$(CXX) $(CXX_FLAGS) $(DRIVER) $< -MM -MP -MT $@ -MF $@.d
DRIVER_TEST_FILES := $(patsubst $(TEST_DIR)/%.cpp,$(RESULT_DIR)/%.test,$(DRIVER_TESTS))
$(DRIVER_TEST_FILES): $(RESULT_DIR)/%.test: $(RESULT_DIR)/%_$(DRIVER_NAME).out

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
