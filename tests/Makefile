# names of the basic deriectories
TEST_DIR = tests
BUILD_DIR = build
DRIVER_DIR = drivers
RESULT_BASE_DIR = results_base
RESULT_DIR = results

#list all source files in TEST_DIR
TEST_FILES  = $(wildcard $(TEST_DIR)/Test**.cpp)
TEST_FILES  += $(wildcard $(TEST_DIR)/**/Test**.cpp)
TEST_FILES  += $(wildcard $(TEST_DIR)/**/**/Test**.cpp)

ifeq ($(OPT), yes)
  CXX_FLAGS := -O3 -Wall -pedantic -std=c++11
else
  CXX_FLAGS := -O0 -g -Wall -pedantic -std=c++11
endif

CXX := g++
ADEPT_DIR := ../../adept-extended
CODI_DIR := ..

INC =  -I$(ADEPT_DIR)/include

TESTS = $(patsubst $(TEST_DIR)/%.cpp,$(RESULT_DIR)/%.test,$(TEST_FILES))

FORWARD_BIN_FILES   = $(patsubst $(TEST_DIR)/%.cpp,$(BUILD_DIR)/%_forward_bin,$(TEST_FILES))
FORWARD_OUT_FILES   = $(patsubst $(TEST_DIR)/%.cpp,$(RESULT_DIR)/%_forward.out,$(TEST_FILES))
FORWARD_2ND_BIN_FILES   = $(patsubst $(TEST_DIR)/%.cpp,$(BUILD_DIR)/%_forward2nd_bin,$(TEST_FILES))
FORWARD_2ND_OUT_FILES   = $(patsubst $(TEST_DIR)/%.cpp,$(RESULT_DIR)/%_forward2nd.out,$(TEST_FILES))
REVERSE_BIN_FILES   = $(patsubst $(TEST_DIR)/%.cpp,$(BUILD_DIR)/%_reverse_bin,$(TEST_FILES))
REVERSE_OUT_FILES   = $(patsubst $(TEST_DIR)/%.cpp,$(RESULT_DIR)/%_reverse.out,$(TEST_FILES))

FORWARD_DRIVER = $(DRIVER_DIR)/forward/forwardDriver.cpp -I$(CODI_DIR)/include -I$(CODI_DIR)/source -I$(DRIVER_DIR)/forward
FORWARD_2ND_DRIVER = $(DRIVER_DIR)/forward2ndOrder/forwardDriver.cpp -I$(CODI_DIR)/include -I$(CODI_DIR)/source -I$(DRIVER_DIR)/forward2ndOrder
REVERSE_DRIVER = $(DRIVER_DIR)/reverse/reverseDriver.cpp -I$(CODI_DIR)/include -I$(DRIVER_DIR)/reverse

.SECONDARY: $(FORWARD_OUT_FILES) $(FORWARD_BIN_FILES) $(FORWARD_2ND_OUT_FILES) $(FORWARD_2ND_BIN_FILES) $(REVERSE_OUT_FILES) $(REVERSE_BIN_FILES)

# rules for test creation
$(BUILD_DIR)/%_forward_bin : $(TEST_DIR)/%.cpp
	@mkdir -p $(@D)
	$(CXX) $(CXX_FLAGS) $< -o $@ $(FORWARD_DRIVER)

$(BUILD_DIR)/%_forward2nd_bin : $(TEST_DIR)/%.cpp
	@mkdir -p $(@D)
	$(CXX) $(CXX_FLAGS) $< -o $@ $(FORWARD_DRIVER)

$(BUILD_DIR)/%_reverse_bin : $(TEST_DIR)/%.cpp
	@mkdir -p $(@D)
	$(CXX) $(CXX_FLAGS) $< -o $@ $(REVERSE_DRIVER)

# rules for generating the text files
$(RESULT_DIR)/%.out : $(BUILD_DIR)/%_bin
	@mkdir -p $(@D)
	$< > $@

# rule for printing the results
$(RESULT_DIR)/%.test : $(RESULT_DIR)/%_forward.out $(RESULT_DIR)/%_forward2nd.out $(RESULT_DIR)/%_reverse.out
	@./compare.sh -n $* -b $(RESULT_BASE_DIR)/$*.out $(RESULT_DIR)/$*_forward.out $(RESULT_DIR)/$*_forward2nd.out $(RESULT_DIR)/$*_reverse.out

all: $(TESTS)
	@mkdir -p $(BUILD_DIR)

Test% : $(RESULT_DIR)/Test%.test
	@echo Test$* build

.PHONY: clean
clean:
	rm -fr $(BUILD_DIR)
	rm -fr $(RESULT_DIR)
