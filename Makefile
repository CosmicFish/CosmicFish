#----------------------------------------------------------------------------------------
#
# This file is part of CosmicFish.
#
# Copyright (C) 2015-2016 by the CosmicFish authors
#
# The CosmicFish code is free software;
# You can use it, redistribute it, and/or modify it under the terms
# of the GNU General Public License as published by the Free Software Foundation;
# either version 3 of the License, or (at your option) any later version.
# The full text of the license can be found in the file LICENSE at
# the top level of the CosmicFish distribution.
#
#----------------------------------------------------------------------------------------

#
# This file contains the main part of the CosmicFish build system.
# Additional build systema macros can be found in the folder make_macros.
#

################### directory structure ########################################

# fix the paths of the project:
PROJECT_DIR=$(abspath .)
SOURCE_DIR=$(PROJECT_DIR)/sources
BUILD_DIR=$(PROJECT_DIR)/build
APPS_DIR=$(PROJECT_DIR)/apps
DOCUMENTATION_DIR=$(PROJECT_DIR)/documentation
EXAMPLES_DIR=$(PROJECT_DIR)/examples
TEST_DIR=$(PROJECT_DIR)/test
BENCHMARK_DIR=$(PROJECT_DIR)/benchmark
BUNDLED_DIR=$(PROJECT_DIR)/bundled
SCRIPT_DIR=$(PROJECT_DIR)/script
MAKE_DIR=$(PROJECT_DIR)/make_macros

# include compiler and options:
include $(MAKE_DIR)/compilers_options.make
# include bundled software:
include $(MAKE_DIR)/bundled_software.make
# include documentation targets:
include $(MAKE_DIR)/documentation.make

# path to the boltzman solvers:
CAMB_PATH=$(PROJECT_DIR)/camb/camb/
EFTCAMB_PATH=$(PROJECT_DIR)/camb/eftcamb/
MGCAMB_PATH=$(PROJECT_DIR)/camb/mgcamb/

# directory structure of camb:
CAMB_OBJ_DIR := $(BUILD_DIR)/fisher_camb/build
CAMB_INC_DIR := $(BUILD_DIR)/fisher_camb/include
CAMB_BIN_DIR := $(BUILD_DIR)/fisher_camb/bin
CAMB_LIB_DIR := $(BUILD_DIR)/fisher_camb/lib
# directory structure of eftcamb:
EFTCAMB_OBJ_DIR := $(BUILD_DIR)/fisher_eftcamb/build
EFTCAMB_INC_DIR := $(BUILD_DIR)/fisher_eftcamb/include
EFTCAMB_BIN_DIR := $(BUILD_DIR)/fisher_eftcamb/bin
EFTCAMB_LIB_DIR := $(BUILD_DIR)/fisher_eftcamb/lib
# directory structure of mgcamb:
MGCAMB_OBJ_DIR := $(BUILD_DIR)/fisher_mgcamb/build
MGCAMB_INC_DIR := $(BUILD_DIR)/fisher_mgcamb/include
MGCAMB_BIN_DIR := $(BUILD_DIR)/fisher_mgcamb/bin
MGCAMB_LIB_DIR := $(BUILD_DIR)/fisher_mgcamb/lib

# linking instructions for camb:
CAMB_F90_LINK=-L $(CAMB_LIB_DIR) -lcamb
CAMB_F90_APP_LINK=-L $(CAMB_LIB_DIR) -lfisher_camb -lcamb
CAMB_F90_INCLUDE=-I $(CAMB_INC_DIR)
# linking instructions for eftcamb:
EFTCAMB_F90_LINK=-L $(EFTCAMB_LIB_DIR) -lcamb
EFTCAMB_F90_APP_LINK=-L $(EFTCAMB_LIB_DIR) -lfisher_eftcamb -lcamb
EFTCAMB_F90_INCLUDE=-I $(EFTCAMB_INC_DIR)
# linking instructions for mgcamb:
MGCAMB_F90_LINK=-L $(MGCAMB_LIB_DIR) -lcamb
MGCAMB_F90_APP_LINK=-L $(MGCAMB_LIB_DIR) -lfisher_mgcamb -lcamb
MGCAMB_F90_INCLUDE=-I $(MGCAMB_INC_DIR)

# flags for camb:
CAMB_FFLAGS = $(FFLAGS) -DCOSMICFISH_CAMB
# flags for eftcamb:
EFTCAMB_FFLAGS = $(FFLAGS) -DCOSMICFISH_EFTCAMB
# flags for mgcamb:
MGCAMB_FFLAGS = $(FFLAGS) -DCOSMICFISH_MGCAMB

################### high level targets #########################################

# cosmicfish target:
cosmicfish: apps test

# default target:
default: cosmicfish

# all target:
all: apps test run_test documentation examples

################### Boltzman solvers targets ###################################

# ISSUES AND TODOS: now the target does not see wether a camb file was changed or not.

# target to create the Boltzman solvers libraries:
bsolvers: camb eftcamb mgcamb

camb: directories_camb
	@cd $(CAMB_PATH) && $(MAKE) all CAMBLIB=$(CAMB_LIB_DIR)/libcamb.a -e OUTPUT_DIR=cosmicfish_build
	@cp $(CAMB_PATH)cosmicfish_build/*.mod $(CAMB_INC_DIR)
	@cp $(CAMB_PATH)cosmicfish_build/*.o $(CAMB_OBJ_DIR)

eftcamb: directories_eftcamb
	@cd $(EFTCAMB_PATH) && $(MAKE) all CAMBLIB=$(EFTCAMB_LIB_DIR)/libcamb.a -e OUTPUT_DIR=cosmicfish_build
	@cp $(EFTCAMB_PATH)cosmicfish_build/*.mod $(EFTCAMB_INC_DIR)
	@cp $(EFTCAMB_PATH)cosmicfish_build/*.o $(EFTCAMB_OBJ_DIR)

mgcamb: directories_mgcamb
	cd $(MGCAMB_PATH) && $(MAKE) all CAMBLIB=$(MGCAMB_LIB_DIR)/libcamb.a -e OUTPUT_DIR=cosmicfish_build
	@cp $(MGCAMB_PATH)cosmicfish_build/*.mod $(MGCAMB_INC_DIR)
	@cp $(MGCAMB_PATH)cosmicfish_build/*.o $(MGCAMB_OBJ_DIR)

################### Fisher libraries targets ###################################

# Fisher libraries source files:
FISHER_SOURCES := $(wildcard $(SOURCE_DIR)/*.f90)
FISHER_TEMP_1 := $(patsubst %.f90,%,$(FISHER_SOURCES))
FISHER_TEMP_2 := $(notdir $(FISHER_TEMP_1))
FISHER_FILES := $(addsuffix .o, $(FISHER_TEMP_2))

# fisher files for camb:
CAMB_FISHER_FILES = $(addprefix $(CAMB_OBJ_DIR)/, $(FISHER_FILES))
# fisher files for eftcamb:
EFTCAMB_FISHER_FILES = $(addprefix $(EFTCAMB_OBJ_DIR)/, $(FISHER_FILES))
# fisher files for mgcamb:
MGCAMB_FISHER_FILES = $(addprefix $(MGCAMB_OBJ_DIR)/, $(FISHER_FILES))

# build the depedency tree for the fisher library:
MAKEDEPEND := $(BUNDLED_DIR)/fort_depend/fort_depend.py

CAMB_DEP_FILE    := $(SOURCE_DIR)/camb_dependency_tree.dep
EFTCAMB_DEP_FILE := $(SOURCE_DIR)/eftcamb_dependency_tree.dep
MGCAMB_DEP_FILE  := $(SOURCE_DIR)/mgcamb_dependency_tree.dep

DEP_OBJECTS := $(notdir $(FISHER_SOURCES))

# Make dependencies target:

# The .dep file depends on the source files, so it automatically gets updated
# when you change your source
$(CAMB_DEP_FILE): $(FISHER_SOURCES)
	@echo "Making cosmicfish camb dependencies"
	@cd $(SOURCE_DIR) && $(MAKEDEPEND) -w -o $(CAMB_DEP_FILE) -f $(DEP_OBJECTS) -p '$$(CAMB_OBJ_DIR)'/ > /dev/null 2>&1

$(EFTCAMB_DEP_FILE): $(FISHER_SOURCES)
	@echo "Making cosmicfish eftcamb dependencies"
	@cd $(SOURCE_DIR) && $(MAKEDEPEND) -w -o $(EFTCAMB_DEP_FILE) -f $(DEP_OBJECTS) -p '$$(EFTCAMB_OBJ_DIR)'/ > /dev/null 2>&1

$(MGCAMB_DEP_FILE): $(FISHER_SOURCES)
	@echo "Making cosmicfish mgcamb dependencies"
	@cd $(SOURCE_DIR) && $(MAKEDEPEND) -w -o $(MGCAMB_DEP_FILE) -f $(DEP_OBJECTS) -p '$$(MGCAMB_OBJ_DIR)'/ > /dev/null 2>&1

include $(CAMB_DEP_FILE)
include $(EFTCAMB_DEP_FILE)
include $(MGCAMB_DEP_FILE)

# the targets:
fisher_lib: fisher_camb fisher_eftcamb fisher_mgcamb

fisher_camb: camb $(CAMB_DEP_FILE) $(CAMB_FISHER_FILES)
	ar -r $(CAMB_LIB_DIR)/lib$@.a $(CAMB_FISHER_FILES)

fisher_eftcamb: eftcamb $(EFTCAMB_DEP_FILE) $(EFTCAMB_FISHER_FILES)
	ar -r $(EFTCAMB_LIB_DIR)/lib$@.a $(EFTCAMB_FISHER_FILES)

fisher_mgcamb: mgcamb $(MGCAMB_DEP_FILE) $(MGCAMB_FISHER_FILES)
	ar -r $(MGCAMB_LIB_DIR)/lib$@.a $(MGCAMB_FISHER_FILES)

# target for the object files in the fisher library:

$(CAMB_OBJ_DIR)/%.o: $(SOURCE_DIR)/%.f90 camb
	$(F90C) $(CAMB_FFLAGS) $(CAMB_F90_LINK) $(LAPACK_LINK) $(CAMB_F90_INCLUDE) $(LAPACK_INCLUDE) $(MODULE_FLAG)$(CAMB_INC_DIR) -c $(SOURCE_DIR)/$*.f90 -o $(CAMB_OBJ_DIR)/$*.o

$(EFTCAMB_OBJ_DIR)/%.o: $(SOURCE_DIR)/%.f90 eftcamb
	$(F90C) $(EFTCAMB_FFLAGS) $(EFTCAMB_F90_LINK) $(LAPACK_LINK) $(EFTCAMB_F90_INCLUDE) $(LAPACK_INCLUDE) $(MODULE_FLAG)$(EFTCAMB_INC_DIR) -c $(SOURCE_DIR)/$*.f90 -o $(EFTCAMB_OBJ_DIR)/$*.o

$(MGCAMB_OBJ_DIR)/%.o: $(SOURCE_DIR)/%.f90 mgcamb
	$(F90C) $(MGCAMB_FFLAGS) $(MGCAMB_F90_LINK) $(LAPACK_LINK) $(MGCAMB_F90_INCLUDE) $(LAPACK_INCLUDE) $(MODULE_FLAG)$(MGCAMB_INC_DIR) -c $(SOURCE_DIR)/$*.f90 -o $(MGCAMB_OBJ_DIR)/$*.o

################################ Apps targets ##################################
# target to create the cosmicfish executables

APPS_SOURCES := $(wildcard $(APPS_DIR)/*.f90)
APPS_NAMES   := $(notdir $(patsubst %.f90,%,$(APPS_SOURCES)))
APPS_BINS    := $(addsuffix .x, $(APPS_NAMES))

CAMB_BINS    := $(addprefix $(CAMB_BIN_DIR)/, $(APPS_BINS))
EFTCAMB_BINS := $(addprefix $(EFTCAMB_BIN_DIR)/, $(APPS_BINS))
MGCAMB_BINS  := $(addprefix $(MGCAMB_BIN_DIR)/, $(APPS_BINS))

.PHONY: apps

apps: camb_apps eftcamb_apps mgcamb_apps

camb_apps: fisher_camb $(CAMB_BINS) 
eftcamb_apps: fisher_eftcamb $(EFTCAMB_BINS) 
mgcamb_apps: fisher_mgcamb $(MGCAMB_BINS) 

$(CAMB_BIN_DIR)/%.x: $(APPS_DIR)/%.f90 $(CAMB_FISHER_FILES) fisher_camb
	$(F90C) $(FFLAGS) $(APPS_DIR)/$*.f90 $(CAMB_F90_APP_LINK) $(LAPACK_LINK) $(CAMB_F90_INCLUDE) $(LAPACK_INCLUDE) -o $(CAMB_BIN_DIR)/$*.x

$(EFTCAMB_BIN_DIR)/%.x: $(APPS_DIR)/%.f90 $(EFTCAMB_FISHER_FILES) fisher_eftcamb
	$(F90C) $(FFLAGS) $(APPS_DIR)/$*.f90 $(EFTCAMB_F90_APP_LINK) $(LAPACK_LINK) $(EFTCAMB_F90_INCLUDE) $(LAPACK_INCLUDE) -o $(EFTCAMB_BIN_DIR)/$*.x

$(MGCAMB_BIN_DIR)/%.x: $(APPS_DIR)/%.f90 $(MGCAMB_FISHER_FILES) fisher_mgcamb
	$(F90C) $(FFLAGS) $(APPS_DIR)/$*.f90 $(MGCAMB_F90_APP_LINK) $(LAPACK_LINK) $(MGCAMB_F90_INCLUDE) $(LAPACK_INCLUDE) -o $(MGCAMB_BIN_DIR)/$*.x

clean_apps:
	@rm -rf $(CAMB_BIN_DIR)/*.x $(EFTCAMB_BIN_DIR)/*.x $(MGCAMB_BIN_DIR)/*.x

################################ Test targets ##################################
# This section takes care of compiling and cleaning the cosmicfish test suite

TEST_SOURCES := $(wildcard $(TEST_DIR)/test_source/*.f90)
TEST_NAMES   := $(notdir $(patsubst %.f90,%,$(TEST_SOURCES)))
TEST_BINS    := $(addsuffix .x, $(TEST_NAMES))

CAMB_TEST    := $(addprefix $(TEST_DIR)/test_build/camb_, $(TEST_BINS))
EFTCAMB_TEST := $(addprefix $(TEST_DIR)/test_build/eftcamb_, $(TEST_BINS))
MGCAMB_TEST  := $(addprefix $(TEST_DIR)/test_build/mgcamb_, $(TEST_BINS))

.PHONY: test

test: camb_test eftcamb_test mgcamb_test python_test

camb_test: fisher_camb $(CAMB_TEST)
eftcamb_test: fisher_eftcamb $(EFTCAMB_TEST)
mgcamb_test: fisher_mgcamb $(MGCAMB_TEST)

python_test:
	@rm -rf .coverage
	@nosetests -v -s $(PROJECT_DIR)/python/cosmicfish_pylib_test --with-coverage --cover-package=cosmicfish_pylib
	
$(TEST_DIR)/test_build/camb_%.x: $(TEST_DIR)/test_source/%.f90  $(CAMB_FISHER_FILES) fisher_camb
	$(F90C) $(CAMB_FFLAGS) $(TEST_DIR)/test_source/$*.f90 $(CAMB_F90_APP_LINK) $(LAPACK_LINK) $(CAMB_F90_INCLUDE) $(LAPACK_INCLUDE) -o $(TEST_DIR)/test_build/camb_$*.x

$(TEST_DIR)/test_build/eftcamb_%.x: $(TEST_DIR)/test_source/%.f90 $(EFTCAMB_FISHER_FILES) fisher_eftcamb
	$(F90C) $(EFTCAMB_FFLAGS) $(TEST_DIR)/test_source/$*.f90 $(EFTCAMB_F90_APP_LINK) $(LAPACK_LINK) $(EFTCAMB_F90_INCLUDE) $(LAPACK_INCLUDE) -o $(TEST_DIR)/test_build/eftcamb_$*.x

$(TEST_DIR)/test_build/mgcamb_%.x: $(TEST_DIR)/test_source/%.f90 $(MGCAMB_FISHER_FILES) fisher_mgcamb
	$(F90C) $(MGCAMB_FFLAGS) $(TEST_DIR)/test_source/$*.f90 $(MGCAMB_F90_APP_LINK) $(LAPACK_LINK) $(MGCAMB_F90_INCLUDE) $(LAPACK_INCLUDE) -o $(TEST_DIR)/test_build/mgcamb_$*.x

run_test: test numdiff
	@cd $(TEST_DIR) && make test

# clean the test build (do not remove fail logs):
clean_test:
	@rm -rf $(TEST_DIR)/test_build/*.x

############################# Benchmark targets ################################


############################# Example targets ################################

.PHONY: examples

examples: apps
	@cd $(EXAMPLES_DIR) && make all
	
clean_examples:	
	@cd $(EXAMPLES_DIR) && make deep_clean

############################## Directory targets ###############################
# This section takes care of creating the directory structure of the cosmicfish
# build.

# create the build directory and its structure:
directories: directories_buid directories_camb directories_eftcamb directories_mgcamb

# create the build directory:
directories_buid:
	@mkdir -p $(BUILD_DIR)
# create directories for cosmicfish built with camb:
directories_camb: directories_buid
	@mkdir -p $(BUILD_DIR)/fisher_camb/lib $(BUILD_DIR)/fisher_camb/include $(BUILD_DIR)/fisher_camb/build $(BUILD_DIR)/fisher_camb/bin
# create directories for cosmicfish built with eftcamb:
directories_eftcamb: directories_buid
	@mkdir -p $(BUILD_DIR)/fisher_eftcamb/lib $(BUILD_DIR)/fisher_eftcamb/include $(BUILD_DIR)/fisher_eftcamb/build $(BUILD_DIR)/fisher_eftcamb/bin
# create directories for cosmicfish built with mgcamb:
directories_mgcamb: directories_buid
	@mkdir -p $(BUILD_DIR)/fisher_mgcamb/lib $(BUILD_DIR)/fisher_mgcamb/include $(BUILD_DIR)/fisher_mgcamb/build $(BUILD_DIR)/fisher_mgcamb/bin

############################## Clean targets ###################################
# This section contains the targets to clean the build of the cosmicfish library

# target to clean the Boltzman solvers:
clean_bsolvers: clean_camb clean_eftcamb clean_mgcamb
# clean camb:
clean_camb:
	@cd $(CAMB_PATH) && make clean OUTPUT_DIR=cosmicfish_build && rm -rf camb
# clean eftcamb:
clean_eftcamb:
	@cd $(EFTCAMB_PATH) && make clean OUTPUT_DIR=cosmicfish_build && rm -rf camb
# clean mgcamb:
clean_mgcamb:
	@cd $(MGCAMB_PATH) && make clean OUTPUT_DIR=cosmicfish_build && rm -rf camb
# tidy up the build directory (cosmicfish library will still be running):
clean: clean_bsolvers
	@rm -rf $(CAMB_OBJ_DIR) $(EFTCAMB_OBJ_DIR) $(MGCAMB_OBJ_DIR)
# remove dependency files:
clean_dep:
	@rm -rf $(CAMB_DEP_FILE) $(EFTCAMB_DEP_FILE) $(MGCAMB_DEP_FILE)
# totally delete everything:
distclean: clean_bsolvers clean_documentation clean_test clean_bundled
	@rm -rf $(BUILD_DIR)
