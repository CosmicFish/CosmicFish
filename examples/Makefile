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
# This file contains the make targets for the examples
#

# structure of the folder:

EXAMPLES_DIR=$(abspath .)
SCRIPT_DIR=$(EXAMPLES_DIR)/script
RAW_RESULTS_DIR=$(EXAMPLES_DIR)/raw_results
RESULTS_DIR=$(EXAMPLES_DIR)/results
PARAMETER_DIR=$(EXAMPLES_DIR)/parameters
PARAMETER_ANALYSIS_DIR=$(EXAMPLES_DIR)/parameters_analysis

# make general targets:

default: all

all: targets fisher_matrices analysis

# Full analysis target:
	
additional_script: targets
	@for i in $(SCRIPT_DIR)/*.py ; do python $$i ; done;

analysis: targets
	@bash $(SCRIPT_DIR)/do_analysis.sh 

fisher_matrices: targets
	@bash $(SCRIPT_DIR)/do_fisher.sh 

# Target to create a makefile that has all single targets:

targets:
	@bash $(SCRIPT_DIR)/create_target.sh 

include $(SCRIPT_DIR)/makefile_target.make

# clean targets:

clean:
	@rm -rf $(RAW_RESULTS_DIR)/*.txt
	@rm -rf $(RAW_RESULTS_DIR)/*.pdf
	@rm -rf $(RAW_RESULTS_DIR)/*.dat
	@rm -rf $(RAW_RESULTS_DIR)/*.paramnames
	@rm -rf $(RAW_RESULTS_DIR)/*.ini

deep_clean: clean
	@rm -rf $(RESULTS_DIR)/*.txt
	@rm -rf $(RESULTS_DIR)/*.pdf
	@rm -rf $(RESULTS_DIR)/*.dat
	@rm -rf $(RESULTS_DIR)/*.paramnames
	@rm -rf $(RESULTS_DIR)/*.ini




	

