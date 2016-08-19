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
# This file contains part of the CosmicFish build system: documentation
# This section takes care of compiling the documentation

############################### global targets #################################

.PHONY: documentation

documentation: doxydoc pydoc

clean_documentation: clean_doxydoc clean_pydoc

############################### doxygen documentation targets ################################

doxydoc: doxygen
ifdef DOXYGEN_COMMAND
	@cd $(DOCUMENTATION_DIR)/doxygen && $(DOXYGEN_COMMAND) CosmicFish_doxyfile
else
	@cd $(DOCUMENTATION_DIR)/doxygen && $(DOXY_APP) CosmicFish_doxyfile
endif

clean_doxydoc:
	@rm -rf $(DOCUMENTATION_DIR)/doxygen/doxygen_doc/*

############################### python documentation targets ################################

ipydoc:
	@ipython nbconvert --ExecutePreprocessor.enabled=True --output $(DOCUMENTATION_DIR)/python/example --to html --template full index.ipynb && rm test.png

pydoc: ipydoc
	@cd $(DOCUMENTATION_DIR)/python && make html SPHINXOPTS=-v

clean_pydoc:
	@rm -rf $(DOCUMENTATION_DIR)/python/_build
	@rm -rf $(DOCUMENTATION_DIR)/python/_static
	@rm -rf $(DOCUMENTATION_DIR)/python/_summaries
	@rm -rf $(DOCUMENTATION_DIR)/python/_templates
