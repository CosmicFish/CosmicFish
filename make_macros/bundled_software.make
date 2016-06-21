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
# This file contains part of the CosmicFish build system: bundled software installer and software dependencies.
# The CosmicFish library needs several external packages, mainly to perform utility duties. This part of the Makefile 
# tries to install them. It might fail because the installation is very naive. If this happens the user has to 
# install manually the packages.
#

############################### global targets #################################

.PHONY: bundled

bundled: numdiff doxygen

clean_bundled: clean_numdiff clean_doxygen

############################### numdiff targets ################################
# test if the system has the numdiff command
NUMDIFF_COMMAND := $(shell which numdiff 2>/dev/null)

# main numdiff target:
.PHONY: numdiff

numdiff: $(BUNDLED_DIR)/numdiff/numdiff

# build the numdiff target if the system one is not found:
$(BUNDLED_DIR)/numdiff/numdiff:
ifdef NUMDIFF_COMMAND
else
	@cd $(BUNDLED_DIR)/numdiff && ./configure && touch ./* && make
endif

# forces building numdiff just to see if everything goes fine...
force_numdiff:
	@cd $(BUNDLED_DIR)/numdiff && ./configure && touch ./* && make

clean_numdiff:
	@rm -rf $(BUNDLED_DIR)/numdiff/config.h
	@rm -rf $(BUNDLED_DIR)/numdiff/config.log
	@rm -rf $(BUNDLED_DIR)/numdiff/config.status
	@rm -rf $(BUNDLED_DIR)/numdiff/Makefile
	@rm -rf $(BUNDLED_DIR)/numdiff/numdiff
	@rm -rf $(BUNDLED_DIR)/numdiff/ndselect
	
############################### doxygen targets ################################
# test if the system has the doxygen command
DOXYGEN_COMMAND := $(shell which doxygen 2>/dev/null)
DOXY_APP := $(BUNDLED_DIR)/doxygen/bin/doxygen

# main doxygen target:
doxygen: $(BUNDLED_DIR)/doxygen/bin/doxygen

# build the doxygen target if needed:
$(BUNDLED_DIR)/doxygen/bin/doxygen:
ifdef DOXYGEN_COMMAND
else
	@cd $(BUNDLED_DIR)/doxygen && ./configure && touch ./* && make
endif

# forces building doxygen just to see if everything goes fine...
force_doxygen:
	@cd $(BUNDLED_DIR)/doxygen && ./configure && touch ./* && make

clean_doxygen:
	@rm -rf $(BUNDLED_DIR)/doxygen/.makeconfig
	@rm -rf $(BUNDLED_DIR)/doxygen/.tmakeconfig
	@rm -rf $(BUNDLED_DIR)/doxygen/Makefile
	@rm -rf $(BUNDLED_DIR)/doxygen/VERSION
	@rm -rf $(BUNDLED_DIR)/doxygen/addon/doxmlparser/examples/metrics/Makefile
	@rm -rf $(BUNDLED_DIR)/doxygen/addon/doxmlparser/examples/metrics/metrics.pro
	@rm -rf $(BUNDLED_DIR)/doxygen/addon/doxmlparser/src/Makefile
	@rm -rf $(BUNDLED_DIR)/doxygen/addon/doxmlparser/src/doxmlparser.pro
	@rm -rf $(BUNDLED_DIR)/doxygen/addon/doxmlparser/test/Makefile
	@rm -rf $(BUNDLED_DIR)/doxygen/addon/doxmlparser/test/xmlparse.pro
	@rm -rf $(BUNDLED_DIR)/doxygen/addon/doxyapp/Makefile
	@rm -rf $(BUNDLED_DIR)/doxygen/addon/doxyapp/doxyapp.pro
	@rm -rf $(BUNDLED_DIR)/doxygen/addon/doxysearch/Makefile
	@rm -rf $(BUNDLED_DIR)/doxygen/addon/doxysearch/doxyindexer.pro
	@rm -rf $(BUNDLED_DIR)/doxygen/addon/doxysearch/doxysearch.pro
	@rm -rf $(BUNDLED_DIR)/doxygen/addon/doxywizard/Makefile
	@rm -rf $(BUNDLED_DIR)/doxygen/addon/doxywizard/doxywizard.pro
	@rm -rf $(BUNDLED_DIR)/doxygen/doc/Makefile
	@rm -rf $(BUNDLED_DIR)/doxygen/examples/Makefile
	@rm -rf $(BUNDLED_DIR)/doxygen/generated_src/
	@rm -rf $(BUNDLED_DIR)/doxygen/libmd5/Makefile
	@rm -rf $(BUNDLED_DIR)/doxygen/libmd5/libmd5.pro
	@rm -rf $(BUNDLED_DIR)/doxygen/qtools/Makefile
	@rm -rf $(BUNDLED_DIR)/doxygen/qtools/qtools.pro
	@rm -rf $(BUNDLED_DIR)/doxygen/src/Makefile
	@rm -rf $(BUNDLED_DIR)/doxygen/src/doxygen.pro
	@rm -rf $(BUNDLED_DIR)/doxygen/src/libdoxycfg.pro
	@rm -rf $(BUNDLED_DIR)/doxygen/src/libdoxycfg.t
	@rm -rf $(BUNDLED_DIR)/doxygen/src/libdoxygen.pro
	@rm -rf $(BUNDLED_DIR)/doxygen/src/libdoxygen.t
	@rm -rf $(BUNDLED_DIR)/doxygen/vhdlparser/Makefile
	@rm -rf $(BUNDLED_DIR)/doxygen/vhdlparser/vhdlparser.pro
	@rm -rf $(BUNDLED_DIR)/doxygen/bin/
	@rm -rf $(BUNDLED_DIR)/doxygen/libmd5/Makefile.libmd5
	@rm -rf $(BUNDLED_DIR)/doxygen/qtools/Makefile.qtools
	@rm -rf $(BUNDLED_DIR)/doxygen/src/Makefile.doxygen
	@rm -rf $(BUNDLED_DIR)/doxygen/src/Makefile.libdoxycfg
	@rm -rf $(BUNDLED_DIR)/doxygen/src/Makefile.libdoxygen
	@rm -rf $(BUNDLED_DIR)/doxygen/vhdlparser/Makefile.vhdlparser