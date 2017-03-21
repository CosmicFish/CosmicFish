#!/bin/bash

#----------------------------------------------------------------------------------------
#
# This file is part of CosmicFish.
#
# Copyright (C) 2015-2017 by the CosmicFish authors
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
# This file contains some definitions that are common to all examples script
#

# get the path of the script:
SCRIPT_PATH="`dirname \"$0\"`"              # relative
SCRIPT_PATH="`( cd \"$SCRIPT_PATH\" && pwd )`"  # absolutized and normalized
if [ -z "$SCRIPT_PATH" ] ; then
  # error; for some reason, the path is not accessible
  # to the script (e.g. permissions re-evaled after suid)
  exit 1  # fail
fi

#------ structure of the path of the examples folder ------

# path of the example directory:
EXAMPLES_DIR=$SCRIPT_PATH/..
# path to the directory containing the parameters for the Fisher matrix calculator:
PARAMETERS_DIR=$EXAMPLES_DIR/parameters
# path to the directory with the analysis parameters:
PARAMETERS_ANALYSIS_DIR=$EXAMPLES_DIR/parameters_analysis
# path to the raw results folder:
RAW_RESULTS_DIR=$EXAMPLES_DIR/raw_results
# path to the results directory:
RESULTS_DIR=$EXAMPLES_DIR/results

#------ where to find the relevant stuff in the CosmicFish library ------

if [ -z "$COSMICFISH_DIR" ]; then
    COSMICFISH_PATH=$SCRIPT_PATH/../..
else
	COSMICFISH_PATH="$COSMICFISH_DIR"
fi  

#------ path to the executable of Fisher matrix:
# 1- CAMB:
FISHER_CAMB_CALCULATOR=$COSMICFISH_PATH/build/fisher_camb/bin/fisher_matrix_calculator.x
# 2- EFTCAMB:
FISHER_EFTCAMB_CALCULATOR=$COSMICFISH_PATH/build/fisher_eftcamb/bin/fisher_matrix_calculator.x
# 3- MGCAMB:
FISHER_MGCAMB_CALCULATOR=$COSMICFISH_PATH/build/fisher_mgcamb/bin/fisher_matrix_calculator.x

# path to the CosmicFish Python library applications:
PYTHON_APPS=$COSMICFISH_PATH/python/apps

#------ path to the CosmicFish bash scripts:
COSMICFISH_SCRIPT_DIR=$COSMICFISH_PATH/script

# import colors:
source $COSMICFISH_SCRIPT_DIR/colors.sh
# import error checking function:
source $COSMICFISH_SCRIPT_DIR/error_check.sh


