#!/bin/bash

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
# This file contains the bash script with the sequence of operations necessary to produce
# the validation results for Planck Pre Launch forecasts
#

#------ common initialization ------
# get the path of the script:
SCRIPT_PATH="`dirname \"$0\"`"              # relative
SCRIPT_PATH="`( cd \"$SCRIPT_PATH\" && pwd )`"  # absolutized and normalized
if [ -z "$SCRIPT_PATH" ] ; then
  # error; for some reason, the path is not accessible
  # to the script (e.g. permissions re-evaled after suid)
  exit 1  # fail
fi

source $SCRIPT_PATH/common.sh

#------ run the analysis using the python application: full_analysis.py ------

for params in $PARAMETERS_ANALYSIS_DIR/*.ini;

do python $PYTHON_APPS/full_analysis.py $params ;

done;

#------ exit without error ------
exit 0
