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
# This file contains the bash script computing Fisher matrices for all the parameters in
# the parameter directory
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

#------ create the Fisher matrices ------

echo
echo -e ${BRed} 'Producing all examples. This will definitely take a while...' ${Color_Off}
echo

for params in $PARAMETERS_DIR/*.ini;

do

# get the file name for the Fisher parameters:
filename=$(basename "$params")
extension="${filename##*.}"
filename="${filename%.*}"
# get the Boltzmann solver:
b_code="${filename##*.}"
filename="${filename%.*}"

if [ "$filename" == "$b_code" ]; then
	SOLVER=$FISHER_CAMB_CALCULATOR;
fi

if [ "$b_code" == 'camb' ]; then
	SOLVER=$FISHER_CAMB_CALCULATOR;
fi

if [ "$b_code" == 'eftcamb' ]; then
	SOLVER=$FISHER_EFTCAMB_CALCULATOR;
fi	

if [ "$b_code" == 'mgcamb' ]; then
	SOLVER=$FISHER_MGCAMB_CALCULATOR;
fi	

echo
echo -e ${BRed}'Doing: ' $filename 'with' $b_code${Color_Off}
echo

$SOLVER $params ;

done;

#------ exit without error ------
exit 0