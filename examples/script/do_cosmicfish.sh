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
	
	# use camb if none is specified:
	if [ "$filename" == "$b_code" ]; then
		BIN_DIR=$FISHER_CAMB_BIN;
	fi
	# use camb if specified:
	if [ "$b_code" == 'camb' ]; then
		BIN_DIR=$FISHER_CAMB_BIN;
	fi
	# use eftcamb if specified:
	if [ "$b_code" == 'eftcamb' ]; then
		BIN_DIR=$FISHER_EFTCAMB_BIN;
	fi	
	# use mgcamb if specified: 
	if [ "$b_code" == 'mgcamb' ]; then
		BIN_DIR=$FISHER_MGCAMB_BIN;
	fi
		
	echo
	echo -e ${BRed}'Doing: ' $filename 'with' $b_code${Color_Off}
	echo
	
	if [ -z "$APP" ]; then
	    USE_APP=( fisher_matrix_calculator )
	else
		if [ "$APP" == "all" ]; then
			USE_APP=(${APPLICATIONS[*]})
		else
			USE_APP=($APP)
		fi	
	fi

	# cycle over applications:	
	for run_app in "${USE_APP[@]}"
	do
		$BIN_DIR$run_app.x $params;
	done

done;

#------ exit without error ------
exit 0