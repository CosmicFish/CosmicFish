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
# This file contains the bash script that creates Fisher matrices and runs full analysis on the results
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

#------ digest the input ------

# get the file name for the Fisher parameters:
filename=$(basename "$1")
extension="${filename##*.}"
filename="${filename%.*}"
# get the Boltzmann solver:
b_code="${filename##*.}"

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

# get the file name for the analysis parameters:
dirname_analysis=$(dirname "$2")
filename_analysis=$(basename "$2")
extension_analysis="${filename_analysis##*.}"
filename_analysis="${filename_analysis%.*}"
b_code_temp="${filename_analysis##*.}"
filename_analysis="${filename_analysis%.*}"

#------ create the Fisher matrix ------

$SOLVER  $1 ;

#------ run analysis ------

python $PYTHON_APPS/full_analysis.py $dirname_analysis/$filename_analysis.$extension_analysis ;

#------ exit without error ------
exit 0