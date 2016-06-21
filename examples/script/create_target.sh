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
# This file contains the bash script that creates a makefile containing the targets for specific parameters
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

#------ create the Makefile with the expanded targets ------

# delete preexisting file:
rm -rf $SCRIPT_PATH/makefile_target.make
# create the new file:
touch  $SCRIPT_PATH/makefile_target.make

full_targets=()

for params in $PARAMETERS_DIR/*.ini;

do

filename=$(basename "$params")
extension="${filename##*.}"
filename="${filename%.*}"

full_targets+=($filename)

echo $filename':' >>$SCRIPT_PATH/makefile_target.make
echo '	-@bash $(SCRIPT_DIR)/run_params.sh $(PARAMETER_DIR)'/$filename.$extension '$(PARAMETER_ANALYSIS_DIR)'/$filename.$extension >> $SCRIPT_PATH/makefile_target.make

done;

# do the global target (usefull to launch in parallel)
printf 'all_targets: ' >>$SCRIPT_PATH/makefile_target.make

for ind in `seq 1 ${#full_targets[@]}`
do
printf ${full_targets[$ind-1]}' ' >>$SCRIPT_PATH/makefile_target.make
done;

#------ exit without error ------
exit 0
