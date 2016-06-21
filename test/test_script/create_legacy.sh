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

# this script takes care of running the test and comparing the output to the
# output that is already provided.

# decide script options:
NUMDIFF_OPTIONS=-q

# parse command line options:
while [[ $# > 0 ]]
do
key="$1"
case $key in
    -v|--verbose)
    NUMDIFF_OPTIONS=-V
    shift
    ;;
    *)
            # unknown option
    ;;
esac
shift
done

# colors:
RED='\033[0;31m'   # Red
GREEN='\033[0;32m' # Green
BLUE='\033[0;34m'  # Blue
NC='\033[0m'       # No Color

# get the path of the script:
SCRIPT_PATH="`dirname \"$0\"`"              # relative
SCRIPT_PATH="`( cd \"$SCRIPT_PATH\" && pwd )`"  # absolutized and normalized
if [ -z "$SCRIPT_PATH" ] ; then
  # error; for some reason, the path is not accessible
  # to the script (e.g. permissions re-evaled after suid)
  exit 1  # fail
fi

# get the path of the test suite:
TEST_DIR=$SCRIPT_PATH/../test_build
TEST_INPUT_DIR=$SCRIPT_PATH/../test_input
TEST_OUTPUT_DIR=$SCRIPT_PATH/../test_output

# some feedback:
printf "${GREEN}********************************************\n"
printf "Creating legacy results for CosmicFish tests:\n"
printf "********************************************\n${NC}"

# go to the test folder:
cd $TEST_DIR

# cycle over the test:
for file in *.x
do

  # get the name of the test:
  TEST_NAME=$(basename "$file")
  TEST_NAME="${TEST_NAME%.*}"

  # feedback:
  echo 'Doing test '$TEST_NAME

  # run the test. If an input file is provided then feed it to the test otherwise not.
  if [ -e "$TEST_INPUT_DIR/$TEST_NAME.input" ]
  then
    ./$file $TEST_INPUT_DIR/$TEST_NAME.input > /dev/null 2>&1
  else
    ./$file > /dev/null 2>&1
  fi

  mv -f output $TEST_OUTPUT_DIR/$TEST_NAME.output

  echo '********************************************'

done

# Finishing:
rm -rf output


#
