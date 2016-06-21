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
NUMDIFF_OPTIONS="-V --absolute-tolerance=0.001 --relative-tolerance=0.001"

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

# test if numdiff works:
if ! type numdiff > /dev/null; then
  NUMDIFF_APP=$SCRIPT_PATH/../../bundled/numdiff/numdiff
else
  NUMDIFF_APP=numdiff
fi

# some feedback:
printf "${GREEN}*********************************************************\n"
printf "Running CosmicFish test suite:\n"
printf "*********************************************************\n${NC}"

# go to the test folder:
cd $TEST_DIR

FAILED_TEST=()
FAILURE_REASON=()
TOTAL_TEST=0

# cycle over the test:
for file in *.x
do

  echo '*********************************************************'

  # get the name of the test:
  TEST_NAME=$(basename "$file")
  TEST_NAME="${TEST_NAME%.*}"

  # feedback:
  echo 'Doing test '$TEST_NAME':'
  TOTAL_TEST=$((TOTAL_TEST+1))

  # remove previous test products:
  rm -rf $TEST_NAME.fail
  rm -rf output
  rm -rf test.log

  # Run the test. If an input file is provided then feed it to the test otherwise not.
  # The script tests to see wether the program runned correctly.

  if [ -e "$TEST_INPUT_DIR/$TEST_NAME.input" ]
  then
    # run the test:
    res1=$(date +%s.%N)
    ./$file $TEST_INPUT_DIR/$TEST_NAME.input > test.log 2>&1
    # test the exit status of the program:
    if [ $? -ne 100 ]; then
        printf "${RED}TEST FAILED:${NC} run failure.\n"
        FAILED_TEST+=($TEST_NAME)
        FAILURE_REASON+=('run failed')
        mv test.log $TEST_NAME.fail
        continue
    fi;
    res2=$(date +%s.%N)
    test_time=$(python -c "print(round(${res2} - ${res1},3))")
  else
    # run the test:
    res1=$(date +%s.%N)
    ./$file > test.log 2>&1
    # test the exit status of the program:
    if [ $? -ne 100 ]; then
        printf "${RED}TEST FAILED:${NC} run failure.\n"
        FAILED_TEST+=($TEST_NAME)
        FAILURE_REASON+=('run failed')
        mv test.log $TEST_NAME.fail
        continue
    fi;
    res2=$(date +%s.%N)
    test_time=$(python -c "print(round(${res2} - ${res1},3))")
  fi;

  # test to see wether there is an output file to compare with and do the comparison:
  if [ -e "$TEST_OUTPUT_DIR/$TEST_NAME.output" ]
    then

    if ( $NUMDIFF_APP $TEST_OUTPUT_DIR/$TEST_NAME.output output $NUMDIFF_OPTIONS >> test.log 2>&1 ) ; then
        printf "${GREEN}OK:${NC} test passed succesfully in ${test_time} (s).\n"
    else
        printf "${RED}TEST FAILED:${NC} diff failure.\n"
        FAILED_TEST+=($TEST_NAME)
        FAILURE_REASON+=('diff failed')
        echo >> test.log ; echo >> test.log ;
        echo 'output file:' >> test.log ;
        echo >> test.log ; echo >> test.log ;
        cat output >> test.log ;
        mv test.log $TEST_NAME.fail
        continue
    fi
  else
    printf "${RED}TEST FAILED:${NC} output file not found.\n"
    FAILED_TEST+=($TEST_NAME)
    FAILURE_REASON+=('output not found')
    echo >> test.log ; echo >> test.log ;
    echo 'Output file not found. Test output:' >> test.log ;
    echo >> test.log ; echo >> test.log ;
    cat output >> test.log ;
    mv test.log $TEST_NAME.fail
    continue
  fi

done

printf "${GREEN}*********************************************************\n"
printf "CosmicFish test report:\n"
printf "*********************************************************\n${NC}"

printf "\n Passed test: "$(( $TOTAL_TEST-${#FAILED_TEST[@]} ))" / "$TOTAL_TEST"\n"

if [ "${#FAILED_TEST[@]}" -gt "0" ]; then
  printf "\n Failed test:\n\n"
  for ind in `seq 1 ${#FAILED_TEST[@]}`
  do
  echo " " ${FAILED_TEST[$ind-1]} ": " ${FAILURE_REASON[$ind-1]}
  done;
  printf "\n"
else
  printf " All test successfully passed.\n\n"
fi

# finishing up:
rm -rf output
rm -rf test.log

# return the test result:
if [ "${#FAILED_TEST[@]}" -gt "0" ]; then
	exit 1
else
	exit 0
fi

