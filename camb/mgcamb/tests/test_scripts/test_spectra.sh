#
# Script that runs camb with all the parameters in the tests/parameters/ folder
#

#!/bin/bash

# go to the camb root folder:
cd ..

# compile eftcamb:
make clean
make camb

# create all the spectra by cycling over all the parameter files:
for file in tests/parameters/*
	do

	./camb $file | tee -a tests/results/Spectra_results/spectra.log
	echo ; echo ; echo ;

done

# clean up:
make clean
rm camb
