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
# This file contains part of the CosmicFish build system: compiler and compiler options
#


# detect the fortran compiler: if ifort is present go for it, otherwise use gfortran.
HAS_IFORT = $(shell which ifort >/dev/null; echo $$?)

ifeq "$(HAS_IFORT)" "0" 
	# use ifort
    export F90C := ifort
    export F90CRLINK := -cxxlib
    # release options:
	FFLAGS_RELEASE   := -qopenmp -mkl=parallel -O3 -W0 -WB -fpp -qopt-report=0 
	# debug options:
	FFLAGS_DEBUG     := -qopenmp -mkl=parallel -fpp -g -qopt-report=0 -fp-stack-check -O0 -traceback -check all -check bounds -check uninit -check noarg_temp_created -DDEBUG
	# module flag:
	MODULE_FLAG      := -module 
	# lapack linking
	LAPACK_LINK      := 
	LAPACK_INCLUDE   :=
else 
	# use gfortran
    export F90C := gfortran
    export ifortErr := 1
	# release options:
	FFLAGS_RELEASE   := -O3 -fopenmp -ffast-math -fmax-errors=4 -ffree-line-length-none -funroll-loops -cpp -ffpe-summary=none
	# debug options:
	FFLAGS_DEBUG     := -O0 -g -fopenmp -fmax-errors=0 -ffree-line-length-none -cpp -fbounds-check -fbacktrace -ffpe-trap=invalid,overflow,zero -DDEBUG
	# module flag:
	MODULE_FLAG := -J
	# lapack linking
	LAPACK_LINK      := -llapack -lblas
	LAPACK_INCLUDE   :=
endif

# wether to compile with debug options or not:
DEBUG := 0
ifeq ($(DEBUG), 1)
    export FFLAGS = $(FFLAGS_DEBUG)
else
    export FFLAGS = $(FFLAGS_RELEASE)
endif
