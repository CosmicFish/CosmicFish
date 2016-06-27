CosmicFish: a cosmology forecasting tool
========================================

[![Build Status](https://travis-ci.org/CosmicFish/CosmicFish.svg?branch=master)](https://travis-ci.org/CosmicFish/CosmicFish)

This library contains a set of tools and applications to
perform forecast of the accuracy of future cosmological
experiments.

If you find this collection useful, feel free to download, use it and suggest
pull requests!

The official distribution is available for download at:

http://cosmicfish.github.io

while the official developers version is available on GitHub, and you can clone the
repository using:

	git clone https://github.com/CosmicFish/CosmicFish

### 1. CosmicFish Requirements:

Currently, the distributed source code requires the intel fortran compiler and some other stuff.
Here's a list of dependencies for a full build:

	ifort or gfortran (5.3.0)
	MKL or BLAS/LAPACK
	numdiff
	doxygen

If you use ifort the code will automatically try to link mkl while if you use gfortran you need BLAS and LAPACK libraries installed. You can do so by ``sudo apt-get install liblapack-dev``. If both ifort and gfortran are installed the code will use ifort. Compilation options can be found in ``make_macros/compilers_options.make``.

Notice that only the 5.3.0 version of gfortran is supported as previous versions missed some OOP features. You should (hopefully) be able to get this version of the code by:

	sudo add-apt-repository ppa:ubuntu-toolchain-r/test
	sudo apt-get update
	sudo apt-get install gcc-5 g++-5 gfortran-5

The python library then requires several other libraries to work. We warmly suggest to install a
bundled package like (https://www.continuum.io/downloads). Just in case here's a list of all python
dependencies:

	python 2.7
	python library matplotlib
	python library os
	python library scipy
	python library numpy
	python library math
	python library itertools
	python library copy

	nosetests for python testing
	with coverage (pip install coverage should work)
	sphinx (http://www.sphinx-doc.org) for python documentation

Some of the requirements of the library are included and will be compiled along with the library if not present on the system.
These can de found in the folder:

	cosmicfish/bundled

### 2. Installation procedure:

After downloading the CosmicFish library it can be compiled by running:

	make cosmicfish

and tested using:

	make run_test

If you want to build the CosmicFish library in parallel, using multiple processors, just add ``-jN`` to the make commands. If you use ``-j`` the Makefile will use all available processors.

### 3. Examples and packages

After the library is built it can be directly used to produce some examples with:

	make examples

This will produce a set of example results in the ``examples`` folder. These are meant to outline a standard workflow with the CosmicFish library.

The structure of the example folder can be used to produce CosmicFish packages that naturally exploit the set up workflow.
To do so just copy (and rename) the example folder somewhere on your system then tell the system where the CosmicFish library is built with:

	export COSMICFISH_DIR=path_to_cosmicfish_build

Then replace the parameters and the analysis parameters and run everything with ``make all`` in your directory.

For further details please check the examples directory readme.

### 4. Extensive documentation:

The mathematical and physical details of what is implemented in the CosmicFish library are described in the implementation notes:

    arXiv:1606.06268 CosmicFish Implementation Notes V1.0

A browsable documentation of the fortran library is available at:

http://cosmicfish.github.io/documentation/CosmicFish/index.html

In addition the documentation of the Python library to manipulate Fisher matrices is available at:

http://cosmicfish.github.io/documentation/CosmicFishPyLib/index.html

### 5. Citing this work

If you use this library, please consider citing the following works:

	arXiv:1606.06273 Information Gain in Cosmology: From the Discovery of Expansion to Future Surveys
	arXiv:1606.06268 CosmicFish Implementation Notes V1.0

If the results were obtained with the CosmicFish code built with CAMB, EFTCAMB or MGCAMB in addition to the previous references the standard references for these codes should be cited as well.

This is the usual, fair way of giving credit to contributors to a scientific result. In addition, it helps us justify our effort in developing the CosmicFish code as an academic undertaking.

### 6. Licence Information

See the LICENSE file in this directory.

### 7. Build system target

The Makefile has several targets. Here's a breakdown of the main ones:

**Main**:

* ``all``: builds all the library, runs all the tests and builds the documentation;
* ``cosmicfish``: builds all the library and the apps, builds the test without running them and builds the documentation;
* ``default``: the default target is cosmicfish;
* ``clean``: removes all built material leaving the structure intact;
* ``distclean``: removes everything;

**Secondary**:

* ``apps``: compiles the CosmicFish apps;
* ``bsolvers``: compiles all the Boltzmann solvers;
* ``bundled``: compiles all bundled software;
* ``directories``: creates the directory structure of the CosmicFish library;
* ``documentation``: builds both the doxygen and the python documentation;
* ``test``: builds the CosmicFish tests;
* ``run_test``: runs both the CosmicFish tests and Python tests;

**Some options**:

Compile with debug options:

	make target DEBUG=1

Change the build directory:

	make target BUILD_DIR=wherever_you_want_to_build
