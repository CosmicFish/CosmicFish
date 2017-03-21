====================
CosmicFish: examples
====================

This folder contains the relevant code to produce a set of examples with the
CosmicFish code.

The idea is to use the experiments used to validate the code to show the user
an example of a pipeline completely exploiting the CosmicFish code.

To produce all examples just issue::

	make all

in the examples folder. This is equivalent to::

  make examples

in the main CosmicFish folder.

1. Creating examples:
=====================

The script in the examples folder have a structure that is flexible enough so
that the user can use it to create his own pipeline.

Just copy the examples folder in the main CosmicFish directory, replace the
parameters and the analysis parameters and run with ``make all``.

You can also copy the examples folder somewhere else on your system. In this
case minimum modifications to the scripts are necessary to instruct the code
to find the CosmicFish library.

If in the system an environment variable called ``$COSMICFISH_DIR`` is found
then no modification is required.

Otherwise just go to ``examples/script/common.sh`` and point ``COSMICFISH_PATH``
to the right direction.

Notice that the example scripts will select the Boltzmann code based on the extension
of the parametere file as follows::

   parameters_name.camb.ini    : will use CAMB
   parameters_name.eftcamb.ini : will use EFTCAMB
   parameters_name.mgcamb.ini  : will use MGCAMB

On the other hand, for maximum flexibility, the analysis parameters do not need the
Boltzmann code extension.

2. Makefile targets:
====================

The Makefile has several targets. Here's a breakdown of the main ones:

**Main**:

* ``all``: runs all examples. Creates the Fisher matrices and then analyses them;
* ``fisher_matrices``: just creates all the Fisher matrices;
* ``analysis``: analyses all the Fisher matrices;
* ``clean``: removes just raw results but leaves analysis results there;
* ``deep_clean``: removes all results;

**Secondary**:

* ``targets``: creates Makefile targets with the names of the parameters. After issuing this one can run only one choice of parameters;
* ``additional_script``: this runs all the python script in the script folder. By default does nothing but the user might use it to include his own python script in his pipeline;
