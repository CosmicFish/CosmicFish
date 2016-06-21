===============
EFTCAMB sources
===============

Thanks for downloading EFTCAMB sources! This is EFTCAMB_Oct15.

In order to build EFTCAMB, please follow these 6 steps:

1. Download CAMB source (Sources_Aug15) from: https://github.com/cmbant/CAMB/tree/CAMB_sources

2. Add these new files into your CAMB-sources folder:

   EFT_def.90
   EFT_designer.f90
   EFT_functions.f90
   EFT_Horndeski.f90
   EFT_main.f90
   equations_EFT.f90

3. Replace the default files with the following files:

   cmbmain.f90
   inidriver.F90
   lensing.f90
   Makefile
   Makefile_main
   modules.f90
   params.ini
   reionization.f90
   subroutines.f90

4. configure the Makefile according to your Fortran compiler.

5. $make

6. $./camb params.ini

--
EFTCAMB team
Apr/2016

====
CAMB
====
:CAMB:  Code for Anisotropies in the Microwave Background, Fortran 95 code
:Homepage: http://camb.info/

Description and installation
============================

For full details see the `ReadMe <http://camb.info/readme.html>`_.

Branches
========

The master branch contains latest changes to the main release version.

CAMB_sources is the updated public `CAMB Sources <http://camb.info/sources/>`_ code.

The vehre_formatted branch is a developement version, which integrates CAMB and CAMB sources, and uses Fortran 2008 (and hence requires gfortran 6).
This branch also has an integrated test suite, which runs on Travis automatically for new commits.
