!----------------------------------------------------------------------------------------
!
! This file is part of CosmicFish.
!
! Copyright (C) 2015-2016 by the CosmicFish authors
!
! The CosmicFish code is free software;
! You can use it, redistribute it, and/or modify it under the terms
! of the GNU General Public License as published by the Free Software Foundation;
! either version 3 of the License, or (at your option) any later version.
! The full text of the license can be found in the file LICENSE at
! the top level of the CosmicFish distribution.
!
!----------------------------------------------------------------------------------------

!> @file fisher_matrix_calculator.f90
!! This file contains a program that, exploiting the CosmicFish library, creates Fisher
!! matrices for several cosmological experiments.

!> @author Marco Raveri and Matteo Martinelli

program fisher_matrix_calculator

    use precision
    use cosmicfish_types
    use init_from_file
    use ModelParams
    use AMLutils
    use IniFile
    use cosmicfish_camb_interface
    use Fisher_calculator_Cls
    use Fisher_calculator_SN
    use Fisher_calculator_RD
    use Fisher_calculator_Derived
    use Fisher_manipulation
    use cosmicfish_utilities

    implicit none

    Type(CAMBparams)        :: P
    Type(cosmicfish_params) :: FP
    character(LEN=Ini_max_string_len) paramfile

    integer  :: num_param, num_derived
    real(dl), allocatable, dimension(:,:) :: Fisher_matrix, Derived_Matrix
    real(dl), allocatable, dimension(:)   :: fiducial_params

    ! feedback:
    call CosmicFish_write_header( name=' CosmicFish Fisher matrix calculator. ' )

    ! read the parameter file name:
    paramfile = ''
    if (GetParamCount() /= 0)  paramfile = GetParam(1)
    if (paramfile == '') stop 'No parameter input file'
    ! initialize parameters:
    call init_cosmicfish_from_file( P, FP, paramfile )
    ! get the number of parameters:
    call num_param_Fisher(P, FP, num_param)
    ! allocation:
    allocate( Fisher_matrix(num_param,num_param) )

    Fisher_matrix = 0._dl

    ! Cls Fisher matrix:
    if ( FP%cosmicfish_want_cls ) then
        ! compute the Fisher matrix:
        call Fisher_Cls( P, FP, num_param, Fisher_Matrix, outroot=FP%outroot )
        ! save it to file:
        call save_Fisher_to_file( P, FP, num_param, Fisher_Matrix, filename=FP%outroot//'fisher_matrix_cls.dat' )
        ! save the parameter names to file:
        call save_paramnames_to_file( P, FP, num_param, filename=FP%outroot//'fisher_matrix_cls.paramnames' )
    end if

    if ( FP%cosmicfish_want_SN  ) then
        ! compute the Fisher matrix:
        call Fisher_SN( P, FP, num_param, Fisher_Matrix, outroot=FP%outroot )
        ! save it to file:
        call save_Fisher_to_file( P, FP, num_param, Fisher_Matrix, filename=FP%outroot//'fisher_matrix_SN.dat' )
        ! save the parameter names to file:
        call save_paramnames_to_file( P, FP, num_param, filename=FP%outroot//'fisher_matrix_SN.paramnames' )
    end if

    if ( FP%cosmicfish_want_RD  ) then
        ! compute the Fisher matrix:
        call Fisher_RD( P, FP, num_param, Fisher_Matrix, outroot=FP%outroot)
        ! save it to file:
        call save_Fisher_to_file( P, FP, num_param, Fisher_Matrix, filename=FP%outroot//'fisher_matrix_RD.dat' )
        ! save the parameter names to file:
        call save_paramnames_to_file( P, FP, num_param, filename=FP%outroot//'fisher_matrix_RD.paramnames' )
    end if

    if ( FP%cosmicfish_want_Mpk ) write(*,*) 'Not yet implemented'

    if (FP%cosmicfish_want_derived ) then
        ! get number of derived parameters:
        call dimension_derived_parameters( P, FP, num_derived )
        ! allocate the derived parameters matrix:
        allocate( Derived_Matrix(num_param, num_derived) )
        ! compute the Fisher matrix:
        call Fisher_derived( P, FP, num_param, num_derived, Derived_Matrix, outroot=FP%outroot )
        ! save it to file:
        call save_derived_Fisher_to_file( P, FP, num_param, num_derived, Derived_Matrix, filename=FP%outroot//'fisher_matrix_derived.dat' )
        ! save the parameter names to file:
        call save_derived_paramnames_to_file( P, FP, num_param, num_derived, filename=FP%outroot//'fisher_matrix_derived.paramnames' )
    end if

end program fisher_matrix_calculator
