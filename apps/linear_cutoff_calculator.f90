!----------------------------------------------------------------------------------------
!
! This file is part of CosmicFish.
!
! Copyright (C) 2015-2017 by the CosmicFish authors
!
! The CosmicFish code is free software;
! You can use it, redistribute it, and/or modify it under the terms
! of the GNU General Public License as published by the Free Software Foundation;
! either version 3 of the License, or (at your option) any later version.
! The full text of the license can be found in the file LICENSE at
! the top level of the CosmicFish distribution.
!
!----------------------------------------------------------------------------------------

!> @file linear_cutoff_calculator.f90
!! This file contains a program that, exploiting the CosmicFish library, computes the
!! maximum multipole to use for linear theory to be valid.

!> @author Marco Raveri

program linear_cutoff_calculator

    use precision
    use ModelParams
    use AMLutils
    use IniFile
    use cosmicfish_utilities
    use cosmicfish_types
    use init_from_file
    use Fisher_calculator_Cls
    use Fisher_calculator_Cls_linearization

    implicit none

    Type(CAMBparams)        :: P
    Type(cosmicfish_params) :: FP
    character(LEN=Ini_max_string_len) paramfile, command_line_k_cutoff

    real(dl) :: k_cutoff
    integer  :: error, cl_dim, AllocateStatus, i
    integer, allocatable, dimension(:) :: l_cut_halofit, l_cut_precision, l_cut_all

    real(dl), parameter :: atol = 0.0_dl
    real(dl), parameter :: rtol = 1.0_dl

    ! feedback:
    call CosmicFish_write_header( name=' CosmicFish linear cutoff calculator. ' )

    ! read the parameter file name:
    paramfile = ''
    if (GetParamCount() /= 0)  paramfile = GetParam(1)
    if (paramfile == '') stop 'No parameter input file'

    ! initialize parameters:
    call init_cosmicfish_from_file( P, FP, paramfile )

    ! read the (optional) k_max cutoff:
    command_line_k_cutoff = ''
    k_cutoff              = 0.1_dl
    if (GetParamCount() > 1)  then
        command_line_k_cutoff = GetParam(2)
        read(command_line_k_cutoff,*) k_cutoff
    end if

    ! check whether the user wants the Cls Fisher:
    if  ( .not. FP%cosmicfish_want_cls ) then
        if ( FP%cosmicfish_feedback > 0 ) then
            write(*,*) 'No window found in the parameter file. Exiting.'
        end if
        stop 0
    end if

    ! get the dimension of the Cls covariance:
    call dimension_cl_covariance( P, FP, cl_dim )

    ! allocate the l cut vector:
    allocate( l_cut_halofit(cl_dim), l_cut_precision(cl_dim), l_cut_all(cl_dim), stat = AllocateStatus)
    if (AllocateStatus /= 0) stop "Allocation failed: l_cut."

    ! call the linearization:
    if ( FP%cosmicfish_feedback >= 1 ) then
        write(*,'(a)')
        write(*,'(a)') 'Halofit linearization:'
        write(*,'(a)')
    end if
    call halofit_linearization( P, FP, rtol, atol, cl_dim, l_cut_halofit, error, outroot=FP%outroot )
    ! print to file the output:
    if ( FP%cosmicfish_feedback >= 1 ) then
        write(*,'(a)')
        call print_linearization( P, FP, l_cut_halofit, outroot=FP%outroot//'halofit_linearization.ini' )
    end if

    ! apply the l cut:
    call apply_l_cut( P, FP, l_cut_halofit )

    ! call the precizion test:
    if ( .False. ) then
        if ( FP%cosmicfish_feedback >= 1 ) then
            write(*,'(a)')
            write(*,'(a)') 'Precision l cut:'
            write(*,'(a)')
        end if
        call precision_settings_l_cut( P, FP, rtol, atol, cl_dim, l_cut_precision, error, outroot=FP%outroot )
        ! print to file the output:
        if ( FP%cosmicfish_feedback >= 1 ) then
            write(*,'(a)')
            call print_linearization( P, FP, l_cut_precision, outroot=FP%outroot//'precision_l_cut.ini' )
        end if
    else
        l_cut_precision = l_cut_halofit
    end if

    ! print feedback:
    if ( FP%cosmicfish_feedback >= 1 ) then
        write(*,'(a)')
        write(*,'(a)') 'Global linearization results:'
        write(*,'(a)')
    end if
    ! get the most conservative:
    do i=1, cl_dim
        l_cut_all(i) = min( l_cut_halofit(i), l_cut_precision(i) )
    end do
    ! print to file and screen the l cut:
    if ( error==0 ) then
        call print_linearization( P, FP, l_cut_all, outroot=FP%outroot//'linearization.ini' )
    end if

end program linear_cutoff_calculator
