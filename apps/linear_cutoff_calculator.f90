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

    implicit none

    Type(CAMBparams)        :: P
    Type(cosmicfish_params) :: FP
    character(LEN=Ini_max_string_len) paramfile, command_line_k_cutoff

    real(dl) :: k_cutoff
    integer  :: error

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

    ! set camb parameters:
    call CAMBParams_Set(P,error)
    ! check the parameter set:
    if ( error /= 0 ) then
        write(*,*) 'ERROR: cannot set parameters. Calculation cannot proceed.'
        stop 1
    end if


end program linear_cutoff_calculator
