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

!> @file fiducial_calculator.f90
!! This file contains a program that, exploiting the CosmicFish library, prints to file
!! the fiducial model for the desired experimental configuration.

!> @author Marco Raveri

program fiducial_calculator

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

    ! feedback:
    call CosmicFish_write_header( name=' CosmicFish fiducial model printer. ' )

    ! read the parameter file name:
    paramfile = ''
    if (GetParamCount() /= 0)  paramfile = GetParam(1)
    if (paramfile == '') stop 'No parameter input file'
    ! initialize parameters:
    call init_cosmicfish_from_file( P, FP, paramfile )

    ! Cls Fisher matrix:
    if ( FP%cosmicfish_want_cls ) then
        write(*,*) 'WARNING: Cls fiducial printing not yet implemented.'
    end if

    if ( FP%cosmicfish_want_SN  ) then
        write(*,*) 'WARNING: SN fiducial printing not yet implemented.'
    end if

    if ( FP%cosmicfish_want_RD  ) then
        write(*,*) 'WARNING: RD fiducial printing not yet implemented.'
    end if

    if ( FP%cosmicfish_want_Mpk ) write(*,*) 'Not yet implemented'

    if (FP%cosmicfish_want_derived ) then
        write(*,*) 'WARNING: Derived Parameters fiducial printing not yet implemented.'
    end if

end program fiducial_calculator
