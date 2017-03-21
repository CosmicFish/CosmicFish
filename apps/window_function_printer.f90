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

!> @file window_function_printer.f90
!! This file contains a program that, exploiting the CosmicFish library, prints the
!! survey dN/dz to file.

!> @author Marco Raveri

program window_function_printer

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
    character(LEN=Ini_max_string_len) paramfile

    integer  :: num_points = 1000
    real(dl) :: z_min      = 0._dl
    real(dl) :: z_max      = 100._dl

    integer :: error, i

    ! feedback:
    call CosmicFish_write_header( name=' CosmicFish survey window function printer. ' )

    ! read the parameter file name:
    paramfile = ''
    if (GetParamCount() /= 0)  paramfile = GetParam(1)
    if (paramfile == '') stop 'No parameter input file'

    ! initialize parameters:
    call init_cosmicfish_from_file( P, FP, paramfile )

    ! set camb parameters:
    call CAMBParams_Set(P,error)
    ! check the parameter set:
    if ( error /= 0 ) then
        write(*,*) 'ERROR: cannot set parameters. Calculation cannot proceed.'
        stop 1
    end if

    ! compute the window functions on a z grid:
    do i = 1, num_points

    end do

    ! print the file header:

    ! open txt file
    ! print the header with non-advancing

    ! print the data:

    ! close the file:

    ! exit with error: the application is yet not fully developed
    write(*,*) 'ERROR: the application is yet not fully developed.'

end program window_function_printer
