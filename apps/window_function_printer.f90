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
    use Fisher_calculator_Cls_windows

    implicit none

    Type(CAMBparams)        :: P
    Type(cosmicfish_params) :: FP
    character(LEN=Ini_max_string_len) paramfile

    integer  :: num_points = 1000
    real(dl) :: z_min      = 0._dl
    real(dl) :: z_max      = 10._dl

    integer  :: error, i, j
    real(dl) :: redshift, winamp
    Type(TRedWin) , pointer :: RedWin

    ! feedback:
    call CosmicFish_write_header( name=' CosmicFish survey window function printer. ' )

    ! read the parameter file name:
    paramfile = ''
    if (GetParamCount() /= 0)  paramfile = GetParam(1)
    if (paramfile == '') stop 'No parameter input file'

    ! initialize parameters:
    call init_cosmicfish_from_file( P, FP, paramfile )

    ! check whether the user wants the window functions:
    if  ( .not. FP%cosmicfish_want_cls .or. &
        & .not. ( FP%fisher_cls%Fisher_want_LSS_lensing .or. FP%fisher_cls%Fisher_want_LSS_counts ).or. &
        & FP%fisher_cls%LSS_number_windows == 0 ) then

        if ( FP%cosmicfish_feedback > 0 ) then
            write(*,*) 'No window found in the parameter file. Exiting.'
        end if

        stop

    end if

    ! initialize the window functions:
    call init_camb_sources_windows( FP )

    ! set camb parameters:
    call CAMBParams_Set(P,error)
    ! check the parameter set:
    if ( error /= 0 ) then
        write(*,*) 'ERROR: cannot set parameters. Calculation cannot proceed.'
        stop 1
    end if

    ! open the output file:
    open(unit=666, FILE=FP%outroot//'windows.dat', ACTION="write", STATUS="replace")

    ! write the header:
    write(666,'(a)') '#'
    write(666,'(a)') '# This file contains the window function of a CosmicFish survey.'
    write(666,'(a,F6.1,a,F6.1,a,I6,a)') '# Windows printed from redshift:', z_min, ' to redshift', z_max, ' on ', num_points, ' points.'
    write(666,'(a)') '#'

    ! write the number of windows:
    write(666,'(a)', advance='no') '# z'//'    '
    do i=1, FP%fisher_cls%LSS_number_windows
        write(666,'(a)', advance='no') 'W_'//TRIM(integer_to_string(i))//'    '
    end do
    write(666,'(a)')


    ! do the actual window function sampling:
    do i=1, num_points

        redshift = z_min +REAL( i -1 )/REAL( num_points -1 )*( z_max -z_min )

        write(666, '(E15.5)', advance='no') redshift
        do j=1, FP%fisher_cls%LSS_number_windows
            RedWin => Redshift_W(j)
            write(666,'(E15.5E3)', advance='no') counts_background_z( Win=RedWin, z=redshift )
        end do
        write(666,'(a)')

    end do

    ! close the file:
    close(666)

end program window_function_printer
