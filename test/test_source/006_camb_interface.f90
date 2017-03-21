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

program test_camb_interface

    use precision
    use cosmicfish_types
    use init_from_file
    use ModelParams
    use AMLutils
    use IniFile
    use cosmicfish_camb_interface

    implicit none

    Type(CAMBparams)       :: P
    Type(cosmicfish_params) :: FP
    character(LEN=Ini_max_string_len) paramfile, param_name

    integer :: num_param, i
    real(dl), dimension(:), allocatable :: params_array, params_array_2



    ! read the parameter file name:
    paramfile = ''
    if (GetParamCount() /= 0)  paramfile = GetParam(1)
    if (paramfile == '') stop 'No parameter input file'
    ! initialize parameters:
    call init_cosmicfish_from_file( P, FP, paramfile )
    ! get the number of parameters:
    call num_param_Fisher(P, FP, num_param)
    ! allocate the parameter vector:
    allocate( params_array(num_param), params_array_2(num_param) )
    ! write the P and FP values of parameters to the parameter vector:
    call CAMB_params_to_params_array( P, FP, num_param, params_array )

    ! change the params_array:
    params_array = 1.1_dl*params_array
    ! call the inverse function a couple of times to make sure that everything
    ! works as expected:
    call params_array_to_CAMB_param( P, FP, num_param, params_array )

    call CAMB_params_to_params_array( P, FP, num_param, params_array_2 )

    do i=1, num_param
        if ( params_array(i) /= params_array_2(i) ) then
            print*, 'Params array to CAMB and its inverse not working for param: ', i
        end if
    end do

    open (unit = 100, file = "output")

    write(100,*) 'Number of parameters: ', num_param
    write(100,*)
    write(100,*) '1.1* fiducial parameters: '
    write(100,*)

    do i=1, num_param

        call Fisher_param_names(P, FP, i, param_name)

        write(100,*) i, trim(param_name), params_array(i), params_array_2(i)

    end do

    write(100,*)
    write(100,*) 'fiducial parameters: '
    write(100,*)

    do i=1, num_param

        call Fisher_param_names(P, FP, i, param_name)

        write(100,*) i, trim(param_name), params_array(i)/(1.1_dl)

    end do

    close(100)

    ! all the successful test should return the following.
    stop 100

end program test_camb_interface
