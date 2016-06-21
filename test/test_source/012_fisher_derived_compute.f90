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

program test_fisher_derived_compute

    use precision
    use cosmicfish_types
    use init_from_file
    use ModelParams
    use AMLutils
    use IniFile
    use cosmicfish_camb_interface
    use Fisher_calculator_Derived

    implicit none

    Type(CAMBparams)       :: P
    Type(cosmicfish_params) :: FP
    character(LEN=Ini_max_string_len) paramfile

    integer  :: num_param, num_derived, i, err
    character(LEN=500) :: param_name
    real(dl), allocatable, dimension(:) :: derived

    ! read the parameter file name:
    paramfile = ''
    if (GetParamCount() /= 0)  paramfile = GetParam(1)
    if (paramfile == '') stop 'No parameter input file'
    ! initialize parameters:
    call init_cosmicfish_from_file( P, FP, paramfile )
    ! get the number of parameters:
    call num_param_Fisher(P, FP, num_param)
    ! get the number of derived parameters:
    call dimension_derived_parameters( P, FP, num_derived )

    allocate( derived(num_derived) )

    open (unit = 100, file = "output")

    call Compute_Derived_Parameters( P, FP, derived, num_derived, err )

    if ( err/=0 ) then
        write(100,*) 'Error with the computation of derived parameters'
        stop
    end if

    do i=1, num_derived
        call Fisher_derived_param_names(P, FP, i, param_name)
        write(100,*) i, trim(param_name), ' = ', derived(i)
    end do


    close(100)

    ! all the successful test should return the following.
    stop 100

end program test_fisher_derived_compute
