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

program test_fisher_SN_read

    use precision
    use cosmicfish_types
    use init_from_file
    use ModelParams
    use AMLutils
    use IniFile
    use cosmicfish_camb_interface
    use Fisher_calculator_Cls
    use Fisher_calculator_SN
    use Fisher_manipulation

    implicit none

    Type(CAMBparams)       :: P
    Type(cosmicfish_params) :: FP
    character(LEN=Ini_max_string_len) paramfile

    integer  :: num_param
    real(dl), allocatable :: Fisher_matrix(:,:)

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

    if ( FP%cosmicfish_want_SN  ) then
        call Fisher_SN( P, FP, num_param, Fisher_Matrix )
        call Fisher_protection_unconstrained( Fisher_Matrix )
        call Fisher_protection_degenerate( Fisher_Matrix )
    end if

    ! have to set this to zero because the values are never the same...
    Fisher_Matrix = 0._dl
    call save_Fisher_to_file( P, FP, num_param, Fisher_Matrix, filename='output' )

    ! all the successful test should return the following.
    stop 100

end program test_fisher_SN_read
