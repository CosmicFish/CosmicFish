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

program test_fisher_compute_mag

    use precision
    use cosmicfish_types
    use init_from_file
    use ModelParams
    use AMLutils
    use IniFile
    use cosmicfish_camb_interface
    use Fisher_calculator_Cls
    use Fisher_calculator_SN

    implicit none

    Type(CAMBparams)       :: P
    Type(cosmicfish_params) :: FP
    character(LEN=Ini_max_string_len) paramfile

    integer :: num_param, number, AllocateStatus, err
    real(dl), dimension(:), allocatable :: param_array
    real(dl) :: initial_step
    real(dl), allocatable, dimension(:,:) :: sn_mock
    real(dl), allocatable, dimension(:)   :: sn_mock_covariance
    real(dl), allocatable, dimension(:)   :: sn_fiducial, sn_derivative, sn_temp
    integer  :: i,j,k

    integer, parameter :: param_test = 1

    ! read the parameter file name:
    paramfile = ''
    if (GetParamCount() /= 0)  paramfile = GetParam(1)
    if (paramfile == '') stop 'No parameter input file'
    ! initialize parameters:
    call init_cosmicfish_from_file( P, FP, paramfile )
    ! get the number of parameters:
    call num_param_Fisher(P, FP, num_param)
    ! allocate a parameter array:
    allocate( param_array(num_param) )
    ! write the parameters in the parameter array:
    call CAMB_params_to_params_array( P, FP, num_param, param_array )

    number = FP%fisher_SN%total_SN_number

    ! generate the mock data:
    allocate( sn_mock(3,number), stat = AllocateStatus )
    if (AllocateStatus /= 0) stop "Allocation failed. Not enough memory to allocate SN mock."

    sn_mock = 0._dl
    ! we need to generate redshifts:
    i = 1
    do j=1, FP%fisher_SN%number_SN_windows
        do k=1, FP%fisher_SN%SN_number(j)
            sn_mock(1,i) = FP%fisher_SN%SN_redshift_start(j) &
                & +REAL(k-1)*(FP%fisher_SN%SN_redshift_end(j)-FP%fisher_SN%SN_redshift_start(j))/REAL(FP%fisher_SN%SN_number(j)-1)
            i = i + 1
        end do
    end do

    ! generate the mock covariance:
    allocate( sn_mock_covariance(number), stat = AllocateStatus )
    if (AllocateStatus /= 0) stop "Allocation failed. Not enough memory to allocate SN mock covariance."

    sn_mock_covariance = 0._dl

    ! compute the fiducial:
    allocate( sn_fiducial(number), stat = AllocateStatus )
    if (AllocateStatus /= 0) stop "Allocation failed. Not enough memory to allocate SN fiducial."

    call Observed_SN_Magnitude( P, FP, sn_fiducial, sn_mock(1,:), sn_mock(2,:), sn_mock(3,:), &
        & number, err )

    call save_sn_mock_to_file( number, sn_mock, sn_mock_covariance, filename='output', magnitude = sn_fiducial )

    ! all the successful test should return the following.
    stop 100

end program test_fisher_compute_mag
