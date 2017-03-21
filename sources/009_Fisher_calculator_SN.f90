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

!> @file 009_Fisher_calculator_SN.f90
!! This file contains the code that computes the Fisher matrix for Supernovae measurements.

!----------------------------------------------------------------------------------------
!> This module contains the code that computes the Fisher matrix for Supernovae measurements.

!> @author Marco Raveri

module Fisher_calculator_SN

    use precision
    use ModelParams
    use cosmicfish_types
    use omp_lib
    use cosmicfish_utilities
    use cosmicfish_camb_interface

#ifdef COSMICFISH_EFTCAMB
    use EFTinitialization
#endif

    implicit none

    private

    public Observed_SN_Magnitude, Generate_SN_Mock_Data, Generate_SN_Mock_Covariance, save_sn_mock_to_file, &
        & Observed_SN_Magnitude_derivative, Fisher_SN

contains

    ! ---------------------------------------------------------------------------------------------
    !> This subroutine returns the observed magnitude of a set of supernovae based on the input
    !! cosmological parameters.
    subroutine Observed_SN_Magnitude( P, FP, sn_magnitude, redshift, color, stretch, num, err )

        implicit none

        Type(CAMBparams)      , intent(in)  :: P            !< Input CAMBparams
        Type(cosmicfish_params), intent(in)  :: FP           !< Input Cosmicfish params
        real(dl), dimension(:), intent(out) :: sn_magnitude !< Output array with the supernovae magnitude
        real(dl), dimension(:), intent(in)  :: redshift     !< Input array with the values of redshift
        real(dl), dimension(:), intent(in)  :: color        !< Input array with the values of color
        real(dl), dimension(:), intent(in)  :: stretch      !< Input array with the values of stretch
        integer               , intent(in)  :: num          !< number of supernovae.
                                                            !< Length of the sn_magnitude,
                                                            !< redshift, color, stretch vectors
        integer, intent(out)    :: err                      !< Output error code:
                                                            !< 0 = all fine
                                                            !< 1 = error in input
                                                            !< 2 = error in computation
                                                            !< 3 = error in quality
                                                            !< 4 = routine specific

        ! other definitions:
        integer  :: i
        real(dl) :: alpha, beta, zero_mag
        logical  :: EFTsuccess
        ! initialization:
        sn_magnitude = 0._dl
        err          = 0

        call CAMBParams_Set(P)

#ifdef COSMICFISH_EFTCAMB
        ! check stability if EFTCAMB is used
        if (P%EFTflag /= 0) then
            call EFTCAMB_initialization(EFTsuccess)
            if (.not. EFTsuccess) then
                write(*,*) 'EFTCAMB: unstable model'
                err = 4
                return
            end if
        end if
#endif

        alpha    = FP%fisher_SN%alpha_SN
        beta     = FP%fisher_SN%alpha_SN
        zero_mag = FP%fisher_SN%M0_SN

        ! compute the supernova magnitude array:
        do i=1, num
            sn_magnitude(i) = 5._dl*log10( LuminosityDistance(redshift(i)) ) &
                & -alpha*stretch(i) +beta*color(i) + zero_mag -25._dl
        end do

    end subroutine Observed_SN_Magnitude

    ! ---------------------------------------------------------------------------------------------
    !> This subroutine generates a mock SN data set compilation based on the input specifications
    !> contained in the CosmicFish parameters.
    subroutine Generate_SN_Mock_Data( P, FP, sn_mock )

        use Random_utilities

        implicit none

        Type(CAMBparams)        , intent(in)  :: P            !< Input CAMBparams
        Type(cosmicfish_params)  , intent(in)  :: FP           !< Input Cosmicfish params
        real(dl), dimension(:,:), intent(out) :: sn_mock      !< Output supernovae mock
                                                              !< This array has to be allocated
                                                              !< with dimension 3 x (FP%fisher_SN%total_SN_number)

        ! other definitions:
        integer  :: i, j, k

        ! initialize the random number generator:
        Rand_Feedback = 0
        call initRandom()

        ! generate the mock data:
        i = 1
        do j=1, FP%fisher_SN%number_SN_windows
            do k=1, FP%fisher_SN%SN_number(j)
                sn_mock(1,i) = random_uniform( FP%fisher_SN%SN_redshift_start(j), FP%fisher_SN%SN_redshift_end(j)) ! mock redshift
                sn_mock(2,i) = random_gaussian( 0._dl, FP%fisher_SN%color_dispersion )                             ! mock color
                sn_mock(3,i) = random_gaussian( 0._dl, FP%fisher_SN%stretch_dispersion )                           ! mock stretch
                i = i + 1
            end do
        end do

    end subroutine Generate_SN_Mock_Data

    ! ---------------------------------------------------------------------------------------------
    !> This subroutine generates a mock SN covariance based on the input specifications
    !> contained in the CosmicFish parameters and the SN mock.
    subroutine Generate_SN_Mock_Covariance( P, FP, sn_mock, sn_mock_covariance )

        use constants, only : c

        implicit none

        Type(CAMBparams)        , intent(in)  :: P            !< Input CAMBparams
        Type(cosmicfish_params)  , intent(in)  :: FP           !< Input Cosmicfish params
        real(dl), dimension(:,:), intent(in)  :: sn_mock      !< Input supernovae mock
                                                              !< This array has to be allocated
                                                              !< with dimension 3 x (FP%fisher_SN%total_SN_number)
        real(dl), dimension(:)  , intent(out) :: sn_mock_covariance !< SN mock covariance array (assumed diagonal)
                                                                    !< of size FP%fisher_SN%total_SN_number

        ! other definitions:
        integer  :: i, j, k

        ! generate the mock covariance first terms:
        do i=1, FP%fisher_SN%total_SN_number
            sn_mock_covariance(i) = (FP%fisher_SN%magnitude_sigma)**2 &
                & +(FP%fisher_SN%sigma_lens_0*sn_mock(1,i))**2 & ! lensing term
                & +(5._dl/sn_mock(1,i)/log(10._dl)*FP%fisher_SN%c_sigmaz/c*1000._dl)**2 ! redshift term
        end do

        ! add the systematic effects to the covariance:
        do i=1, FP%fisher_SN%total_SN_number
            sn_mock_covariance(i) = sn_mock_covariance(i) &
                & + FP%fisher_SN%beta_SN**2*(FP%fisher_SN%dcolor_offset+FP%fisher_SN%dcolor_zcorr*sn_mock(1,i)**2)**2  &
                & + FP%fisher_SN%alpha_SN**2*(FP%fisher_SN%dshape_offset+FP%fisher_SN%dshape_zcorr*sn_mock(1,i)**2)**2  &
                & + 2._dl*FP%fisher_SN%alpha_SN*(FP%fisher_SN%cov_ms_offset+FP%fisher_SN%cov_ms_zcorr*sn_mock(1,i)**2)     &
                & - 2._dl*FP%fisher_SN%beta_SN*(FP%fisher_SN%cov_mc_offset+FP%fisher_SN%cov_mc_zcorr*sn_mock(1,i)**2)     &
                & - 2._dl*FP%fisher_SN%alpha_SN*FP%fisher_SN%beta_SN*(FP%fisher_SN%cov_sc_offset+FP%fisher_SN%cov_sc_zcorr*sn_mock(1,i)**2)
        end do

    end subroutine Generate_SN_Mock_Covariance

    ! ---------------------------------------------------------------------------------------------
    !> This subroutine outputs to file a SN covariance.
    subroutine Save_Sn_Mock_To_File( mock_size, sn_mock_data, sn_mock_covariance, filename, magnitude )

        implicit none

        integer                 , intent(in) :: mock_size           !< Total length of the SN mock data set
        real(dl), dimension(:,:), intent(in) :: sn_mock_data        !< SN mock data, of dimension 3 x mock_size
        real(dl), dimension(:)  , intent(in) :: sn_mock_covariance  !< SN mock covariance, diagonal, of dimension mock_size
        character(len=*)        , intent(in) :: filename            !< Name of the file where the SN mock will be saved

        real(dl), dimension(:)  , intent(in), optional :: magnitude !< Optionally the array with the magnitude of the SN catalog
                                                                    !< can be written to file

        ! other definitions:
        integer  :: i

        ! write the mock:
        open(unit=666, FILE=filename, ACTION="write", STATUS="replace")

        if ( present(magnitude) ) then

            write(666,'(a)') "# number of the SN, redshift, magnitude, color, stretch, covariance "
            do i = 1, mock_size
                write(666,'(I10,5E15.5)') i, sn_mock_data(1,i), magnitude(i), sn_mock_data(2,i), sn_mock_data(3,i), sn_mock_covariance(i)
            end do

        else

            write(666,'(a)') "# number of the SN, redshift, color, stretch, covariance "
            do i = 1, mock_size
                write(666,'(I10,4E15.5)') i, sn_mock_data(1,i), sn_mock_data(2,i), sn_mock_data(3,i), sn_mock_covariance(i)
            end do

        end if

        close(666)

    end subroutine save_sn_mock_to_file

    ! ---------------------------------------------------------------------------------------------
    !> This subroutine computes the derivative of the observed magnitude of a set of supernovae.
    !! The calculation is carried out with a fixed stepsize finite difference formula or by
    !! the Richardson’s deferred approach to the limit.
    subroutine Observed_SN_Magnitude_derivative( P, FP, param_num, &
        & mag_initial, mag_derivative, mag_der_error, &
        & sn_mock, &
        & initial_step, &
        & num_param, &
        & error_code  )

        implicit none

        Type(CAMBparams)       , intent(in) :: P       !< Input CAMBparams
        Type(cosmicfish_params) , intent(in) :: FP      !< Input Cosmicfish params

        integer , intent(in)  :: param_num  !< Input number of the parameter wrt the derivative is computed
        integer , intent(in)  :: num_param  !< Input number of parameters involved in the calculation

        real(dl), dimension(:), intent(in)  :: mag_initial     !< SN magnitude computed for h=0
        real(dl), dimension(:), intent(out) :: mag_derivative  !< Output returns the derivative of the SN magnitude vector
        real(dl), dimension(:), intent(out) :: mag_der_error   !< Output returns the error on the megnitude derivative

        real(dl), dimension(:,:), intent(in) :: sn_mock        !< Input supernovae mock. Three columns with SN redshift, color and stretch

        real(dl), intent(in)  :: initial_step  !< Input initial value of h. Does not need to be small for the adaptive algorithm.
        integer , intent(out) :: error_code    !< Output error code. When run the code will not stop on error but report it
                                               !<  0 = everything was fine computation exited correctly
                                               !<  1 = error on input parameters
                                               !<  2 = error on computation
                                               !<  3 = error on quality of the results

        ! definitions:
        real(dl), dimension(:), allocatable :: sn_temp_1, sn_temp_2, errt
        integer , dimension(:), allocatable :: accept_mask
        real(dl), dimension(:,:,:), allocatable :: a
        real(dl), dimension(num_param) :: param_vector, fiducial_param_vector, temp_param_vector
        ! parameters of the algorithm:
        real(dl), parameter :: con  = 1.4_dl  ! parameter that decides the decrease in the stepsize h
        real(dl), parameter :: big  = 1.d+30  ! a big number used to fill the error array initially
        real(dl), parameter :: safe = 2._dl   ! parameter used to define the stopping cryterium
        integer , parameter :: ntab = 30      ! maximum number of Cls computation
        real(dl) :: con2 = con*con
        integer  :: i, j, k, l, m, side, err, AllocateStatus
        real(dl) :: fac, hh
        real(dl) :: temp1, temp2, temp3
        Type(CAMBparams) :: P_temp
        Type(cosmicfish_params) :: FP_temp
        integer  :: accepted
        real(dl) :: acceptance
        integer  :: mock_sn_size
        real(dl), dimension(:), allocatable :: redshift, color, stretch

        ! initialize output variables:
        error_code     = 0
        mag_derivative = 0._dl
        mag_der_error  = 0._dl

        ! check input values are correct
        if ( initial_step == 0._dl ) then
            stop 'Observed_SN_Magnitude_derivative initial stepsize is zero'
        end if

        mock_sn_size = FP%fisher_SN%total_SN_number

        ! allocate the arrays:
        if ( FP%adaptivity ) then
            allocate( sn_temp_1(mock_sn_size),  &
                & sn_temp_2(mock_sn_size),  &
                & errt(mock_sn_size),       &
                & accept_mask(mock_sn_size),&
                & stat = AllocateStatus )
            if (AllocateStatus /= 0) stop "Observed_SN_Magnitude_derivative: Allocation failed. Not enough memory"
            allocate( a(ntab, ntab, mock_sn_size) , stat = AllocateStatus )
            if (AllocateStatus /= 0) stop "Observed_SN_Magnitude_derivative: Allocation failed. Not enough memory for the derivative array"
        else
            allocate( sn_temp_1(mock_sn_size),  &
                & sn_temp_2(mock_sn_size),  &
                & stat = AllocateStatus )
            if (AllocateStatus /= 0) stop "Observed_SN_Magnitude_derivative: Allocation failed. Not enough memory"
            allocate( a(1,1, mock_sn_size) , stat = AllocateStatus )
            if (AllocateStatus /= 0) stop "Observed_SN_Magnitude_derivative: Allocation failed. Not enough memory for the derivative array"
        end if

        allocate( redshift(mock_sn_size),       &
            & color(mock_sn_size),       &
            & stretch(mock_sn_size),     &
            & stat = AllocateStatus )

        ! copy the data from the sn mock:
        redshift(:) = sn_mock(1,:)
        color(:)    = sn_mock(2,:)
        stretch(:)  = sn_mock(3,:)

        ! initialize the parameter array:
        P_temp  = P
        FP_temp = FP

        call CAMB_params_to_params_array( P_temp, FP, num_param, param_vector )

        fiducial_param_vector = param_vector
        temp_param_vector     = param_vector

        ! initialize the computation:
        hh            = initial_step
        mag_der_error = big
        side          = 0 ! both sides
        if ( FP%adaptivity ) accept_mask = 0

        ! do the initial step on the right:
        param_vector(param_num)      = fiducial_param_vector(param_num) + hh
        call params_array_to_CAMB_param( P_temp, FP_temp, num_param, param_vector )
        call Observed_SN_Magnitude( P_temp, FP_temp, sn_temp_1, redshift, color, stretch, mock_sn_size, err )
        ! if on the right we cannot go we go to the left:
        if ( err /= 0 ) side = 1
        ! do the initial step on the left:
        temp_param_vector(param_num) = fiducial_param_vector(param_num) - hh
        call params_array_to_CAMB_param( P_temp, FP_temp, num_param, temp_param_vector )
        call Observed_SN_Magnitude( P_temp, FP_temp, sn_temp_2, redshift, color, stretch, mock_sn_size, err )

        ! if on the left we cannot go decide what to do:
        if ( err /= 0 ) then
            if ( side == 1 ) then ! we cannot go to the left nor to the right. Return derivative zero with error code.
                error_code      = 4
                mag_der_error   = 0.0_dl
                mag_derivative   = 0.0_dl
                return
            else if ( side == 0 ) then ! we cannot go to the left but we can go to the right:
                side = 2
            end if
        end if

        ! store the results based on the side we decided to go:
        if (side == 0) then ! both sides
            a(1,1,:) = (sn_temp_1(:) - sn_temp_2(:))/(2.0_dl*hh)
        else if (side == 1) then ! left side: increment negative
            a(1,1,:) = (sn_temp_2(:) - mag_initial(:))/(-hh)
        else if (side == 2) then ! right side: increment positive
            a(1,1,:) = (sn_temp_1(:) - mag_initial(:))/(hh)
        end if

        ! fixed stepsize derivative if adaptivity is not wanted:
        if ( .not. FP%adaptivity ) then
            mag_derivative(:) = a(1,1,:)
            mag_der_error     = 0.0_dl
            deallocate( sn_temp_1, sn_temp_2, a )
            return
        end if

        ! start the derivative computation:
        do i = 2, ntab

            ! feedback for the adaptive derivative:
            if ( FP%cosmicfish_feedback >= 2 ) then
                accepted   = sum(accept_mask)
                acceptance = REAL(accepted)/REAL(mock_sn_size)
                print*, 'iteration:', i
                print*, '  stepsize:', hh
                print*, '  accepted elements', accepted, 'accepted ratio', acceptance
            end if

            ! decrease the stepsize:
            hh = hh/con

            ! compute the finite difference:
            if (side == 0) then
                param_vector(param_num)      = fiducial_param_vector(param_num) + hh
                temp_param_vector(param_num) = fiducial_param_vector(param_num) - hh
                ! compute for + hh
                call params_array_to_CAMB_param( P_temp, FP_temp, num_param, param_vector )
                call Observed_SN_Magnitude( P_temp, FP_temp, sn_temp_1, redshift, color, stretch, mock_sn_size, err )
                if ( err /= 0 ) then
                    error_code     = 4
                    mag_der_error  = 0.0_dl
                    mag_derivative = 0.0_dl
                    return
                end if
                ! compute for - hh
                call params_array_to_CAMB_param( P_temp, FP_temp, num_param, temp_param_vector )
                call Observed_SN_Magnitude( P_temp, FP_temp, sn_temp_2, redshift, color, stretch, mock_sn_size, err )
                if ( err /= 0 ) then
                    error_code     = 4
                    mag_der_error  = 0.0_dl
                    mag_derivative = 0.0_dl
                    return
                end if
                ! store the result
                a(1,i,:) = (sn_temp_1(:) - sn_temp_2(:))/(2.0_dl*hh)
            else if (side==1) then
                param_vector(param_num)      = fiducial_param_vector(param_num) - hh
                call params_array_to_CAMB_param( P_temp, FP_temp, num_param, param_vector )
                call Observed_SN_Magnitude( P_temp, FP_temp, sn_temp_1, redshift, color, stretch, mock_sn_size, err )
                if ( err /= 0 ) then
                    error_code     = 4
                    mag_der_error  = 0.0_dl
                    mag_derivative = 0.0_dl
                    return
                end if
                ! store the result
                a(1,i,:) = (sn_temp_1(:) - mag_initial(:))/(-hh)
            else if (side==2) then
                param_vector(param_num)      = fiducial_param_vector(param_num) + hh
                call params_array_to_CAMB_param( P_temp, FP_temp, num_param, param_vector )
                call Observed_SN_Magnitude( P_temp, FP_temp, sn_temp_1, redshift, color, stretch, mock_sn_size, err )
                if ( err /= 0 ) then
                    error_code     = 4
                    mag_der_error  = 0.0_dl
                    mag_derivative = 0.0_dl
                    return
                end if
                ! store the result
                a(1,i,:) = (sn_temp_1(:) - mag_initial(:))/(hh)
            end if

            ! now do the extrapolation table:
            fac=CON2
            do j = 2, i

                forall ( k = 1:mock_sn_size, accept_mask(k) == 0 )
                    ! compute the interpolation table:
                    a(j,i,k) = (a(j-1,i,k)*fac-a(j-1,i-1,k))/(fac-1._dl)
                    ! estimate the error:
                    errt(k)  = max( abs(a(j,i,k)-a(j-1,i,k)), &
                        & abs(a(j,i,k)-a(j-1,i-1,k)))
                end forall

                fac=CON2*fac

                ! if there is an improvement save the result:

                forall ( k = 1:mock_sn_size, &
                    & errt(k) .le. mag_der_error(k) .and. accept_mask(k) == 0 )

                    mag_der_error(k) = errt(k)
                    mag_derivative(k) = a(j,i,k)
                end forall

            end do
            ! accept a value of the derivative if the value of the extrapolation is stable
            ! Store this information into the mask matrix
            ! so that the value contained in cl_error and cl_derivative is not overwritten.
            forall ( k = 1:mock_sn_size, &
                & accept_mask(k) == 0 .and. &
                & (a(i,i,k) - a(i-1,i-1,k)) .ge. safe*mag_der_error(k) )

                accept_mask(k) = 1
            end forall

            ! if all the entries match the required accuracy exit the loop
            if ( sum(accept_mask) == mock_sn_size ) exit

        end do

        ! quality check on the results:
        if ( sum(accept_mask) /= mock_sn_size ) then
            error_code = 3
        end if

        ! deallocate the arrays:
        if ( FP%adaptivity ) then
            deallocate( sn_temp_1, sn_temp_2, errt, accept_mask, a )
        end if

        return

    end subroutine Observed_SN_Magnitude_derivative

    ! ---------------------------------------------------------------------------------------------
    !> This subroutine computes the SN Fisher matrix.
    subroutine Fisher_SN( P, FP, num_param, Fisher_Matrix, outroot )

        implicit none

        Type(CAMBparams)       , intent(in) :: P       !< Input CAMBparams
        Type(cosmicfish_params) , intent(in) :: FP      !< Input Cosmicfish params

        integer , intent(in)  :: num_param                                        !< Input dimension of the Fisher matrix
        real(dl), dimension(num_param,num_param), intent(out) :: Fisher_Matrix    !< Output Fisher matrix
        character(len=*), intent(in), optional :: outroot                         !< Optional input: filename that will be used to dump to file optional output if feedback is greater than 1.

        real(dl) :: time_1, time_2, time_1_MC, time_2_MC
        real(dl) :: initial_step, temp, temp1, temp2
        integer  :: AllocateStatus, err, ind, ind2
        integer  :: i, j, k, l, number
        real(dl), allocatable :: param_array(:)
        real(dl), allocatable, dimension(:,:) :: sn_mock
        real(dl), allocatable, dimension(:)   :: sn_mock_covariance
        real(dl), allocatable, dimension(:)   :: sn_fiducial, sn_derivative, sn_temp
        real(dl), allocatable, dimension(:,:) :: sn_derivative_array

        integer  :: MC_i
        real(dl), dimension(num_param,num_param) :: MC_Fisher_Matrix

        ! allocate a parameter array:
        allocate( param_array(num_param) )
        ! write the parameters in the parameter array:
        call CAMB_params_to_params_array( P, FP, num_param, param_array )

        number = FP%fisher_SN%total_SN_number

        ! allocation:
        allocate( sn_mock(3,number), stat = AllocateStatus )
        if (AllocateStatus /= 0) stop "Allocation failed. Not enough memory to allocate SN mock."
        allocate( sn_mock_covariance(number), stat = AllocateStatus )
        if (AllocateStatus /= 0) stop "Allocation failed. Not enough memory to allocate SN mock covariance."
        allocate( sn_fiducial(number), sn_derivative(number), sn_temp(number), stat = AllocateStatus )
        if (AllocateStatus /= 0) stop "Allocation failed. Not enough memory to allocate SN fiducial."
        allocate( sn_derivative_array( num_param, number ), stat = AllocateStatus)
        if (AllocateStatus /= 0) stop "Allocation failed. Not enough memory to allocate sn_derivative_array"

        ! do the Monte Carlo:

        MC_Fisher_Matrix = 0._dl

        time_1_MC = omp_get_wtime()
        do MC_i=1, FP%fisher_SN%SN_Fisher_MC_samples

            ! generate the SN mock:
            call Generate_SN_Mock_Data( P, FP, sn_mock )

            ! generate the mock covariance:
            call Generate_SN_Mock_Covariance( P, FP, sn_mock, sn_mock_covariance )

            ! compute the fiducial:
            time_1 = omp_get_wtime()
            call Observed_SN_Magnitude( P, FP, sn_fiducial, sn_mock(1,:), sn_mock(2,:), sn_mock(3,:), &
                & number, err )
            time_2 = omp_get_wtime() - time_1

            if ( present(outroot) ) call save_sn_mock_to_file( number, sn_mock, sn_mock_covariance, &
                & filename=TRIM(outroot)//'sn_fiducial_'//trim(adjustl(integer_to_string(MC_i)))//'.dat', magnitude = sn_fiducial )

            if ( err == 0 ) then
            else if ( err == 4 ) then
                write(*,*) 'Fiducial model unstable. Calculation cannot proceed.'
                stop
            else
                write(*,*) 'Error in computing the fiducial. Calculation cannot proceed.'
                stop
            end if

            ! print some feedback
            if ( MC_i==1 .and. FP%cosmicfish_feedback >= 1 ) then
                write(*,'(a)')         '**************************************************************'
                write(*,'(a,i10)')         'Number of supernovae         : ', number
                write(*,'(a,i10)')         'Number of parameters         : ', num_param
                write(*,'(a,i10)')         'Number of MC repetitions     : ', FP%fisher_SN%SN_Fisher_MC_samples
                write(*,'(a,f10.3,a)')     'Time taken for SN calculation: ', time_2, ' (s)'
                write(*,'(a,i10)')         'Number of omp processes      : ', OMP_GET_MAX_THREADS()
                if ( FP%adaptivity ) then
                    write(*,'(a,f10.3,a)') 'Estimated completion time    : ', FP%fisher_SN%SN_Fisher_MC_samples*time_2*20._dl*num_param, ' (s)'
                else
                    write(*,'(a,f10.3,a)') 'Estimated completion time    : ', FP%fisher_SN%SN_Fisher_MC_samples*time_2*2._dl*num_param, ' (s)'
                end if
                write(*,'(a)')         '**************************************************************'
            end if

            if ( MC_i>1 .and. FP%cosmicfish_feedback >= 1 ) then
                write(*,'(a)') '**************************************************************'
                write(*,'(a,i10)') 'Doing MC computation of SN Fisher sample: ', MC_i
                write(*,'(a)') '**************************************************************'
            end if

            ! compute the derivatives of the Cls:
            if ( FP%cosmicfish_feedback >= 1 ) then
                write(*,'(a)') 'Computation of the SN magnitude derivatives'
            end if

            time_1 = omp_get_wtime()
            do ind = 1, num_param

                if ( FP%cosmicfish_feedback >= 1 ) then
                    write(*,'(a,i5,a,E15.5)') 'Doing parameter number: ', ind, ' value: ', param_array(ind)
                end if

                if ( param_array(ind) /= 0._dl ) then
                    initial_step = param_array(ind)*3.0_dl/100._dl
                else
                    initial_step = 3.0_dl/100._dl
                end if

                call Observed_SN_Magnitude_derivative( P, FP, ind, &
                    & sn_fiducial, sn_derivative, sn_temp, &
                    & sn_mock, &
                    & initial_step, &
                    & num_param, &
                    & err  )

                if (err /= 0) print*, 'WARNING: something went wrong with the SN derivative of parameter:', ind, 'Error code:', err

                if ( MC_i==1 .and. present(outroot) .and. FP%cosmicfish_feedback >= 2 ) then
                    call save_sn_mock_to_file( number, sn_mock, sn_mock_covariance, &
                        & filename=TRIM(outroot)//'sn_derivative_'//trim(adjustl(integer_to_string(ind)))//'.dat', magnitude = sn_derivative )
                end if

                sn_derivative_array(ind,:) = sn_derivative(:)

            end do
            time_2 = omp_get_wtime() - time_1

            if ( FP%cosmicfish_feedback >= 1 ) then
                write(*,*)
                write(*,'(a,f8.3,a)') 'Time to compute all the derivatives:', time_2, ' (s)'
            end if

            ! compute the Fisher:
            do ind = 1, num_param
                do ind2 = 1, num_param

                    temp2 = 0._dl

                    do k = 1, number
                        temp2 = temp2 + sn_derivative_array(ind,k)*sn_derivative_array(ind2,k)/sn_mock_covariance(k)
                    end do

                    MC_Fisher_Matrix(ind, ind2) = MC_Fisher_Matrix(ind, ind2) + temp2

                end do
            end do
            time_2 = omp_get_wtime() - time_1

            if ( FP%cosmicfish_feedback >= 1 ) then
                write(*,'(a,f8.3,a)') 'Time to compute the Fisher matrix:', time_2, ' (s)'
            end if

        end do
        time_2_MC = omp_get_wtime() - time_1

        Fisher_Matrix = MC_Fisher_Matrix/FP%fisher_SN%SN_Fisher_MC_samples

        if ( FP%cosmicfish_feedback >= 1 ) then
            write(*,'(a)') '**************************************************************'
            write(*,'(a,f8.3,a)') 'Time to compute the MC Fisher matrix:', time_2_MC, ' (s)'
            write(*,'(a)') '**************************************************************'
        end if

    end subroutine Fisher_SN

end module Fisher_calculator_SN

