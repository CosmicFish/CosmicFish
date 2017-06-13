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

!> @file 013_Fisher_calculator_varalpha.f90
!! This file contains the code that computes the Fisher matrix for alpha
!! variation.
!! Reference paper: http://arxiv.org/abs/1311.5841

!----------------------------------------------------------------------------------------
!> This module contains the code that computes the Fisher matrix for
!> alpha variation (in redshift)

!> @author Matteo Martinelli, Marco Raveri

module Fisher_calculator_alpha

    use precision
    use ModelParams
    use cosmicfish_types
    use CAMBmain
    use CAMB
    use constants
    use cosmicfish_camb_interface
    use lensnoise
    use omp_lib
    use cosmicfish_utilities
    use Matrix_utilities

#ifdef COSMICFISH_EFTCAMB
    use EFTinitialization
#endif

    implicit none

    private

    public alpha_theo_Calc, save_alpha_mock_to_file, alpha_DerivativeCalc, Fisher_alpha

contains

    ! ---------------------------------------------------------------------------------------------
    !> This subroutine returns the theoretical prediction of redshift variation of alpha given the input
    !! cosmological parameters and the observation interval given in specifications.
    subroutine alpha_theo_Calc( P, FP, varalpha, err )
        implicit none

        Type(CAMBparams)      , intent(in)                                        :: P            !< Input CAMBparams
        Type(cosmicfish_params), intent(in)                                       :: FP           !< Input Cosmicfish params
        real(dl), dimension(FP%fisher_alpha%number_alpha_redshifts), intent(out)  :: varalpha     !< Output array with the
                                                                                                  !< theoretical alpha variation
        integer, intent(out)                                                      :: err          !< Output error code:
                                                                                                  !< 0 = all fine
                                                                                                  !< 1 = error in input
                                                                                                  !< 2 = error in computation
                                                                                                  !< 3 = error in quality
                                                                                                  !< 4 = routine specific

        ! other definitions:
        integer  :: i
        logical  :: EFTsuccess
!        real(dl), parameter :: conv_time=3.15e7, conv_mpc=1/(3.09e19) !years to seconds and 1/Mpc to 1/km
!        real(dl) ::deltaT
        ! initialization:
        varalpha      = 0._dl
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




        ! compute the supernova magnitude array:
        do i=1, FP%fisher_alpha%number_alpha_redshifts
            varalpha(i) = FP%fisher_alpha%alpha_coupling !put equation 9 (generalized) here
            write(*,*) FP%fisher_alpha%alpha_redshift(i), varalpha(i), FP%fisher_alpha%alpha_error(i)
        end do

    end subroutine alpha_theo_Calc


    ! ---------------------------------------------------------------------------------------------
    !> This subroutine outputs to file a SN covariance.
    subroutine save_alpha_mock_to_file(FP, mock_size, alpha_theory, alpha_error, filename)
        implicit none

        Type(cosmicfish_params), intent(in)        :: FP             !< Input Cosmicfish params
        integer, intent(in)                        :: mock_size      !< number of mock points
        real(dl), dimension(mock_size), intent(in) :: alpha_theory   !< fiducial cosmology redshift drift
        real(dl), dimension(mock_size), intent(in) :: alpha_error    !< error on delta v
        character(len=*), intent(in)               :: filename       !< Name of the file where the alpha variation mock will be saved
        integer :: i

        open(unit=42, FILE=filename, ACTION="write", STATUS="replace")

        write(42,'(a)') "# redshift, Delta alpha/alpha, sigma(Delta alpha/alpha) "
        do i = 1, mock_size
            write(42,'(5E15.5)') FP%fisher_alpha%alpha_redshift(i), alpha_theory(i), alpha_error(i)
        end do

        close(42)



    end subroutine save_alpha_mock_to_file



    ! ---------------------------------------------------------------------------------------------
    !> This subroutine computes the derivative of the observed alpha variation signal.
    !! The calculation is carried out with a fixed stepsize finite difference formula or by
    !! the Richardson’s deferred approach to the limit.
    subroutine alpha_DerivativeCalc( P, FP, param_num, alpha_derivative, alpha_der_error, alpha_initial, initial_step, num_param, error_code)
        implicit none

        Type(CAMBparams)       , intent(in) :: P       !< Input CAMBparams
        Type(cosmicfish_params) , intent(in) :: FP      !< Input Cosmicfish params

        integer , intent(in)  :: param_num  !< Input number of the parameter wrt the derivative is computed
        integer , intent(in)  :: num_param  !< Input number of parameters involved in the calculation

        real(dl), dimension(:), intent(in)  :: alpha_initial     !< deltav computed for h=0
        real(dl), dimension(:), intent(out) :: alpha_derivative  !< Output returns the derivative of the deltav vector
        real(dl), dimension(:), intent(out) :: alpha_der_error   !< Output returns the error on the deltav derivative


        real(dl), intent(in)  :: initial_step  !< Input initial value of h. Does not need to be small for the adaptive algorithm.
        integer , intent(out) :: error_code    !< Output error code. When run the code will not stop on error but report it
                                               !<  0 = everything was fine computation exited correctly
                                               !<  1 = error on input parameters
                                               !<  2 = error on computation
                                               !<  3 = error on quality of the results


        integer :: mock_alpha_size
        real(dl), dimension(:), allocatable :: alpha_temp_1, alpha_temp_2, errt
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
        real(dl), dimension(:), allocatable :: redshift, varalpha



        ! initialize output variables:
        error_code     = 0
        alpha_derivative = 0._dl
        alpha_der_error  = 0._dl


        ! check input values are correct
        if ( initial_step == 0._dl ) then
            stop 'alpha_derivative initial stepsize is zero'
        end if

        mock_alpha_size = FP%fisher_alpha%number_alpha_redshifts

        ! allocate the arrays:
        if ( FP%adaptivity ) then
           write(*,*) 'Adaptivity not yet implemented with alpha'
        !THIS IS FOR MARCO TO ENJOY
        !            allocate( sn_temp_1(mock_sn_size),  &
        !                & sn_temp_2(mock_sn_size),  &
        !                & errt(mock_sn_size),       &
        !                & accept_mask(mock_sn_size),&
        !                & stat = AllocateStatus )
        !            if (AllocateStatus /= 0) stop "Observed_SN_Magnitude_derivative: Allocation failed. Not enough memory"
        !            allocate( a(ntab, ntab, mock_sn_size) , stat = AllocateStatus )
        !            if (AllocateStatus /= 0) stop "Observed_SN_Magnitude_derivative: Allocation failed. Not enough memory for the derivative array"
        else
            allocate( alpha_temp_1(mock_alpha_size),  &
                & alpha_temp_2(mock_alpha_size),  &
                & stat = AllocateStatus )
            if (AllocateStatus /= 0) stop "alpha_derivative: Allocation failed. Not enough memory"
            allocate( a(1,1, mock_alpha_size) , stat = AllocateStatus )
            if (AllocateStatus /= 0) stop "alpha_derivative: Allocation failed. Not enough memory for the derivative array"
        end if

        allocate( redshift(mock_alpha_size),       &
            & varalpha(mock_alpha_size),       &
            & stat = AllocateStatus )

        ! copy the data from the sn mock:
        redshift(:) = FP%fisher_alpha%alpha_redshift(:)!sn_mock(1,:)
        !        deltav(:)    = sn_mock(2,:)

        ! initialize the parameter array:
        P_temp  = P
        FP_temp = FP

        call CAMB_params_to_params_array( P_temp, FP, num_param, param_vector )

        fiducial_param_vector = param_vector
        temp_param_vector     = param_vector

        ! initialize the computation:
        hh            = initial_step
        alpha_der_error = big
        side          = 0 ! both sides
        if ( FP%adaptivity ) accept_mask = 0

        ! do the initial step on the right:
        param_vector(param_num)      = fiducial_param_vector(param_num) + hh
        call params_array_to_CAMB_param( P_temp, FP_temp, num_param, param_vector )
        call alpha_Theo_Calc(P_temp, FP_temp, alpha_temp_1, err)
        ! if on the right we cannot go we go to the left:
        if ( err /= 0 ) side = 1
        ! do the initial step on the left:
        temp_param_vector(param_num) = fiducial_param_vector(param_num) - hh
        call params_array_to_CAMB_param( P_temp, FP_temp, num_param, temp_param_vector )
        call alpha_Theo_Calc(P_temp, FP_temp, alpha_temp_2, err)

        ! if on the left we cannot go decide what to do:
        if ( err /= 0 ) then
            if ( side == 1 ) then ! we cannot go to the left nor to the right. Return derivative zero with error code.
                error_code      = 4
                alpha_der_error   = 0.0_dl
                alpha_derivative   = 0.0_dl
                return
            else if ( side == 0 ) then ! we cannot go to the left but we can go to the right:
                side = 2
            end if
        end if

        ! store the results based on the side we decided to go:
        if (side == 0) then ! both sides
            a(1,1,:) = (alpha_temp_1(:) - alpha_temp_2(:))/(2.0_dl*hh)
        else if (side == 1) then ! left side: increment negative
            a(1,1,:) = (alpha_temp_2(:) - alpha_initial(:))/(-hh)
        else if (side == 2) then ! right side: increment positive
            a(1,1,:) = (alpha_temp_1(:) - alpha_initial(:))/(hh)  
        end if

        ! fixed stepsize derivative if adaptivity is not wanted:
        if ( .not. FP%adaptivity ) then
            alpha_derivative(:) = a(1,1,:)
            alpha_der_error     = 0.0_dl
            deallocate( alpha_temp_1, alpha_temp_2, a )
            return
        end if

        ! start the derivative computation:
        do i = 2, ntab

            ! feedback for the adaptive derivative:
            if ( FP%cosmicfish_feedback >= 2 ) then
                accepted   = sum(accept_mask)
                acceptance = REAL(accepted)/REAL(mock_alpha_size)
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
                call alpha_Theo_Calc(P_temp, FP_temp, alpha_temp_1, err)
                if ( err /= 0 ) then
                    error_code     = 4
                    alpha_der_error  = 0.0_dl
                    alpha_derivative = 0.0_dl
                    return
                end if
                ! compute for - hh
                call params_array_to_CAMB_param( P_temp, FP_temp, num_param, temp_param_vector )
                call alpha_Theo_Calc(P_temp, FP_temp, alpha_temp_2, err)
                if ( err /= 0 ) then
                    error_code     = 4
                    alpha_der_error  = 0.0_dl
                    alpha_derivative = 0.0_dl
                    return
                end if
                ! store the result
                a(1,i,:) = (alpha_temp_1(:) - alpha_temp_2(:))/(2.0_dl*hh)
            else if (side==1) then
                param_vector(param_num)      = fiducial_param_vector(param_num) - hh
                call params_array_to_CAMB_param( P_temp, FP_temp, num_param, param_vector )
                call alpha_Theo_Calc(P_temp, FP_temp, alpha_temp_1, err)
                if ( err /= 0 ) then
                    error_code     = 4
                    alpha_der_error  = 0.0_dl
                    alpha_derivative = 0.0_dl
                    return
                end if
                ! store the result
                a(1,i,:) = (alpha_temp_1(:) - alpha_initial(:))/(-hh)
            else if (side==2) then
                param_vector(param_num)      = fiducial_param_vector(param_num) + hh
                call params_array_to_CAMB_param( P_temp, FP_temp, num_param, param_vector )
                call alpha_Theo_Calc(P_temp, FP_temp, alpha_temp_1, err)
                if ( err /= 0 ) then
                    error_code     = 4
                    alpha_der_error  = 0.0_dl
                    alpha_derivative = 0.0_dl
                    return
                end if
                ! store the result
                a(1,i,:) = (alpha_temp_1(:) - alpha_initial(:))/(hh)
            end if

            ! now do the extrapolation table:
            fac=CON2
            do j = 2, i

                forall ( k = 1:mock_alpha_size, accept_mask(k) == 0 )
                    ! compute the interpolation table:
                    a(j,i,k) = (a(j-1,i,k)*fac-a(j-1,i-1,k))/(fac-1._dl)
                    ! estimate the error:
                    errt(k)  = max( abs(a(j,i,k)-a(j-1,i,k)), &
                        & abs(a(j,i,k)-a(j-1,i-1,k)))
                end forall

                fac=CON2*fac

                ! if there is an improvement save the result:

                forall ( k = 1:mock_alpha_size, &
                    & errt(k) .le. alpha_der_error(k) .and. accept_mask(k) == 0 )

                    alpha_der_error(k) = errt(k)
                    alpha_derivative(k) = a(j,i,k)
                end forall

                ! accept a value of the derivative if the value of the extrapolation is stable
                ! Store this information into the mask matrix
                ! so that the value contained in rd_error and rd_derivative is not overwritten.
                forall ( k = 1:mock_alpha_size, &
                    & accept_mask(k) == 0 .and. &
                    & (a(i,i,k) - a(i-1,i-1,k)) .ge. safe*alpha_der_error(k) )

                    accept_mask(k) = 1
                end forall

                ! if all the entries match the required accuracy exit the loop
                if ( sum(accept_mask) == mock_alpha_size ) exit

            end do

            ! quality check on the results:
            if ( sum(accept_mask) /= mock_alpha_size ) then
                error_code = 3
            end if

            ! deallocate the arrays:
            if ( FP%adaptivity ) then
            !MARCO!
            !                deallocate( sn_temp_1, sn_temp_2, errt, accept_mask, a )
            end if

            return

        end do

    end subroutine alpha_DerivativeCalc


    ! ---------------------------------------------------------------------------------------------
    !> This subroutine computes the Fisher matrix for alpha variations.
    subroutine Fisher_alpha( P, FP, num_param, Fisher_Matrix, outroot )
        implicit none

        Type(CAMBparams)      , intent(in)  :: P                    !< Input CAMBparams
        Type(cosmicfish_params), intent(in)  :: FP                   !< Input Cosmicfish params
        real(dl), allocatable, dimension(:) :: varalpha,alphasigma      !< Output array with the theoretical alpha variation

        integer , intent(in)  :: num_param                                        !< Input dimension of the Fisher matrix
        real(dl), dimension(num_param,num_param), intent(out) :: Fisher_Matrix    !< Output Fisher matrix
        character(len=*), intent(in), optional :: outroot                         !< Optional input: filename that will be used to dump to file optional output if feedback is greater than 1.

        real(dl) :: time_1, time_2
        real(dl) :: initial_step, temp, temp1, temp2
        integer  :: err, AllocateStatus, ind, ind2
        integer  :: i, j, k, l
        real(dl), allocatable :: param_array(:)
        real(dl), dimension(:), allocatable :: alpha_fiducial, alpha_derivative, alpha_temp, alpha_error
        real(dl), dimension(:,:), allocatable :: alpha_derivative_array


        ! allocate a parameter array:
        allocate( param_array(num_param) )
        ! write the parameters in the parameter array:
        call CAMB_params_to_params_array( P, FP, num_param, param_array )



        ! allocation:
        allocate(alpha_fiducial(FP%fisher_alpha%number_alpha_redshifts), stat = AllocateStatus)
        if (AllocateStatus /= 0) stop "Allocation failed. Not enough memory to allocate alpha_fiducial"
        allocate(alpha_temp(FP%fisher_alpha%number_alpha_redshifts), stat = AllocateStatus)
        if (AllocateStatus /= 0) stop "Allocation failed. Not enough memory to allocate alpha_temp"
        allocate(alpha_derivative(FP%fisher_alpha%number_alpha_redshifts), stat = AllocateStatus)
        if (AllocateStatus /= 0) stop "Allocation failed. Not enough memory to allocate alpha_derivative"
        allocate(alpha_error(FP%fisher_alpha%number_alpha_redshifts), stat = AllocateStatus)
        if (AllocateStatus /= 0) stop "Allocation failed. Not enough memory to allocate alpha_derivative"

        ! computation of the fiducial model:
        time_1 = omp_get_wtime()
        call alpha_Theo_Calc(P, FP, alpha_fiducial, err)
        alpha_error(:) = FP%fisher_alpha%alpha_error(:)
        time_2 = omp_get_wtime() - time_1

        if ( present(outroot) ) call save_alpha_mock_to_file(FP, FP%fisher_alpha%number_alpha_redshifts, alpha_fiducial, alpha_error, filename=TRIM(outroot)//'fiducial.dat')

        if ( err == 0 ) then
        else if ( err == 4 ) then
            write(*,*) 'Fiducial model unstable. Calculation cannot proceed.'
            stop
        else
            write(*,*) 'Error in computing the fiducial. Calculation cannot proceed.'
            stop
        end if

        ! print some feedback
        if ( FP%cosmicfish_feedback >= 1 ) then
            write(*,'(a)')         '**************************************************************'
            write(*,'(a,i10)')         'Some more feedback to be added '
            write(*,'(a,i10)')         'Number of parameters         : ', num_param
            write(*,'(a,f10.3,a)')     'Time taken for alpha variation calculation: ', time_2, ' (s)'
            write(*,'(a,i10)')         'Number of omp processes      : ', OMP_GET_MAX_THREADS()
            if ( FP%adaptivity ) then
                write(*,'(a,f10.3,a)') 'Estimated completion time    : ', time_2*20._dl*num_param, ' (s)'
            else
                write(*,'(a,f10.3,a)') 'Estimated completion time    : ', time_2*2._dl*num_param, ' (s)'
            end if
            write(*,'(a)')         '**************************************************************'
        end if

        ! compute the derivatives of the Cls:
        if ( FP%cosmicfish_feedback >= 1 ) then
            write(*,'(a)') '**************************************************************'
            write(*,'(a)') 'Computation of the delta alpha/alpha derivatives'
            write(*,'(a)') '**************************************************************'
        end if

        allocate(alpha_derivative_array(num_param,FP%fisher_alpha%number_alpha_redshifts), stat = AllocateStatus)
        if (AllocateStatus /= 0) stop "Allocation failed. Not enough memory to allocate rd_derivative_array"

        time_1 = omp_get_wtime()
        do ind = 1, num_param

            if ( FP%cosmicfish_feedback >= 1 ) then
                write(*,'(a,i5,a,E15.5)') 'Doing parameter number: ', ind, ' value: ', param_array(ind)
            end if

            if ( param_array(ind) /= 0._dl ) then
                initial_step = param_array(ind)*3.0_dl/100._dl !MM: shouldn't this be a user defined quantity
            else
                initial_step = 3.0_dl/100._dl
            end if

            call alpha_DerivativeCalc( P, FP, ind, alpha_derivative, alpha_temp, alpha_fiducial, initial_step, num_param, err)

            if (err /= 0) print*, 'WARNING: something went wrong with the delta v derivative of parameter:', ind

            if ( present(outroot) .and. FP%cosmicfish_feedback >= 2 ) then
                call save_alpha_mock_to_file(FP, FP%fisher_alpha%number_alpha_redshifts, alpha_derivative, alpha_temp, &
                filename=TRIM(outroot)//'varalpha_derivative_'//trim(adjustl(integer_to_string(ind)))//'.dat')
            end if

            alpha_derivative_array(ind,:) = alpha_derivative(:)

        end do
        time_2 = omp_get_wtime() - time_1

        if ( FP%cosmicfish_feedback >= 1 ) then
            write(*,*)
            write(*,'(a,f8.3,a)') 'Time to compute all the derivatives:', time_2, ' (s)'
        end if


        ! compute the Fisher matrix:
        do ind = 1, num_param
            do ind2 = 1, num_param

                temp2 = 0._dl

                do i = 1, FP%fisher_alpha%number_alpha_redshifts
                    temp = (1./alpha_error(i)**2.)*alpha_derivative_array(ind,i)*alpha_derivative_array(ind2,i)
                    temp2 = temp2 + temp
                end do

                Fisher_matrix(ind, ind2) = temp2

            end do
        end do
        time_2 = omp_get_wtime() - time_1

        if ( FP%cosmicfish_feedback >= 1 ) then
            write(*,'(a,f8.3,a)') 'Time to compute the Fisher matrix:', time_2, ' (s)'
        end if


    end subroutine Fisher_alpha

end module Fisher_calculator_alpha
