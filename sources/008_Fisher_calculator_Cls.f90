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

!> @file 008_Fisher_calculator_Cls.f90
!! This file contains the code that computes the Fisher matrix for Cls.

!----------------------------------------------------------------------------------------
!> This module contains the code that computes the Fisher matrix for Cls.

!> @author Marco Raveri

module Fisher_calculator_Cls

    use precision
    use ModelParams
    use cosmicfish_types
    use CAMBmain
    use CAMB
    use constants, only: const_pi
    use cosmicfish_camb_interface
    use lensnoise
    use Fisher_calculator_Cls_windows
    use omp_lib
    use cosmicfish_utilities
    use Matrix_utilities

#ifdef COSMICFISH_EFTCAMB
    use EFTinitialization
#endif

    implicit none

    private

    public dimension_cl_covariance, cl_covariance_out, save_cl_covariance_to_file, &
        & cls_derivative, add_noise_cls_to_fiducial, Fisher_Cls, Fiducial_Cls

contains

    ! ---------------------------------------------------------------------------------------------
    !> This subroutine computes the number of entrance of the cls covariance matrix
    !! The matrix of the Cls is made up by T, E, CMBWL, W1,..., Wn, B.
    subroutine dimension_cl_covariance( P, FP, num )

        implicit none

        Type(CAMBparams)      , intent(in)  :: P         !< Input CAMBparams
        Type(cosmicfish_params), intent(in)  :: FP        !< Input Cosmicfish params
        integer               , intent(out) :: num       !< Output size of the Cls covariance matrix

        integer  :: i

        num = 0

        if ( FP%fisher_cls%Fisher_want_CMB_T       ) num = num + 1
        if ( FP%fisher_cls%Fisher_want_CMB_E       ) num = num + 1
        if ( FP%fisher_cls%Fisher_want_CMB_lensing ) num = num + 1

        do i = 1, num_redshiftwindows
            if ( Redshift_w(i)%kind == window_counts .and. FP%fisher_cls%Fisher_want_LSS_counts ) then
                num = num + 1
            else if ( Redshift_w(i)%kind == window_lensing .and. FP%fisher_cls%Fisher_want_LSS_lensing ) then
                num = num + 1
            end if
        end do

        if ( FP%fisher_cls%Fisher_want_CMB_B       ) num = num + 1

    end subroutine dimension_cl_covariance

    ! ---------------------------------------------------------------------------------------------
    !> This subroutine that gives the covariance matrix of the Cls for the desired model.
    !! Calls internally CAMB_get results.
    subroutine cl_covariance_out( P, FP, cl_dim, l_min, l_max, cl_out, err )

        implicit none

        Type(CAMBparams)         :: P                                        !< Input CAMBparams
        Type(cosmicfish_params)  :: FP                                       !< Input Cosmicfish params
        integer, intent(in)      :: cl_dim                                   !< Input dimension of the Cl matrix
        integer, intent(in)      :: l_min                                    !< Input min l computed
        integer, intent(in)      :: l_max                                    !< Input max l computed

        real(dl)                 :: cl_out(cl_dim, cl_dim, l_max-l_min+1)    !< Output matrix containing the cls
        integer, intent(out)     :: err                                      !< Output error code:
                                                                             !< 0 = all fine
                                                                             !< 1 = error in input
                                                                             !< 2 = error in computation
                                                                             !< 3 = error in quality
                                                                             !< 4 = routine specific

        ! definitions:
        logical  :: EFTsuccess
        real(dl) :: ell
        real(dl), dimension(:,:), allocatable :: outarr
        integer , dimension(:)  , allocatable :: array_select
        integer , dimension(cl_dim) :: l_maxes, l_mins

        integer  :: i, j, k, l, m, temp, stat, ind
        real(dl) :: temp_factor

        ! initialization:
        cl_out = 0.0_dl
        err    = 0

        ! set the window function:
        if ( FP%fisher_cls%LSS_number_windows > 0 ) call init_camb_sources_windows( FP )

        ! set the camb parameters:
        call CAMBParams_Set(P)

#ifdef COSMICFISH_EFTCAMB
        ! check stability if EFTCAMB is used:
        if (P%EFTflag /= 0) then
            call EFTCAMB_initialization(EFTsuccess)
            if (.not. EFTsuccess) then
                write(*,*) 'EFTCAMB: unstable model'
                err = 4
                return
            end if
        end if
#endif

        ! call the CAMB_GetResults to compute the Cls.
        call CAMB_GetResults(P)

        ! test for errors in the computation
        if (global_error_flag/=0) then
            write (*,*) 'Error result '//trim(global_error_message)
            err = 2
            return
        endif

        ! use the scal Cls array as a starting point for the Cls covariance:

        allocate(outarr(1:3+num_redshiftwindows,1:3+num_redshiftwindows))
        allocate(array_select(3+num_redshiftwindows))

        ! decide which columns and rows of the scalar matrix to mantain:
        array_select = 0

        do i=1, 3+num_redshiftwindows

            if ( FP%fisher_cls%Fisher_want_CMB_T       ) array_select(1) = 1
            if ( FP%fisher_cls%Fisher_want_CMB_E       ) array_select(2) = 1
            if ( FP%fisher_cls%Fisher_want_CMB_lensing ) array_select(3) = 1

            do j = 1, num_redshiftwindows
                if ( Redshift_w(j)%kind == window_counts .and. FP%fisher_cls%Fisher_want_LSS_counts ) then
                    array_select(3+j) = 1
                else if ( Redshift_w(j)%kind == window_lensing .and. FP%fisher_cls%Fisher_want_LSS_lensing ) then
                    array_select(3+j) = 1
                end if
            end do

        end do

        ! cycle over the scalar Cls matrix.
        do i=l_min, l_max

            ell = REAL(i)
            temp_factor = 2._dl*pi/(ell*(ell+1._dl))

            outarr = temp_factor*Cl_scalar_array(i,1,1:3+num_redshiftwindows,1:3+num_redshiftwindows)
            outarr(1:2,:)=sqrt(FP%output_factor)*outarr(1:2,:)
            outarr(:,1:2)=sqrt(FP%output_factor)*outarr(:,1:2)

            l = 1
            m = 1

            do j=1, 3+num_redshiftwindows
                if ( array_select(j)==1 ) then
                    m = 1
                    do k=1, 3+num_redshiftwindows
                        if ( array_select(k)==1 ) then
                            cl_out(l, m, i-l_min+1) = outarr( j, k )
                            m = m + 1
                        end if
                    end do
                    l = l + 1
                end if
            end do

        end do

        ! now replace specific ones with lensed or add B mode power spectrum
        if ( P%WantTensors .and. .not. P%DoLensing) then
            ! if want tensor not lensing out = total
            do i = l_min, min( CP%Max_l_tensor, l_max )
                ell = REAL(i)
                temp_factor = 2*const_pi/(ell*(ell+1))
                j           = i-l_min+1

                ! TT
                if ( FP%fisher_cls%Fisher_want_CMB_T ) cl_out(1,1,j) = temp_factor*FP%output_factor*(Cl_scalar(i, 1, C_Temp)+ Cl_tensor(i, 1, C_Temp))
                ! EE and TE
                if ( FP%fisher_cls%Fisher_want_CMB_T .and. FP%fisher_cls%Fisher_want_CMB_E ) then
                    cl_out(2,2,j) = temp_factor*FP%output_factor*(Cl_scalar(i, 1, C_E)+ Cl_tensor(i, 1, C_E))
                    cl_out(1,2,j) = temp_factor*FP%output_factor*(Cl_scalar(i, 1, C_Cross) + Cl_tensor(i, 1, CT_Cross))
                end if
                ! if only EE
                if ( .not. FP%fisher_cls%Fisher_want_CMB_T .and. FP%fisher_cls%Fisher_want_CMB_E ) then
                    cl_out(1,1,j) = temp_factor*FP%output_factor*(Cl_scalar(i, 1, C_E)+ Cl_tensor(i, 1, C_E))
                end if

                ! BB
                if ( FP%fisher_cls%Fisher_want_CMB_B ) cl_out( cl_dim, cl_dim, j ) = temp_factor*FP%output_factor*Cl_tensor(i, 1, CT_B)

            end do
            do i = CP%Max_l_tensor+1,l_max
                ell = REAL(i)
                temp_factor = 2*const_pi/(ell*(ell+1))
                j           = i-l_min+1

                ! TT
                if ( FP%fisher_cls%Fisher_want_CMB_T ) cl_out(1,1,j) = temp_factor*FP%output_factor*Cl_scalar(i, 1, C_Temp)
                ! EE and TE
                if ( FP%fisher_cls%Fisher_want_CMB_T .and. FP%fisher_cls%Fisher_want_CMB_E ) then
                    cl_out(2,2,j) = temp_factor*FP%output_factor*Cl_scalar(i, 1, C_E)
                    cl_out(1,2,j) = temp_factor*FP%output_factor*Cl_scalar(i, 1, C_Cross)
                end if
                ! if only EE
                if ( .not. FP%fisher_cls%Fisher_want_CMB_T .and. FP%fisher_cls%Fisher_want_CMB_E ) then
                    cl_out(1,1,j) = temp_factor*FP%output_factor*Cl_scalar(i, 1, C_E)
                end if

                ! BB
                if ( FP%fisher_cls%Fisher_want_CMB_B ) cl_out( cl_dim, cl_dim, j ) = 0._dl

            end do
        else if (.not. P%WantTensors .and. P%DoLensing) then
            ! if want lensing and not tensor out = lensed
            do i = l_min, min( lmax_lensed, l_max )
                ell = REAL(i)
                temp_factor = 2*const_pi/(ell*(ell+1))
                j           = i-l_min+1

                ! TT
                if ( FP%fisher_cls%Fisher_want_CMB_T ) cl_out(1,1,j) = temp_factor*FP%output_factor*Cl_lensed(i, 1, CT_Temp)
                ! EE and TE
                if ( FP%fisher_cls%Fisher_want_CMB_T .and. FP%fisher_cls%Fisher_want_CMB_E ) then
                    cl_out(2,2,j) = temp_factor*FP%output_factor*Cl_lensed(i, 1, CT_E)
                    cl_out(1,2,j) = temp_factor*FP%output_factor*Cl_lensed(i, 1, CT_Cross)
                end if
                ! if only EE
                if ( .not. FP%fisher_cls%Fisher_want_CMB_T .and. FP%fisher_cls%Fisher_want_CMB_E ) then
                    cl_out(1,1,j) = temp_factor*FP%output_factor*Cl_lensed(i, 1, CT_E)
                end if

                ! BB
                if ( FP%fisher_cls%Fisher_want_CMB_B ) cl_out( cl_dim, cl_dim, j ) = temp_factor*FP%output_factor*Cl_lensed(i, 1, CT_B)
            end do
        else if (P%WantTensors .and. P%DoLensing) then
            ! if want lensing and want tensor out = total lensed
            do i = l_min, min( CP%Max_l_tensor, lmax_lensed, l_max )
                ell = REAL(i)
                temp_factor = 2*const_pi/(ell*(ell+1))
                j           = i-l_min+1

                ! TT
                if ( FP%fisher_cls%Fisher_want_CMB_T ) cl_out(1,1,j) = temp_factor*FP%output_factor*(Cl_lensed(i, 1, CT_Temp) +Cl_tensor(i, 1, CT_Temp))
                ! EE and TE
                if ( FP%fisher_cls%Fisher_want_CMB_T .and. FP%fisher_cls%Fisher_want_CMB_E ) then
                    cl_out(2,2,j) = temp_factor*FP%output_factor*(Cl_lensed(i, 1, CT_E) +Cl_tensor(i, 1, CT_E))
                    cl_out(1,2,j) = temp_factor*FP%output_factor*(Cl_lensed(i, 1, CT_Cross) +Cl_tensor(i, 1, CT_Cross))
                end if
                ! if only EE
                if ( .not. FP%fisher_cls%Fisher_want_CMB_T .and. FP%fisher_cls%Fisher_want_CMB_E ) then
                    cl_out(1,1,j) = temp_factor*FP%output_factor*(Cl_lensed(i, 1, CT_E) +Cl_tensor(i, 1, CT_E))
                end if

                ! BB
                if ( FP%fisher_cls%Fisher_want_CMB_B ) cl_out( cl_dim, cl_dim, j ) = temp_factor*FP%output_factor*(Cl_lensed(i, 1, CT_B) +Cl_tensor(i, 1, CT_B))
            end do
            do i = min( CP%Max_l_tensor, lmax_lensed, l_max )+1, min( lmax_lensed, l_max )
                ell = REAL(i)
                temp_factor = 2*const_pi/(ell*(ell+1))
                j           = i-l_min+1

                ! TT
                if ( FP%fisher_cls%Fisher_want_CMB_T ) cl_out(1,1,j) = temp_factor*FP%output_factor*Cl_lensed(i, 1, CT_Temp)
                ! EE and TE
                if ( FP%fisher_cls%Fisher_want_CMB_T .and. FP%fisher_cls%Fisher_want_CMB_E ) then
                    cl_out(2,2,j) = temp_factor*FP%output_factor*Cl_lensed(i, 1, CT_E)
                    cl_out(1,2,j) = temp_factor*FP%output_factor*Cl_lensed(i, 1, CT_Cross)
                end if
                ! if only EE
                if ( .not. FP%fisher_cls%Fisher_want_CMB_T .and. FP%fisher_cls%Fisher_want_CMB_E ) then
                    cl_out(1,1,j) = temp_factor*FP%output_factor*Cl_lensed(i, 1, CT_E)
                end if

                ! BB
                if ( FP%fisher_cls%Fisher_want_CMB_B ) cl_out( cl_dim, cl_dim, j ) = temp_factor*FP%output_factor*Cl_lensed(i, 1, CT_B)
            end do
        end if

        ! prepare the l_max cut:

        l_maxes = 1
        ind     = 1

        if ( FP%fisher_cls%Fisher_want_CMB_T ) then
            l_maxes(ind) = FP%fisher_cls%l_max_TT
            ind = ind + 1
        end if

        if ( FP%fisher_cls%Fisher_want_CMB_E ) then
            l_maxes(ind) = FP%fisher_cls%l_max_EE
            ind = ind + 1
        end if

        if ( FP%fisher_cls%Fisher_want_CMB_lensing ) then
            l_maxes(ind) = min( FP%fisher_cls%l_max_TT, &
                & FP%fisher_cls%l_max_EE, &
                & FP%fisher_cls%l_max_BB  )
            ind = ind + 1
        end if

        if ( FP%fisher_cls%Fisher_want_CMB_B ) then
            l_maxes(cl_dim) = FP%fisher_cls%l_max_BB
        end if

        do k = 1, FP%fisher_cls%LSS_number_windows
            l_maxes(ind) = FP%fisher_cls%LSS_lmax(k)
            ind = ind + 1
        end do

        ! prepare the l_min cut:

         l_mins  = 1
         ind     = 1

         if ( FP%fisher_cls%Fisher_want_CMB_T ) then
             l_mins(ind) = FP%fisher_cls%l_min_TT
             ind = ind + 1
         end if

         if ( FP%fisher_cls%Fisher_want_CMB_E ) then
             l_mins(ind) = FP%fisher_cls%l_min_EE
             ind = ind + 1
         end if

         if ( FP%fisher_cls%Fisher_want_CMB_lensing ) then
             l_mins(ind) = max( FP%fisher_cls%l_min_TT, &
                 & FP%fisher_cls%l_min_EE, &
                 & FP%fisher_cls%l_min_BB  )
             ind = ind + 1
         end if

         if ( FP%fisher_cls%Fisher_want_CMB_B ) then
             l_mins(cl_dim) = FP%fisher_cls%l_min_BB
         end if

         do k = 1, FP%fisher_cls%LSS_number_windows
             l_mins(ind) = FP%fisher_cls%LSS_lmin(k)
             ind = ind + 1
         end do

        ! apply the l cut:
        do i=l_min, l_max
            do j=1, cl_dim
                do k=1, j

                    if ( i <= min(l_maxes(j), l_maxes(k)) .and. i >= max(l_mins(j), l_mins(k)) ) then
                        cl_out( k, j, i-l_min+1) = cl_out( k, j, i-l_min+1)
                    else
                        cl_out( k, j, i-l_min+1) = 0._dl
                    end if

                end do
            end do
        end do

        ! ensure the fact that the covariance is symmetric:
        do i=l_min, l_max
            do j=1, cl_dim
                do k=1, j

                    cl_out( j, k, i-l_min+1) = cl_out( k, j, i-l_min+1)

                end do
            end do
        end do

        ! remove cross correlation if it is not wanted:
        if ( .not. FP%fisher_cls%Fisher_want_XC ) then
            do i=l_min, l_max
                do j=1, cl_dim
                    do k=1, cl_dim

                        if ( j/=k ) then
                            cl_out( j, k, i-l_min+1) = 0._dl
                        end if

                    end do
                end do
            end do
        end if

    end subroutine cl_covariance_out

    ! ---------------------------------------------------------------------------------------------
    !> This subroutine saves to file the input Cls covariance matrix
    subroutine save_cl_covariance_to_file( cl_cov_in, cl_dim, l_min, l_max, filename )

        implicit none

        integer         , intent(in)  :: cl_dim         !< Number of Cls in a row of the matrix
        integer         , intent(in)  :: l_min          !< Minimum l of the Cls
        integer         , intent(in)  :: l_max          !< Maximum l of the Cls
        real(dl), dimension(cl_dim, cl_dim, l_max-l_min+1), intent(in) :: cl_cov_in  !< Input Cls covariance matrix
        character(len=*), intent(in)  :: filename       !< Name of the file where the Cls covariance will be saved

        integer  :: i, j, k
        real(dl) :: factor
        real(dl), dimension(cl_dim, cl_dim) :: dummy

        open(unit=666, FILE=filename, ACTION="write", STATUS="replace")
        ! cycle over l
        do i = 1, l_max-l_min+1
            ! change this factor if want to compensate for factors l()...
            factor     = 1._dl
            dummy(:,:) = cl_cov_in(:,:,i)
            write(666,'(E15.5)',advance='no') REAL(i+l_min-1)
            ! cycle over the Cl cov matrix at a fixed l
            do j = 1, cl_dim
                do k = 1, cl_dim
                    write(666,'(ES20.5E3)',advance='no') factor*dummy(j,k)
                end do
            end do
            ! go to newline
            write(666,'(a)') ' '
        end do
        close(666)

    end subroutine

    ! ---------------------------------------------------------------------------------------------
    !> This subroutine computes the derivative of the Cls covariance matrix with a fixed stepsize
    !> finite difference formula or by the Richardson’s deferred approach to the limit.
    subroutine cls_derivative( P, FP, param_num, &
        & cl_initial, cl_derivative, cl_error, &
        & initial_step, &
        & out_dim, l_min, l_max, num_param, &
        & error_code )

        implicit none

        Type(CAMBparams)       , intent(in) :: P       !< Input CAMBparams
        Type(cosmicfish_params) , intent(in) :: FP      !< Input Cosmicfish params

        integer , intent(in)  :: param_num  !< Input number of the parameter wrt the derivative is computed
        integer , intent(in)  :: out_dim    !< Input dimension of the Cls matrix
        integer , intent(in)  :: l_min      !< Input minimum value of l
        integer , intent(in)  :: l_max      !< Input maximum value of l
        integer , intent(in)  :: num_param  !< Input number of parameters involved in the calculation

        real(dl), dimension(out_dim, out_dim, l_max-l_min+1), intent(in)  :: cl_initial     !< Cls covariance matrix computed for h=0
        real(dl), dimension(out_dim, out_dim, l_max-l_min+1), intent(out) :: cl_derivative  !< Output returns the derivative of the Cls
        real(dl), dimension(out_dim, out_dim, l_max-l_min+1), intent(out) :: cl_error       !< Output returns the error on the Cls derivative at a given l

        real(dl), intent(in)  :: initial_step  !< Input initial value of h. Does not need to be small for the adaptive algorithm.
        integer , intent(out) :: error_code    !< Output error code. When run the code will not stop on error but report it
                                               !<  0 = everything was fine computation exited correctly
                                               !<  1 = error on input parameters
                                               !<  2 = error on computation
                                               !<  3 = error on quality of the results

        ! definitions:
        real(dl), dimension(:,:,:), allocatable :: cl_temp_1, cl_temp_2, errt
        integer , dimension(:,:,:), allocatable :: accept_mask
        real(dl), dimension(:,:,:,:,:), allocatable :: a
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

        ! initialize output variables:
        error_code    = 0
        cl_derivative = 0._dl
        cl_error      = 0._dl

        ! check input values are correct
        if ( initial_step == 0._dl ) then
            stop 'cls_derivative initial stepsize is zero'
        end if

        ! allocate the arrays:

        if ( FP%adaptivity ) then
            allocate( cl_temp_1(out_dim, out_dim, l_max-l_min+1),  &
                & cl_temp_2(out_dim, out_dim, l_max-l_min+1),  &
                & errt(out_dim, out_dim, l_max-l_min+1),       &
                & accept_mask(out_dim, out_dim, l_max-l_min+1),&
                & stat = AllocateStatus )
            if (AllocateStatus /= 0) stop "cls_derivative: Allocation failed. Not enough memory"
            allocate( a(ntab, ntab, out_dim, out_dim, l_max-l_min+1) , stat = AllocateStatus )
            if (AllocateStatus /= 0) stop "cls_derivative: Allocation failed. Not enough memory for the derivative array"
        else
            allocate( cl_temp_1(out_dim, out_dim, l_max-l_min+1),  &
                & cl_temp_2(out_dim, out_dim, l_max-l_min+1),  &
                & stat = AllocateStatus )
            if (AllocateStatus /= 0) stop "cls_derivative: Allocation failed. Not enough memory"
            allocate( a(1,1, out_dim, out_dim, l_max-l_min+1) , stat = AllocateStatus )
            if (AllocateStatus /= 0) stop "cls_derivative: Allocation failed. Not enough memory for the derivative array"
        end if

        ! initialize the parameter array:
        P_temp  = P
        FP_temp = FP

        call CAMB_params_to_params_array( P_temp, FP, num_param, param_vector )

        fiducial_param_vector = param_vector
        temp_param_vector     = param_vector

        ! initialize the computation:
        hh          = initial_step
        cl_error    = big
        side        = 0 ! both sides
        if ( FP%adaptivity ) accept_mask = 0

        ! do the initial step on the right:
        param_vector(param_num)      = fiducial_param_vector(param_num) + hh
        call params_array_to_CAMB_param( P_temp, FP_temp, num_param, param_vector )
        call cl_covariance_out( P_temp, FP_temp, out_dim, l_min, l_max, cl_temp_1, err )
        ! if on the right we cannot go we go to the left:
        if ( err /= 0 ) side = 1
        ! do the initial step on the left:
        temp_param_vector(param_num) = fiducial_param_vector(param_num) - hh
        call params_array_to_CAMB_param( P_temp, FP_temp, num_param, temp_param_vector )
        call cl_covariance_out( P_temp, FP_temp, out_dim, l_min, l_max, cl_temp_2, err )

        ! if on the left we cannot go decide what to do:
        if ( err /= 0 ) then
            if ( side == 1 ) then
                ! we cannot go to the left nor to the right. Return derivative zero with error code.
                error_code     = 4
                cl_error       = 0.0_dl
                cl_derivative  = 0.0_dl
                return
            else if ( side == 0 ) then
                ! we cannot go to the left but we can go to the right:
                side = 2
            end if
        end if

        ! store the results based on the side we decided to go:
        if (side == 0) then ! both sides
            a(1,1,:,:,:) = (cl_temp_1(:,:,:) - cl_temp_2(:,:,:))/(2.0_dl*hh)
        else if (side == 1) then ! left side: increment negative
            a(1,1,:,:,:) = (cl_temp_2(:,:,:) - cl_initial(:,:,:))/(-hh)
        else if (side == 2) then ! right side: increment positive
            a(1,1,:,:,:) = (cl_temp_1(:,:,:) - cl_initial(:,:,:))/(hh)
        end if

        ! fixed stepsize derivative if adaptivity is not wanted:
        if ( .not. FP%adaptivity ) then
            cl_derivative(:,:,:) = a(1,1,:,:,:)
            cl_error             = 0.0_dl
            deallocate( cl_temp_1, cl_temp_2, a )
            return
        end if

        ! start the derivative computation:
        do i = 2, ntab

            ! feedback for the adaptive derivative:
            if ( FP%cosmicfish_feedback >= 2 ) then
                accepted   = sum(accept_mask)
                acceptance = REAL(accepted)/REAL(out_dim*out_dim*(l_max-l_min+1))
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
                call cl_covariance_out( P_temp, FP_temp, out_dim, l_min, l_max, cl_temp_1, err )
                if ( err /= 0 ) then
                    error_code     = 4
                    cl_error       = 0.0_dl
                    cl_derivative  = 0.0_dl
                    return
                end if
                ! compute for - hh
                call params_array_to_CAMB_param( P_temp, FP_temp, num_param, temp_param_vector )
                call cl_covariance_out( P_temp, FP_temp, out_dim, l_min, l_max, cl_temp_2, err )
                if ( err /= 0 ) then
                    error_code     = 4
                    cl_error       = 0.0_dl
                    cl_derivative  = 0.0_dl
                    return
                end if
                ! store the result
                a(1,i,:,:,:) = (cl_temp_1(:,:,:) - cl_temp_2(:,:,:))/(2.0_dl*hh)
            else if (side==1) then
                param_vector(param_num)      = fiducial_param_vector(param_num) - hh
                call params_array_to_CAMB_param( P_temp, FP_temp, num_param, param_vector )
                call cl_covariance_out( P_temp, FP_temp, out_dim, l_min, l_max, cl_temp_1, err )
                if ( err /= 0 ) then
                    error_code     = 4
                    cl_error       = 0.0_dl
                    cl_derivative  = 0.0_dl
                    return
                end if
                ! store the result
                a(1,i,:,:,:) = (cl_temp_1(:,:,:) - cl_initial(:,:,:))/(-hh)
            else if (side==2) then
                param_vector(param_num)      = fiducial_param_vector(param_num) + hh
                call params_array_to_CAMB_param( P_temp, FP_temp, num_param, param_vector )
                call cl_covariance_out( P_temp, FP_temp, out_dim, l_min, l_max, cl_temp_1, err )
                if ( err /= 0 ) then
                    error_code     = 4
                    cl_error       = 0.0_dl
                    cl_derivative  = 0.0_dl
                    return
                end if
                ! store the result
                a(1,i,:,:,:) = (cl_temp_1(:,:,:) - cl_initial(:,:,:))/(hh)
            end if

            ! now do the extrapolation table:
            fac=CON2
            do j = 2, i

                forall ( k = 1:out_dim, l = 1:out_dim, m = 1:l_max-l_min+1, &
                    & accept_mask(k,l,m) == 0 )
                    ! compute the interpolation table:
                    a(j,i,k,l,m) = (a(j-1,i,k,l,m)*fac-a(j-1,i-1,k,l,m))/(fac-1._dl)
                    ! estimate the error:
                    errt(k,l,m)  = max( abs(a(j,i,k,l,m)-a(j-1,i,k,l,m)), &
                        & abs(a(j,i,k,l,m)-a(j-1,i-1,k,l,m)))
                end forall

                fac=CON2*fac

                ! if there is an improvement save the result:

                forall ( k = 1:out_dim, l = 1:out_dim, m = 1:l_max-l_min+1, &
                    & errt(k,l,m) .le. cl_error(k,l,m) .and. accept_mask(k,l,m) == 0 )

                    cl_error(k,l,m) = errt(k,l,m)
                    cl_derivative(k,l,m) = a(j,i,k,l,m)
                end forall

            end do

            ! accept a value of the derivative if the value of the extrapolation is stable
            ! Store this information into the mask matrix
            ! so that the value contained in cl_error and cl_derivative is not overwritten.
            forall ( k = 1:out_dim, l = 1:out_dim, m = 1:l_max-l_min+1, &
                & accept_mask(k,l,m) == 0 .and. &
                & (a(i,i,k,l,m) - a(i-1,i-1,k,l,m)) .ge. safe*cl_error(k,l,m) )

                accept_mask(k,l,m) = 1
            end forall

            ! if all the entries match the required accuracy exit the loop
            if ( sum(accept_mask) == out_dim*out_dim*(l_max-l_min+1) ) exit

        end do

        ! quality check on the results:
        if ( sum(accept_mask) /= out_dim*out_dim*(l_max-l_min+1) ) then
            error_code = 3
        end if

        ! deallocate the arrays:
        if ( FP%adaptivity ) then
            deallocate( cl_temp_1, cl_temp_2, errt, accept_mask, a )
        end if

        return

    end subroutine cls_derivative

    ! ---------------------------------------------------------------------------------------------
    !> This subroutine computes and adds the noise cls to the fiducial cls.
    !> On output the fiducial cls will be modified.
    subroutine add_noise_cls_to_fiducial( P, FP, cl_fiducial, out_dim, l_min, l_max, error_code)

        implicit none

        Type(CAMBparams)       , intent(in) :: P       !< Input CAMBparams
        Type(cosmicfish_params) , intent(in) :: FP      !< Input Cosmicfish params

        integer , intent(in)  :: out_dim    !< Input dimension of the Cls matrix
        integer , intent(in)  :: l_min      !< Input minimum value of l
        integer , intent(in)  :: l_max      !< Input maximum value of l

        real(dl), dimension(out_dim, out_dim, l_max-l_min+1), intent(inout)  :: cl_fiducial     !< Fiducial cls covariance. On output this contains the fiducial plus the noise cls.

        integer , intent(out) :: error_code    !< Output error code. When run the code will not stop on error but report it
                                               !<  0 = everything was fine computation exited correctly
                                               !<  1 = error on input parameters
                                               !<  2 = error on computation
                                               !<  3 = error on quality of the results

        ! definitions:

        integer  :: jjmaxTmax = 200  !< Maximum value for FUTURCMB. Should work for reasonable l_max.

        integer  :: ind, i, j, k
        real(dl) :: temp1, temp2, ell, factor
        real(dl) :: Compute_CMB_TT_Noise, Compute_CMB_EE_Noise, Compute_CMB_BB_Noise
        real(dl), allocatable :: NoiseVar(:), NoiseVarPol(:), sigma2(:)
        real(dl), allocatable :: incls(:,:), inlcls(:,:), clnoise(:,:), nldd(:)
        real(dl), allocatable :: fskyes(:)
        logical  :: EFTsuccess

        ! initialization:
        error_code = 0
        ind        = 1

        ! initialize CMB stuff:
        if ( FP%fisher_cls%Fisher_want_CMB_T .or. FP%fisher_cls%Fisher_want_CMB_E .or. &
            & FP%fisher_cls%Fisher_want_CMB_B .or. FP%fisher_cls%Fisher_want_CMB_lensing ) then

            allocate( NoiseVar(FP%fisher_cls%CMB_n_channels), &
                & NoiseVarPol(FP%fisher_cls%CMB_n_channels), &
                & sigma2(FP%fisher_cls%CMB_n_channels) )

            temp1 = P%Tcmb*const_pi/180._dl/60._dl
            temp2 = 180._dl*60._dl*sqrt(8._dl*log(2._dl))/const_pi

            forall( i = 1:FP%fisher_cls%CMB_n_channels )
                NoiseVar(i)    = ( FP%fisher_cls%CMB_temp_sens(i)*FP%fisher_cls%CMB_fwhm(i)*temp1 )**2
                NoiseVarPol(i) = ( FP%fisher_cls%CMB_pol_sens(i)*FP%fisher_cls%CMB_fwhm(i)*temp1  )**2
                sigma2(i)      = -1._dl*( FP%fisher_cls%CMB_fwhm(i)/temp2 )**2
            end forall

        end if

        ! add TT noise:
        if ( FP%fisher_cls%Fisher_want_CMB_T ) then

            do j = l_min, l_max
                ! compute noise:
                ell = REAL(j)
                Compute_CMB_TT_Noise = 0._dl
                do i=1, FP%fisher_cls%CMB_n_channels
                    Compute_CMB_TT_Noise = Compute_CMB_TT_Noise &
                        & + 1._dl/NoiseVar(i)*exp(ell*(ell+1)*sigma2(i) )
                end do
                Compute_CMB_TT_Noise = 1._dl/Compute_CMB_TT_Noise
                ! add noise to fiducial:
                k = j - l_min + 1
                cl_fiducial( ind, ind, k ) = cl_fiducial( ind, ind, k ) + Compute_CMB_TT_Noise
            end do
            ind = ind + 1

        end if
        ! add EE noise:
        if ( FP%fisher_cls%Fisher_want_CMB_E ) then

            do j = l_min, l_max
                ! compute noise:
                ell = REAL(j)
                Compute_CMB_EE_Noise = 0._dl
                do i=1, FP%fisher_cls%CMB_n_channels
                    Compute_CMB_EE_Noise = Compute_CMB_EE_Noise &
                        & + 1._dl/NoiseVarPol(i)*exp(ell*(ell+1)*sigma2(i) )
                end do
                Compute_CMB_EE_Noise = 1._dl/Compute_CMB_EE_Noise
                ! add noise to fiducial:
                k = j - l_min + 1
                cl_fiducial( ind, ind, k ) = cl_fiducial( ind, ind, k ) + Compute_CMB_EE_Noise
            end do
            ind = ind + 1

        end if
        ! add lensing noise:
        if ( FP%fisher_cls%Fisher_want_CMB_lensing ) then

            ! we need scalar and lensed scalar Cls. We have to call camb and we want the model to be stable
#ifdef COSMICFISH_EFTCAMB
            ! check stability if EFTCAMB is used:
            if (P%EFTflag /= 0) then
                call CAMBParams_Set(P)
                call EFTCAMB_initialization(EFTsuccess)
                if (.not. EFTsuccess) then
                    write(*,*) 'EFTCAMB: unstable fiducial model. Calculation cannot proceed.'
                    stop
                end if
            end if
#endif
            call CAMB_GetResults(P)

            ! test for errors in the computation:
            if (global_error_flag/=0) then
                write (*,*) 'Error result '//trim(global_error_message)
                write (*,*) 'add_noise_cls_to_fiducial: calculation of fiducial model failed. Calculation cannot proceed.'
                stop
            endif

            allocate(incls  (5, 1:l_max))
            allocate(inlcls (4, 1:l_max))
            allocate(clnoise(2, 1:l_max))
            allocate(nldd(1:l_max))
            ! initialization of the Cls matrices:
            incls     = 0._dl
            inlcls    = 0._dl
            clnoise   = 0._dl
            nldd      = 0._dl

            ! write the Cls in the proper places:
            do i = l_min, l_max

                ell    = REAL(i)
                factor = 2*const_pi/(ell*(ell+1))
                ! write the cls in the format that FUTURECMB likes:
                incls(1,i) = FP%output_factor*factor*Cl_scalar_array(i,1,1,1)         ! TT
                incls(2,i) = FP%output_factor*factor*Cl_scalar_array(i,1,2,2)         ! EE
                incls(3,i) = FP%output_factor*factor*Cl_scalar_array(i,1,1,2)         ! TE
                incls(4,i) = 2*const_pi*Cl_scalar_array(i,1,3,3)                      ! dd
                incls(5,i) = sqrt(FP%output_factor)*sqrt(ell*(ell+1))*factor*Cl_scalar_array(i,1,1,3)        ! Td

                inlcls(1,i) = FP%output_factor*factor*Cl_lensed(i, 1, CT_Temp)        ! TT ! micro K^2
                inlcls(2,i) = FP%output_factor*factor*Cl_lensed(i, 1, CT_E)           ! EE ! micro K^2
                inlcls(3,i) = FP%output_factor*factor*Cl_lensed(i, 1, CT_Cross)       ! TE ! micro K^2
                inlcls(4,i) = FP%output_factor*factor*Cl_lensed(i, 1, CT_B)           ! BB ! micro K^2

                ell = REAL(i)
                Compute_CMB_TT_Noise = 0._dl
                do j=1, FP%fisher_cls%CMB_n_channels
                    Compute_CMB_TT_Noise = Compute_CMB_TT_Noise &
                        & + 1._dl/NoiseVar(j)*exp(ell*(ell+1)*sigma2(j) )
                end do
                Compute_CMB_TT_Noise = 1._dl/Compute_CMB_TT_Noise
                clnoise(1,i) = Compute_CMB_TT_Noise

                Compute_CMB_EE_Noise = 0._dl
                do j=1, FP%fisher_cls%CMB_n_channels
                    Compute_CMB_EE_Noise = Compute_CMB_EE_Noise &
                        & + 1._dl/NoiseVarPol(j)*exp(ell*(ell+1)*sigma2(j) )
                end do
                Compute_CMB_EE_Noise = 1._dl/Compute_CMB_EE_Noise
                clnoise(2,i) = Compute_CMB_EE_Noise

            end do

            ! compute the lensing noise with FUTURCMB:
            call calc_Nldd(incls, inlcls, clnoise, nldd, jjmaxTmax, l_max)

            ! add noise to the fiducial:
            do i = l_min, l_max
                ell    = REAL(i)
                factor = 1._dl/(ell*(ell+1))
                k = i - l_min + 1
                cl_fiducial( ind, ind, k ) = cl_fiducial( ind, ind, k ) + factor*nldd(i)
            end do

            ind = ind + 1

        end if
        ! add BB noise:
        if ( FP%fisher_cls%Fisher_want_CMB_b ) then

            do j = l_min, l_max
                ! compute noise:
                ell = REAL(j)
                Compute_CMB_BB_Noise = 0._dl
                do i=1, FP%fisher_cls%CMB_n_channels
                    Compute_CMB_BB_Noise = Compute_CMB_BB_Noise &
                        & + 1._dl/NoiseVarPol(i)*exp(ell*(ell+1)*sigma2(i) )
                end do
                Compute_CMB_BB_Noise = 1._dl/Compute_CMB_BB_Noise
                ! add noise to fiducial:
                k = j - l_min + 1
                cl_fiducial( out_dim, out_dim, k ) = cl_fiducial( out_dim, out_dim, k ) + Compute_CMB_BB_Noise
            end do

        end if

        ! add counts and lensing windows noises:
        j = 1
        do k = 1, num_redshiftwindows
            ! counts noise
            if ( Redshift_w(k)%kind == window_counts .and. FP%fisher_cls%Fisher_want_LSS_counts ) then
                do i = 1, l_max-l_min+1
                    cl_fiducial( ind, ind, i ) = &
                        & cl_fiducial( ind, ind, i ) + 1._dl/FP%fisher_cls%LSS_num_galaxies(j)
                end do
                ind = ind + 1
                j   = j + 1
            ! lensing noise
            else if ( Redshift_w(k)%kind == window_lensing .and. FP%fisher_cls%Fisher_want_LSS_lensing ) then
                do i = 1, l_max-l_min+1
                    cl_fiducial( ind, ind, i ) = &
                        & cl_fiducial( ind, ind, i ) &
                        & +FP%fisher_cls%LSS_intrinsic_ellipticity(j)**2/FP%fisher_cls%LSS_num_galaxies(j)
                end do
                ind = ind + 1
                j   = j + 1
            end if
        end do

        ! multiply by f_sky:

        ! prepare an f_sky vector:
        allocate( fskyes(out_dim) )
        fskyes = 1._dl

        ind = 1

        if ( FP%fisher_cls%Fisher_want_CMB_T ) then
            fskyes(ind) = FP%fisher_cls%CMB_TT_fsky
            ind = ind + 1
        end if

        if ( FP%fisher_cls%Fisher_want_CMB_E ) then
            fskyes(ind) = FP%fisher_cls%CMB_EE_fsky
            ind = ind + 1
        end if

        if ( FP%fisher_cls%Fisher_want_CMB_lensing ) then
            fskyes(ind) = min( FP%fisher_cls%CMB_TT_fsky, &
                & FP%fisher_cls%CMB_EE_fsky, &
                & FP%fisher_cls%CMB_BB_fsky )
            ind = ind + 1
        end if

        if ( FP%fisher_cls%Fisher_want_CMB_B ) then
            fskyes(out_dim) = FP%fisher_cls%CMB_BB_fsky
        end if

        do k = 1, FP%fisher_cls%LSS_number_windows
            fskyes(ind) = FP%fisher_cls%LSS_fsky(k)
            ind = ind + 1
        end do

        ! multiply in the appropriate manner the fsky vector to the Cls matrix:
        do j=1, out_dim
            do k=1, out_dim

                cl_fiducial(j, k, :) = cl_fiducial(j, k, :)/sqrt(sqrt(fskyes(j)*fskyes(k)))

            end do
        end do

    end subroutine add_noise_cls_to_fiducial

    ! ---------------------------------------------------------------------------------------------
    !> This subroutine computes the Fisher matrix for Cls.
    subroutine Fisher_Cls( P, FP, num_param, Fisher_Matrix, outroot )

        implicit none

        Type(CAMBparams)       , intent(in) :: P       !< Input CAMBparams
        Type(cosmicfish_params) , intent(in) :: FP      !< Input Cosmicfish params

        integer , intent(in)  :: num_param                                        !< Input dimension of the Fisher matrix
        real(dl), dimension(num_param,num_param), intent(out) :: Fisher_Matrix    !< Output Fisher matrix
        character(len=*), intent(in), optional :: outroot                         !< Optional input: filename that will be used to dump to file optional output if feedback is greater than 1.

        real(dl) :: time_1, time_2
        real(dl) :: initial_step, temp, temp1, temp2
        integer  :: cl_dim, l_min, l_max, err, AllocateStatus, ind, ind2
        integer  :: i, j, k, l
        real(dl), allocatable :: param_array(:)
        real(dl), dimension(:,:,:), allocatable :: cl_fiducial, cl_derivative, cl_temp
        real(dl), dimension(:,:,:,:), allocatable :: cl_derivative_array
        real(dl), dimension(:,:)    , allocatable :: dummy_matrix_1, dummy_matrix_2
        real(dl), dimension(:,:)    , allocatable :: dummy_matrix_3, dummy_matrix_4, dummy_matrix_5

        ! allocate a parameter array:
        allocate( param_array(num_param) )
        ! write the parameters in the parameter array:
        call CAMB_params_to_params_array( P, FP, num_param, param_array )
        ! get the size of the cls covariance matrix:
        call dimension_cl_covariance( P, FP, cl_dim )
        ! protect against errors:
        if ( cl_dim == 0 ) then
            write(*,*) 'WARNING: the Cls matrix is of zero dimension'
            Fisher_matrix = 0._dl
            return
        end if

        ! decide l_min and l_max:
        l_min = 2
        l_max = maxval( FP%fisher_cls%LSS_lmax, FP%fisher_cls%LSS_number_windows )
        l_max = max( l_max, FP%fisher_cls%l_max_TT, &
            & FP%fisher_cls%l_max_EE, FP%fisher_cls%l_max_BB )

        ! allocation:
        allocate(cl_fiducial(cl_dim, cl_dim, l_max-l_min+1), stat = AllocateStatus)
        if (AllocateStatus /= 0) stop "Allocation failed. Not enough memory to allocate cl_fiducial"
        allocate(cl_temp(cl_dim, cl_dim, l_max-l_min+1), stat = AllocateStatus)
        if (AllocateStatus /= 0) stop "Allocation failed. Not enough memory to allocate cl_temp"
        allocate(cl_derivative(cl_dim, cl_dim, l_max-l_min+1), stat = AllocateStatus)
        if (AllocateStatus /= 0) stop "Allocation failed. Not enough memory to allocate cl_derivative"

        ! computation of the fiducial model:
        time_1 = omp_get_wtime()
        call cl_covariance_out( P, FP, cl_dim, l_min, l_max, cl_fiducial, err )
        time_2 = omp_get_wtime() - time_1

        if ( present(outroot) ) call save_cl_covariance_to_file( cl_fiducial, cl_dim, l_min, l_max, filename=TRIM(outroot)//'fiducial.dat' )

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
            write(*,'(a,i10)')         'Dimension of the Cls matrix  : ', cl_dim
            if ( FP%fisher_cls%Fisher_want_XC ) then
                write(*,'(a)')         'Using cross correlation.'
            else
                write(*,'(a)')         'Not using cross correlation.'
            end if
            write(*,'(a,i10)')         'Minimum l                    : ', l_min
            write(*,'(a,i10)')         'Maximum l                    : ', l_max
            write(*,'(a,i10)')         'Number of parameters         : ', num_param
            write(*,'(a,f10.3,a)')     'Time taken for Cl calculation: ', time_2, ' (s)'
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
            write(*,'(a)') 'Computation of the Cls derivatives'
            write(*,'(a)') '**************************************************************'
        end if

        allocate(cl_derivative_array(num_param,cl_dim, cl_dim, l_max-l_min+1), stat = AllocateStatus)
        if (AllocateStatus /= 0) stop "Allocation failed. Not enough memory to allocate cl_derivative_array"

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

            call cls_derivative( P, FP, ind, &
                & cl_fiducial, cl_derivative, cl_temp, &
                & initial_step, &
                & cl_dim, l_min, l_max, num_param, &
                & err )

            if (err /= 0) print*, 'WARNING: something went wrong with the Cls derivative of parameter:', ind

            if ( present(outroot) .and. FP%cosmicfish_feedback >= 2 ) then
                call save_cl_covariance_to_file( cl_derivative, cl_dim, l_min, l_max, filename=TRIM(outroot)//'derivative_'//trim(adjustl(integer_to_string(ind)))//'.dat' )
            end if

            cl_derivative_array(ind,:,:,:) = cl_derivative(:,:,:)

        end do
        time_2 = omp_get_wtime() - time_1

        if ( FP%cosmicfish_feedback >= 1 ) then
            write(*,*)
            write(*,'(a,f8.3,a)') 'Time to compute all the derivatives:', time_2, ' (s)'
        end if

        ! add the noise to the fiducial model:
        if ( FP%cosmicfish_feedback >= 1 ) then
            write(*,'(a)') '**************************************************************'
            write(*,'(a)') 'Adding experimental noise and computing Fisher:'
            write(*,'(a)') '**************************************************************'
        end if

        time_1 = omp_get_wtime()
        ! save the fiducial
        cl_temp = cl_fiducial
        ! add noise to the Cls
        call add_noise_cls_to_fiducial( P, FP, cl_fiducial, cl_dim, l_min, l_max, err)
        if ( err /= 0 ) stop 'Something went wrong with the computation of the cls noise'

        time_2 = omp_get_wtime() - time_1

        if ( present(outroot) .and. FP%cosmicfish_feedback >= 2 ) then
            call save_cl_covariance_to_file( cl_fiducial, cl_dim, l_min, l_max, filename=TRIM(outroot)//'fid_noise.dat' )
        end if

        if ( FP%cosmicfish_feedback >= 1 ) then
            write(*,'(a,f8.3,a)') 'Time to compute the noise Cls:', time_2, ' (s)'
        end if

        ! compute the Fisher matrix:
        allocate( dummy_matrix_1(cl_dim, cl_dim), dummy_matrix_2(cl_dim, cl_dim) )
        allocate( dummy_matrix_3(cl_dim, cl_dim), dummy_matrix_4(cl_dim, cl_dim) )
        allocate( dummy_matrix_5(cl_dim, cl_dim) )

        ! invert the fiducial + noise:
        do ind = 1, l_max-l_min+1

            forall ( k = 1:cl_dim, l = 1:cl_dim )
                dummy_matrix_1(k,l) = cl_fiducial(k,l, ind)
            end forall

            call Sym_Matrix_Inversion( dummy_matrix_1, err )
            if ( err /= 0 ) then
                write(*,*) 'Matrix inversion failed. Calculation cannot proceed. Error at l= ', ind
                stop
            end if

            forall ( k = 1:cl_dim, l = 1:cl_dim )
                cl_fiducial(k,l, ind) = dummy_matrix_1(k,l)
            end forall

        end do

        ! compute the Fisher:
        do ind = 1, num_param
            do ind2 = 1, num_param

                temp2 = 0._dl

                do i = 1, l_max-l_min+1
                    ! Cosmic variance term:
                    temp = REAL(i+l_min) -0.5_dl
                    ! Initialize all the matrix to use:
                    forall ( j = 1:cl_dim, k = 1:cl_dim )
                        dummy_matrix_1(j, k) = cl_derivative_array(ind, j, k, i)
                        dummy_matrix_2(j, k) = cl_derivative_array(ind2, j, k, i)
                        dummy_matrix_3(j, k) = cl_fiducial(j, k, i)
                    end forall
                    ! do the products:
                    dummy_matrix_1 = MATMUL(dummy_matrix_1, dummy_matrix_3)
                    dummy_matrix_1 = MATMUL(dummy_matrix_3, dummy_matrix_1)
                    dummy_matrix_4 = MATMUL(dummy_matrix_2, dummy_matrix_1)
                    ! take the trace of the result:
                    temp1 = 0._dl
                    do j = 1, cl_dim
                        temp1 = temp1 + dummy_matrix_4(j,j)
                    end do
                    ! sum to build up the Fisher matrix entry:
                    temp2 = temp2 + temp*temp1

                end do

                Fisher_matrix(ind, ind2) = temp2

            end do
        end do
        time_2 = omp_get_wtime() - time_1

        if ( FP%cosmicfish_feedback >= 1 ) then
            write(*,'(a,f8.3,a)') 'Time to compute the Fisher matrix:', time_2, ' (s)'
        end if

    end subroutine Fisher_Cls

    ! ---------------------------------------------------------------------------------------------
    !> This subroutine prints to file the fiducial cls.
    subroutine Fiducial_Cls( P, FP, outroot )

        implicit none

        Type(CAMBparams)        , intent(in) :: P       !< Input CAMBparams
        Type(cosmicfish_params) , intent(in) :: FP      !< Input Cosmicfish params
        character(len=*)        , intent(in) :: outroot !< Input filename that will be used to dump to file the fiducial Cls.

        real(dl) :: time_1, time_2
        !real(dl) :: initial_step, temp, temp1, temp2
        integer  :: cl_dim, l_min, l_max, err, AllocateStatus !, ind, ind2
        !integer  :: i, j, k, l
        !real(dl), allocatable :: param_array(:)
        real(dl), dimension(:,:,:), allocatable :: cl_fiducial

        ! get the size of the cls covariance matrix:
        call dimension_cl_covariance( P, FP, cl_dim )

        ! protect against errors:
        if ( cl_dim == 0 ) then
            write(*,*) 'WARNING: the Cls matrix is of zero dimension'
            return
        end if

        ! decide l_min and l_max:
        l_min = 2
        l_max = maxval( FP%fisher_cls%LSS_lmax, FP%fisher_cls%LSS_number_windows )
        l_max = max( l_max, FP%fisher_cls%l_max_TT, &
            & FP%fisher_cls%l_max_EE, FP%fisher_cls%l_max_BB )

        ! allocation:
        allocate(cl_fiducial(cl_dim, cl_dim, l_max-l_min+1), stat = AllocateStatus)
        if (AllocateStatus /= 0) stop "Allocation failed. Not enough memory to allocate cl_fiducial"

        ! computation of the fiducial model:
        time_1 = omp_get_wtime()
        call cl_covariance_out( P, FP, cl_dim, l_min, l_max, cl_fiducial, err )
        time_2 = omp_get_wtime() - time_1

        ! save to file:
        call save_cl_covariance_to_file( cl_fiducial, cl_dim, l_min, l_max, filename=TRIM(outroot)//'fiducial.dat' )

        if ( err == 0 ) then

        else if ( err == 4 ) then
            write(*,*) 'Fiducial model unstable.'
            stop
        else
            write(*,*) 'Error in computing the fiducial.'
            stop
        end if

        ! print some feedback
        if ( FP%cosmicfish_feedback >= 1 ) then
            write(*,'(a)')         '**************************************************************'
            write(*,'(a,i10)')         'Dimension of the Cls matrix  : ', cl_dim
            write(*,'(a,i10)')         'Minimum l                    : ', l_min
            write(*,'(a,i10)')         'Maximum l                    : ', l_max
            write(*,'(a,f10.3,a)')     'Time taken for Cl calculation: ', time_2, ' (s)'
            write(*,'(a,i10)')         'Number of omp processes      : ', OMP_GET_MAX_THREADS()
            write(*,'(a)')         '**************************************************************'
        end if

    end subroutine Fiducial_Cls

    ! ---------------------------------------------------------------------------------------------

end module Fisher_calculator_Cls

!----------------------------------------------------------------------------------------
