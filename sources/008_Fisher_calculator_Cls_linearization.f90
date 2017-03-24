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

!> @file 008_Fisher_calculator_Cls_linearization.f90
!! This file contains the tools to compute data-set linearization in different ways.

!----------------------------------------------------------------------------------------
!> This module contains the tools to compute data-set linearization in different ways.

!> @author Marco Raveri

module Fisher_calculator_Cls_linearization

    use precision
    use ModelParams
    use NonLinear
    use Transfer
    use omp_lib
    use cosmicfish_utilities
    use cosmicfish_types
    use Fisher_calculator_Cls

    implicit none

    private

    public compute_l_cut_cls, halofit_linearization, precision_settings_l_cut, &
        & apply_l_cut, print_linearization

contains

    ! ---------------------------------------------------------------------------------------------
    !> This subroutine computes the l cut based on the difference between two different Cls.
    !! The first one is usad as a reference and the second one is used to estabish the cut.
    !! The cut is determined by the Cosmic Variance weigthed difference with relative and absolute
    !! tollerances.
    subroutine compute_l_cut_cls( Cls_reference, Cls_2, rtol, atol, cl_dim, l_min, l_max, l_cut )

        implicit none

        real(dl), intent(in)     :: rtol     !< Input relative tollerance
        real(dl), intent(in)     :: atol     !< Input absolute tollerance
        integer , intent(in)     :: l_min    !< Input min l computed
        integer , intent(in)     :: l_max    !< Input max l computed
        integer , intent(in)     :: cl_dim   !< Input dimension of the Cl matrix

        real(dl), intent(in)     :: Cls_reference(cl_dim, cl_dim, l_max-l_min+1) !< Input reference Cls
        real(dl), intent(in)     :: Cls_2(cl_dim, cl_dim, l_max-l_min+1)         !< Second Cls to use to get l max

        integer , intent(out)    :: l_cut(cl_dim)    !< Output vector containing the estimated l-space cut

        integer  :: i, l, l_index
        real(dl) :: test_value

        ! initialization:
        l_cut   = l_max

        ! cycle over the diagonal:
        do i = 1, cl_dim
            ! cycle over the Cls matrix:
            do l = l_min, l_max

                ! compute the l index:
                l_index = l-l_min+1

                ! if the Cls is zero we might suspect to have reached the cutoff.
                ! Let's check for that:
                if ( Cls_reference(i,i,l_index)==0._dl ) then
                    ! check that all the subsequent values are zeroes:
                    if ( all( Cls_reference(i,i,l_index:l_max)==0._dl ) ) then
                        l_cut(i) = l_index
                        exit
                    end if
                end if

                ! test the difference between the linear and non linear cls in units of cosmic variance:
                test_value = abs( Cls_reference(i,i,l_index) -Cls_2(i,i,l_index) )
                test_value = test_value/( rtol*sqrt( 2._dl/(2._dl*real(l)+1._dl) )*abs(Cls_reference(i,i,l_index)) +atol )

                ! decide what to do:
                if ( test_value > 1._dl ) then
                    l_cut(i) = l_index
                    exit
                end if

            end do
        end do

    end subroutine compute_l_cut_cls

    ! ---------------------------------------------------------------------------------------------
    !> This subroutine computes the l-space data linearization using halofit.
    !! It compares the fiducial model with and without halofit and establishes where to cut
    !! based on a given relative and absoulte tollerance.
    subroutine halofit_linearization( P, FP, rtol, atol, cl_dim, l_cut, err, outroot )

        implicit none

        Type(CAMBparams)         :: P                !< Input CAMBparams
        Type(cosmicfish_params)  :: FP               !< Input Cosmicfish params
        real(dl), intent(in)     :: rtol             !< Input relative tollerance
        real(dl), intent(in)     :: atol             !< Input absolute tollerance
        integer , intent(in)     :: cl_dim           !< Input dimension of the Cl matrix

        character(len=*), intent(in), optional :: outroot   !< Optional input: filename that will be used to dump to file optional output if feedback is greater than 1.

        integer , intent(out)    :: l_cut(cl_dim)    !< Output vector containing the estimated l-space cut
        integer , intent(out)    :: err              !< Output error code:
                                                     !< 0 = all fine
                                                     !< 1 = error in input
                                                     !< 2 = error in computation
                                                     !< 3 = error in quality
                                                     !< 4 = routine specific

        Type(CAMBparams)         :: temp_P
        Type(cosmicfish_params)  :: temp_FP
        integer                  :: l_min, l_max, l_index, i, l, AllocateStatus
        real(dl)                 :: test_value, time_1, time_2
        real(dl), dimension(:,:,:), allocatable :: cl_linear_fiducial, cl_non_linear_fiducial

        ! check the input:
        if ( rtol<0._dl .or. atol<0._dl .or. cl_dim<0 ) then
            err = 1
            return
        end if

        ! initialization:
        l_cut   = 0
        err     = 0
        temp_P  = P
        temp_FP = FP

        ! get l_min and l_max:
        l_min = 2
        l_max = maxval( temp_FP%fisher_cls%LSS_lmax, temp_FP%fisher_cls%LSS_number_windows )
        l_max = max( l_max, temp_FP%fisher_cls%l_max_TT, &
            & temp_FP%fisher_cls%l_max_EE, temp_FP%fisher_cls%l_max_BB )

        ! allocation:
        allocate(cl_linear_fiducial(cl_dim, cl_dim, l_max-l_min+1), stat = AllocateStatus)
        if (AllocateStatus /= 0) stop "Allocation failed. Not enough memory to allocate cl_linear_fiducial"
        allocate(cl_non_linear_fiducial(cl_dim, cl_dim, l_max-l_min+1), stat = AllocateStatus)
        if (AllocateStatus /= 0) stop "Allocation failed. Not enough memory to allocate cl_non_linear_fiducial"

        ! get the linear fiducial:
        if ( FP%cosmicfish_feedback >= 2 ) then
            write(*,'(a)') 'Computing the linear fiducial.'
        end if

        temp_P%NonLinear = NonLinear_none

        time_1 = omp_get_wtime()
        call cl_covariance_out( temp_P, temp_FP, cl_dim, l_min, l_max, cl_linear_fiducial    , err )
        time_2 = omp_get_wtime() - time_1
        if ( err/=0 ) return
        if ( present(outroot) ) call save_cl_covariance_to_file( cl_linear_fiducial, cl_dim, l_min, l_max, filename=TRIM(outroot)//'fiducial_linear.dat' )
        if ( FP%cosmicfish_feedback >= 2 ) then
            write(*,'(a,f10.3,a)') 'Time taken for Cl calculation: ', time_2, ' (s)'
        end if

        ! get the non-linear fiducial:
        if ( FP%cosmicfish_feedback >= 2 ) then
            write(*,'(a)') 'Computing the non-linear fiducial.'
        end if

        temp_P%NonLinear = NonLinear_both

        temp_P%transfer%high_precision = .True.
        temp_P%WantTransfer            = .True.
        call Transfer_SetForNonlinearLensing(temp_P%Transfer)
        call Transfer_SortAndIndexRedshifts(temp_P%Transfer)

        time_1 = omp_get_wtime()
        call cl_covariance_out( temp_P, temp_FP, cl_dim, l_min, l_max, cl_non_linear_fiducial, err )
        time_2 = omp_get_wtime() - time_1
        if ( err/=0 ) return
        if ( present(outroot) ) call save_cl_covariance_to_file( cl_non_linear_fiducial, cl_dim, l_min, l_max, filename=TRIM(outroot)//'fiducial_non_linear_halofit.dat' )
        if ( FP%cosmicfish_feedback >= 2 ) then
            write(*,'(a,f10.3,a)') 'Time taken for Cl calculation: ', time_2, ' (s)'
        end if

        ! now get the l cut:
        call compute_l_cut_cls( cl_linear_fiducial, cl_non_linear_fiducial, rtol, atol, cl_dim, l_min, l_max, l_cut )

    end subroutine halofit_linearization

    ! ---------------------------------------------------------------------------------------------
    !> This subroutine studies the numerical stability of the results to precision settings
    !! and finds an l-cut based on the input one.
    subroutine precision_settings_l_cut( P, FP, rtol, atol, cl_dim, l_cut, err, outroot )

        implicit none

        Type(CAMBparams)         :: P                !< Input CAMBparams
        Type(cosmicfish_params)  :: FP               !< Input Cosmicfish params
        real(dl), intent(in)     :: rtol             !< Input relative tollerance
        real(dl), intent(in)     :: atol             !< Input absolute tollerance
        integer , intent(in)     :: cl_dim           !< Input dimension of the Cl matrix

        character(len=*), intent(in), optional :: outroot   !< Optional input: filename that will be used to dump to file optional output if feedback is greater than 1.

        integer , intent(out)    :: l_cut(cl_dim)    !< Output vector containing the estimated l-space cut
        integer , intent(out)    :: err              !< Output error code:
                                                     !< 0 = all fine
                                                     !< 1 = error in input
                                                     !< 2 = error in computation
                                                     !< 3 = error in quality
                                                     !< 4 = routine specific

        Type(CAMBparams)         :: temp_P
        Type(cosmicfish_params)  :: temp_FP
        integer                  :: l_min, l_max, l_index, AllocateStatus
        integer                  :: InputAccuracyBoost, InputlAccuracyBoost, InputlSampleBoost
        real(dl)                 :: time_1, time_2
        real(dl), dimension(:,:,:), allocatable :: cl_1, cl_2

        ! check the input:
        if ( rtol<0._dl .or. atol<0._dl .or. cl_dim<0 ) then
            err = 1
            return
        end if

        ! initialization:
        l_cut   = 0
        err     = 0
        temp_P  = P
        temp_FP = FP

        ! get l_min and l_max:
        l_min = 2
        l_max = maxval( temp_FP%fisher_cls%LSS_lmax, temp_FP%fisher_cls%LSS_number_windows )
        l_max = max( l_max, temp_FP%fisher_cls%l_max_TT, &
            & temp_FP%fisher_cls%l_max_EE, temp_FP%fisher_cls%l_max_BB )

        ! allocation:
        allocate(cl_1(cl_dim, cl_dim, l_max-l_min+1), stat = AllocateStatus)
        if (AllocateStatus /= 0) stop "Allocation failed. Not enough memory to allocate cl_1"
        allocate(cl_2(cl_dim, cl_dim, l_max-l_min+1), stat = AllocateStatus)
        if (AllocateStatus /= 0) stop "Allocation failed. Not enough memory to allocate cl_2"

        ! get the fiducial:
        if ( FP%cosmicfish_feedback >= 2 ) then
            write(*,'(a)') 'Computing the standard precision fiducial'
        end if
        time_1 = omp_get_wtime()
        call cl_covariance_out( temp_P, temp_FP, cl_dim, l_min, l_max, cl_1, err )
        time_2 = omp_get_wtime() - time_1
        if ( err/=0 ) return
        if ( present(outroot) ) call save_cl_covariance_to_file( cl_1, cl_dim, l_min, l_max, filename=TRIM(outroot)//'fiducial_precision.dat' )
        if ( FP%cosmicfish_feedback >= 2 ) then
            write(*,'(a,f10.3,a)') 'Time taken for Cl calculation: ', time_2, ' (s)'
        end if

        ! get the high precision fiducial:
        temp_P%window_kmax_boost = 10._dl*temp_P%window_kmax_boost
        ! save input precision settings:
        InputAccuracyBoost   = AccuracyBoost
        InputlAccuracyBoost  = lAccuracyBoost
        InputlSampleBoost    = lSampleBoost
        ! increase by a factor two precision settings:
        AccuracyBoost   = 2*AccuracyBoost
        lAccuracyBoost  = 2*lAccuracyBoost
        lSampleBoost    = 2*lSampleBoost
        ! do the calculation:
        if ( FP%cosmicfish_feedback >= 2 ) then
            write(*,'(a)') 'Computing the high precision fiducial. This might take a while...'
        end if
        time_1 = omp_get_wtime()
        call cl_covariance_out( temp_P, temp_FP, cl_dim, l_min, l_max, cl_2, err )
        time_2 = omp_get_wtime() - time_1
        if ( err/=0 ) return
        if ( present(outroot) ) call save_cl_covariance_to_file( cl_2, cl_dim, l_min, l_max, filename=TRIM(outroot)//'fiducial_high_precision.dat' )
        if ( FP%cosmicfish_feedback >= 2 ) then
            write(*,'(a,f10.3,a)') 'Time taken for Cl calculation: ', time_2, ' (s)'
        end if

        ! now get the l cut:
        call compute_l_cut_cls( cl_1, cl_2, rtol, atol, cl_dim, l_min, l_max, l_cut )

        ! restore precision settings to normal:
        AccuracyBoost  = InputAccuracyBoost
        lAccuracyBoost = InputlAccuracyBoost
        lSampleBoost   = InputlSampleBoost

    end subroutine precision_settings_l_cut

    ! ---------------------------------------------------------------------------------------------
    !> This subroutine applies the l cut transferring it to parameters.
    subroutine apply_l_cut( P, FP, l_cut )

        implicit none

        Type(CAMBparams)         :: P                !< Input CAMBparams
        Type(cosmicfish_params)  :: FP               !< Input Cosmicfish params
        integer         , intent(in) :: l_cut(:)     !< Input vector containing the estimated l-space cut

        integer :: ind, total_elements, k

        ! check the input:
        total_elements = size(l_cut)

        ! cycle over the l_cut vector and store the results:
        ind = 1

        if ( FP%fisher_cls%Fisher_want_CMB_T ) then
            FP%fisher_cls%l_max_TT = l_cut(ind)
            ind = ind + 1
        end if

        if ( FP%fisher_cls%Fisher_want_CMB_E ) then
            FP%fisher_cls%l_max_EE = l_cut(ind)
            ind = ind + 1
        end if

        if ( FP%fisher_cls%Fisher_want_CMB_lensing ) then
            ind = ind + 1
        end if

        if ( FP%fisher_cls%Fisher_want_CMB_B ) then
            FP%fisher_cls%l_max_BB = l_cut(total_elements)
        end if

        do k = 1, FP%fisher_cls%LSS_number_windows
            FP%fisher_cls%LSS_lmax(k) = l_cut(ind)
            ind = ind + 1
        end do

    end subroutine apply_l_cut

    ! ---------------------------------------------------------------------------------------------
    !> This subroutine prints to screen and file the linearization parameters so that they
    !! can be later used.
    subroutine print_linearization( P, FP, l_cut, outroot )

        implicit none

        Type(CAMBparams)         :: P                !< Input CAMBparams
        Type(cosmicfish_params)  :: FP               !< Input Cosmicfish params
        integer         , intent(in) :: l_cut(:)     !< Input vector containing the estimated l-space cut
        character(len=*), intent(in) :: outroot      !< Input filename that will be used to dump to file

        integer :: ind, total_elements, k

        ! check the input:
        total_elements = size(l_cut)

        ! open the output file:
        open(unit=666, FILE=outroot, ACTION="write", STATUS="replace")

        ! cycle over the l_cut vector:
        ind = 1

        if ( FP%fisher_cls%Fisher_want_CMB_T ) then
            if ( FP%cosmicfish_feedback >= 1 ) then
                write(*,'(a)') 'l_max_TT = '//TRIM(integer_to_string(l_cut(ind)))
            end if
            write(666,'(a)') 'l_max_TT = '//TRIM(integer_to_string(l_cut(ind)))
            ind = ind + 1
        end if

        if ( FP%fisher_cls%Fisher_want_CMB_E ) then
            if ( FP%cosmicfish_feedback >= 1 ) then
                write(*,'(a)') 'l_max_EE = '//TRIM(integer_to_string(l_cut(ind)))
            end if
            write(666,'(a)') 'l_max_EE = '//TRIM(integer_to_string(l_cut(ind)))
            ind = ind + 1
        end if

        if ( FP%fisher_cls%Fisher_want_CMB_lensing ) then
            ind = ind + 1
        end if

        if ( FP%fisher_cls%Fisher_want_CMB_B ) then
            if ( FP%cosmicfish_feedback >= 1 ) then
                write(*,'(a)') 'l_max_BB = '//TRIM(integer_to_string(l_cut(total_elements)))
            end if
            write(666,'(a)') 'l_max_BB = '//TRIM(integer_to_string(l_cut(total_elements)))
        end if

        do k = 1, FP%fisher_cls%LSS_number_windows
            if ( FP%cosmicfish_feedback >= 1 ) then
                write(*,'(a)') 'LSS_lmax('//TRIM(integer_to_string(k))//') = '//TRIM(integer_to_string(l_cut(ind)))
            end if
            write(666,'(a)') 'LSS_lmax('//TRIM(integer_to_string(k))//') = '//TRIM(integer_to_string(l_cut(ind)))
            ind = ind + 1
        end do

        ! close output file:
        close(666)

    end subroutine print_linearization

    ! ---------------------------------------------------------------------------------------------

end module Fisher_calculator_Cls_linearization

!----------------------------------------------------------------------------------------
