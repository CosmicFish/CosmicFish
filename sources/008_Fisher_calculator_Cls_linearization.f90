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
    use cosmicfish_utilities
    use cosmicfish_types
    use Fisher_calculator_Cls

    implicit none

    private

    public halofit_linearization, print_linearization

contains

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
        real(dl)                 :: test_value
        real(dl), dimension(:,:,:), allocatable :: cl_linear_fiducial, cl_non_linear_fiducial

        ! check the input:
        if ( rtol<0._dl ) then
            err = 1
            return
        end if
        if ( atol<0._dl ) then
            err = 1
            return
        end if
        if ( cl_dim<0 ) then
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
        temp_P%NonLinear = NonLinear_none

        call cl_covariance_out( temp_P, temp_FP, cl_dim, l_min, l_max, cl_linear_fiducial    , err )
        if ( err/=0 ) return
        if ( present(outroot) ) call save_cl_covariance_to_file( cl_linear_fiducial, cl_dim, l_min, l_max, filename=TRIM(outroot)//'fiducial_linear.dat' )

        ! get the non-linear fiducial:
        temp_P%NonLinear = NonLinear_both

        temp_P%transfer%high_precision = .True.
        temp_P%WantTransfer            = .True.
        call Transfer_SetForNonlinearLensing(temp_P%Transfer)
        call Transfer_SortAndIndexRedshifts(temp_P%Transfer)

        call cl_covariance_out( temp_P, temp_FP, cl_dim, l_min, l_max, cl_non_linear_fiducial, err )
        if ( err/=0 ) return
        if ( present(outroot) ) call save_cl_covariance_to_file( cl_non_linear_fiducial, cl_dim, l_min, l_max, filename=TRIM(outroot)//'fiducial_non_linear_halofit.dat' )

        ! now get the l cut:
        ! initialize to lmax:
        l_cut   = l_max
        ! cycle over the diagonal:
        do i = 1, cl_dim
            ! cycle over the Cls matrix:
            do l = l_min, l_max

                ! compute the l index:
                l_index = l-l_min+1

                ! if the Cls is zero we might suspect to have reached the cutoff.
                ! Let's check for that:
                if ( cl_linear_fiducial(i,i,l_index)==0._dl ) then
                    ! check that all the subsequent values are zeroes:
                    if ( all( cl_linear_fiducial(i,i,l_index:l_max)==0._dl ) ) then
                        l_cut(i) = l_index
                        exit
                    end if
                end if

                ! test the difference between the linear and non linear cls in units of cosmic variance:
                test_value = abs( cl_linear_fiducial(i,i,l_index) -cl_non_linear_fiducial(i,i,l_index) )
                test_value = test_value/( rtol*sqrt( 2._dl/(2._dl*real(l)+1._dl) )*abs(cl_linear_fiducial(i,i,l_index)) +atol )

                ! decide what to do:
                if ( test_value > 1._dl ) then
                    l_cut(i) = l_index
                    exit
                end if

            end do

        end do

    end subroutine halofit_linearization

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
        open(unit=666, FILE=outroot//'linearization.ini', ACTION="write", STATUS="replace")

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
