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
    use cosmicfish_utilities
    use cosmicfish_types
    use Fisher_calculator_Cls

    implicit none

    private

    public halofit_linearization

contains

    ! ---------------------------------------------------------------------------------------------
    !> This subroutine computes the l-space data linearization using halofit.
    !! It compares the fiducial model with and without halofit and establishes where to cut
    !! based on a given relative and absoulte tollerance.
    subroutine halofit_linearization( P, FP, rtol, atol, cl_dim, l_cut, err )

        implicit none

        Type(CAMBparams)         :: P                !< Input CAMBparams
        Type(cosmicfish_params)  :: FP               !< Input Cosmicfish params
        real(dl), intent(in)     :: rtol             !< Input relative tollerance
        real(dl), intent(in)     :: atol             !< Input absolute tollerance
        integer , intent(in)     :: cl_dim           !< Input dimension of the Cl matrix

        integer , intent(out)    :: l_cut(cl_dim)    !< Output vector containing the estimated l-space cut
        integer , intent(out)    :: err              !< Output error code:
                                                     !< 0 = all fine
                                                     !< 1 = error in input
                                                     !< 2 = error in computation
                                                     !< 3 = error in quality
                                                     !< 4 = routine specific

        Type(CAMBparams)         :: temp_P
        Type(cosmicfish_params)  :: temp_FP
        integer                  :: l_min, l_max, i, l, AllocateStatus
        real(dl), dimension(:,:,:), allocatable :: cl_linear_fiducial, cl_non_linear_fiducial

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
        temp_P%NonLinear = .False.
        call cl_covariance_out( temp_P, temp_FP, cl_dim, l_min, l_max, cl_linear_fiducial    , err )
        if ( err/=0 ) return

        ! get the non-linear fiducial:
        temp_P%NonLinear = .True.
        call cl_covariance_out( temp_P, temp_FP, cl_dim, l_min, l_max, cl_non_linear_fiducial, err )
        if ( err/=0 ) return

        ! now get the l cut:
        do i = 1, cl_dim
            ! cycle over the Cls matrix:
            do l = 1, l_max-l_min+1




            end do

            l_cut(i) = 0

        end do

    end subroutine halofit_linearization

    ! ---------------------------------------------------------------------------------------------

end module Fisher_calculator_Cls_linearization
