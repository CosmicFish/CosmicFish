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

!> @file 001_Tensor_derivatives_2pt.f90
!! This file contains the two point tensor derivatives algorithm used later by Fisher matrix calculators.

!----------------------------------------------------------------------------------------
!> This module contains the two point tensor derivatives algorithm used later by Fisher matrix calculators.

!> @author Marco Raveri

module Tensor_Derivatives_2pt

    use precision

    implicit none

    private

    ! Interface block with the tensor functions:
    interface
        ! function R->R^0
        subroutine subroutine_scalar( value, error, result )
            use precision
            implicit none
            integer , intent(in ) :: value
            integer , intent(out) :: error
            real(dl), intent(out) :: result
        end subroutine  subroutine_scalar
        ! function R->R^1
        subroutine subroutine_1D( value, error, result )
            use precision
            implicit none
            integer , intent(in ) :: value
            integer , intent(out) :: error
            real(dl), intent(out), dimension(:) :: result
        end subroutine  subroutine_1D
        ! function R->R^2
        subroutine subroutine_2D( value, error, result )
            use precision
            implicit none
            integer , intent(in ) :: value
            integer , intent(out) :: error
            real(dl), intent(out), dimension(:,:) :: result
        end subroutine  subroutine_2D
        ! function R->R^3
        subroutine subroutine_3D( value, error, result )
            use precision
            implicit none
            integer , intent(in ) :: value
            integer , intent(out) :: error
            real(dl), intent(out), dimension(:,:,:) :: result
        end subroutine  subroutine_3D
        ! function R->R^4
        subroutine subroutine_4D( value, error, result )
            use precision
            implicit none
            integer , intent(in ) :: value
            integer , intent(out) :: error
            real(dl), intent(out), dimension(:,:,:,:) :: result
        end subroutine  subroutine_4D
    end interface

contains

    ! ---------------------------------------------------------------------------------------------
    !> This subroutine computes the fixed stepsize derivative of a scalar function.
    subroutine Tensor_Derivatives_2pt_scalar( function, stepsize, output, error )

        implicit none

        procedure(subroutine_scalar) :: function !< Input function to compute the derivative.
        real(dl), intent(in )        :: stepsize !< Input value of the stepsize to use.
        real(dl), intent(out)        :: output   !< Output value of the derivative.
        integer , intent(out)        :: error    !< Output error code. error=0 everything is fine; error=1 calculation failed.

    end subroutine Tensor_Derivatives_2pt_scalar

    !----------------------------------------------------------------------------------------

end module Tensor_Derivatives_2pt

!----------------------------------------------------------------------------------------
