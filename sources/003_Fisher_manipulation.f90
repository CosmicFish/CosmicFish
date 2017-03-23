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

!> @file 003_Fisher_manipulation.f90
!! This file contains several subroutine and functions that can be used to
!! manipulate and modify a Fisher matrix.

!----------------------------------------------------------------------------------------
!> This module contains the subroutine and functions to carry out Fisher
!! matrix manipulations.

!> @author Marco Raveri

module Fisher_manipulation

    use precision
    use Matrix_utilities
    use cosmicfish_utilities

    implicit none

    private

    public Fisher_delete_column_rows, confidence_coefficient, confidence_interval_params, &
        & marginalise_Fisher, Fisher_PCA, Fisher_FOM, Fisher_protection_unconstrained,   &
        & Fisher_protection_degenerate

contains

    ! ---------------------------------------------------------------------------------------------
    !> Thie subroutine eliminates the rows and columns in indexes from Fisher_in and returns
    !! a reduced Fisher matrix in Fisher_out
    subroutine Fisher_delete_column_rows( Fisher_in, Fisher_out, indexes )

        implicit none

        real(dl), dimension(:,:), intent(in)  :: Fisher_in  !< Input Fisher matrix
        real(dl), dimension(:,:), intent(out) :: Fisher_out !< Output Fisher matrix
        integer , dimension(:)  , intent(in)  :: indexes    !< indexes of the columns and rows to delete

        Fisher_out = 0._dl

        stop 'Fisher_delete_column_rows not yet implemented'

    end subroutine Fisher_delete_column_rows

    ! ---------------------------------------------------------------------------------------------
    !> This function computes the coefficient correspondent to a given confidence level.
    !! Uses MKL to compute the inverse error function.
    function confidence_coefficient( confidence_level )

        implicit none

        real(dl), intent(in) :: confidence_level    !< Required confidence level
        real(dl) :: confidence_coefficient          !< Output coefficient corresponding to the confidence level

        confidence_coefficient = sqrt(2._dl)*dierf(confidence_level)

    end function confidence_coefficient

    ! ---------------------------------------------------------------------------------------------
    !> This function computes the error bounds on the parameters starting from a Fisher matrix
    subroutine confidence_interval_params( Fisher_in, confidence_level, error_out )

        implicit none

        real(dl), intent(in)   :: Fisher_in(:,:)    !< Input Fisher matrix
        real(dl), intent(in)   :: confidence_level  !< Required confidence level
        real(dl), intent(out)  :: error_out(:)      !< Vector with the confidence bounds on parameters

        error_out = 0._dl

        stop 'confidence_interval_params not yet implemented'

    end subroutine confidence_interval_params

    ! ---------------------------------------------------------------------------------------------
    !> This function marginalises the Fisher_in Fisher matrix over the parameters in index
    !! and returns the marginalized Fisher matrix. The calculation is carried out in the
    !! numerically stable way.
    subroutine marginalise_Fisher( Fisher_in, indexes, Fisher_out )

        implicit none

        real(dl), dimension(:,:), intent(in)  :: Fisher_in  !< Input Fisher matrix
        real(dl), dimension(:,:), intent(out) :: Fisher_out !< Ouput Fisher matrix
        integer , dimension(:)  , intent(in)  :: indexes    !< Indexes to marginalise

        Fisher_out = 0._dl

        stop 'marginalise_Fisher not yet implemented'

    end subroutine marginalise_Fisher

    ! ---------------------------------------------------------------------------------------------
    !> This subroutine performs the principal component analysis on the Fisher matrix Fisher_in
    !! returning the ordered eigenvalues and eigenvectors of the matrix.
    subroutine Fisher_PCA( Fisher_in, eigenvalues, eigenvectors )

        implicit none

        real(dl), dimension(:,:), intent(in)  :: Fisher_in    !< Input Fisher matrix
        real(dl), dimension(:)  , intent(out) :: eigenvalues  !< Output eigenvalues
        real(dl), dimension(:,:), intent(out) :: eigenvectors !< Output eigenvectors

        eigenvalues  = 0._dl
        eigenvectors = 0._dl

        stop 'Fisher_PCA not yet implemented'

    end subroutine Fisher_PCA

    ! ---------------------------------------------------------------------------------------------
    !> This subroutine computes the Figure of Merit from the Fisher matrix Fisher_in
    subroutine Fisher_FOM( Fisher_in, FOM )

        implicit none

        real(dl), dimension(:,:), intent(in)  :: Fisher_in  !< Input Fisher matrix
        real(dl), intent(out) :: FOM                        !< Output Figure of Merit

        FOM = 0._dl

        stop 'Fisher_FOM not yet implemented'

    end subroutine Fisher_FOM

    ! ---------------------------------------------------------------------------------------------
    !> This subroutine applies a safeguard against unconstrained parameters.
    !! If a parameter is unconstrained the corresponding column and row of the Fisher matrix
    !! will be exactly zero. This will make the Fisher matrix not invertible.
    !! This subroutine adds a very small quantity to the diagonal element of the unconstrained
    !! parameter ensuring that the matrix is invertible and that the unconstrained parameter
    !! does not impact the constarints on the others.
    subroutine Fisher_protection_unconstrained( Fisher_in )

        implicit none

        real(dl), dimension(:,:), intent(inout)  :: Fisher_in !< Input and output Fisher matrix

        integer :: fisher_dim, i, j
        logical :: is_zero

        fisher_dim = size( Fisher_in, 1 )

        ! cycle over the Fisher matrix to see wether there are columns or row that are
        ! exactly zero.
        do i=1, fisher_dim

            is_zero = .true.

            do j=1, fisher_dim
                if ( Fisher_in(i,j) /= 0._dl ) is_zero = .false.
            end do

            if ( is_zero ) then
                Fisher_in(i,i) = 1.d-16
            end if

        end do

    end subroutine Fisher_protection_unconstrained

    ! ---------------------------------------------------------------------------------------------
    !> This subroutine applies a safeguard against totally degenerate parameters.
    !! If a parameter is totally degenerate with another the Fisher matrix will be not invertible.
    !! This subroutine adds a very small quantity on the diagonal of the Fisher matrix ensuring
    !! that the resulting matrix will be invertible at a cost of a negligible change in the
    !! parameters bounds.
    subroutine Fisher_protection_degenerate( Fisher_in )

        implicit none

        real(dl), dimension(:,:), intent(inout)  :: Fisher_in !< Input and output Fisher matrix

        integer :: fisher_dim, i

        fisher_dim = size( Fisher_in, 1 )

        ! cycle over the Fisher matrix
        do i=1, fisher_dim
            Fisher_in(i,i) = Fisher_in(i,i) + 1.d-16
        end do

    end subroutine Fisher_protection_degenerate

    ! ---------------------------------------------------------------------------------------------

end module Fisher_manipulation

!----------------------------------------------------------------------------------------
