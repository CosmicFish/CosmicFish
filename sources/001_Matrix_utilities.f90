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

!> @file 001_Matrix_utilities.f90
!! This file contains several subroutine and functions that can be used to
!! manipulate and modify matrices.
!! For many of them it is simply a wrapper to LAPACK.

!----------------------------------------------------------------------------------------
!> This module contains the subroutine and functions to manipulate and modify matrices.

!> @author Marco Raveri

module Matrix_utilities

    use precision

    implicit none

    private

    public Sym_Matrix_Inversion, matrix_determinant, matrix_eigenvalues_eigenvectors_sym

contains

    ! ---------------------------------------------------------------------------------------------
    !> This subroutine computes the inverse of a square symmetric matrix by means of LAPACK
    !! LU decomposition and inversion.
    subroutine Sym_Matrix_Inversion(A, error_code)

        implicit none

        real(dl), dimension(:,:), intent(inout) :: A !< Input and output matrix. Inversion happends in place.
        integer , intent(out) :: error_code          !< Error code: 0 means everything ok
                                                     !! 1 means singular matrix
                                                     !! 2 means inversion failed for some reason

        real(dl), dimension(size(A,1)) :: work   ! work array for LAPACK
        integer , dimension(size(A,1)) :: ipiv   ! pivot indices
        integer :: n, info

        ! External procedures defined in LAPACK
        external dgetrf
        external dgetri

        error_code = 0
        n = size(A,1)

        ! dgetrf computes an LU factorization of a general M-by-N matrix A
        ! using partial pivoting with row interchanges.
        call dgetrf(n, n, A, n, ipiv, info)

        if (info /= 0) then
            error_code = 1
            write(*,*) 'Matrix is numerically singular!'
            return
        end if

        ! dgetri computes the inverse of a matrix using the LU factorization
        ! computed by dgetrf.
        call dgetri(n, A, n, ipiv, work, n, info)

        if (info /= 0) then
            error_code = 2
            write(*,*) 'Matrix inversion failed!'
            return
        end if

    end subroutine Sym_Matrix_Inversion

    ! ---------------------------------------------------------------------------------------------
    !> This subroutine computes the determinant of a matrix.
    subroutine matrix_determinant(A, det, error_code)

        implicit none

        real(dl), dimension(:,:), intent(in) :: A !< Input matrix.
        real(dl), intent(out) :: det              !< Output determinant of the matrix.
        integer , intent(out) :: error_code       !< Error code: 0 means everything ok

        integer  :: N, i, info, AllocateStatus
        real(dl) :: sgn
        integer , dimension(:), allocatable :: ipiv   ! pivot indices
        ! External procedures defined in LAPACK
        external dgetrf

        N          = size(A,1)
        error_code = 0

        allocate( ipiv(N), stat = AllocateStatus )
        if (AllocateStatus /= 0) stop "matrix_determinant: allocation failed. Not enough memory."

        ipiv=0

        ! dgetrf computes an LU factorization of a general M-by-N matrix A
        ! using partial pivoting with row interchanges.
        call dgetrf(n, n, A, n, ipiv, info)

        if (info /= 0) then
            error_code = 1
            write(*,*) 'Matrix is numerically singular!'
            stop
        end if

        det = 1._dl

        do i=1, N
            det = det*A(i,i)
        end do

        sgn = 1._dl

        do i=1, N
            if(ipiv(i)/=i)then
                sgn=-sgn
            end if
        end do

        det=sgn*det

        return

    end subroutine matrix_determinant

    ! ---------------------------------------------------------------------------------------------
    !> This subroutine computes the eigenvalues and eigenvectors of a matrix.
    subroutine matrix_eigenvalues_eigenvectors_sym(A, eigenvalues, eigenvectors, error_code)

        implicit none

        real(dl), dimension(:,:), intent(in)  :: A             !< Input matrix.
        real(dl), dimension(:)  , intent(out) :: eigenvalues   !< Output eigenvalues
        real(dl), dimension(:,:), intent(out) :: eigenvectors  !< Output eigenvectors
        integer , intent(out) :: error_code                    !< Error code: 0 means everything ok

        eigenvectors = 0._dl
        eigenvalues  = 0._dl
        error_code   = 0

        stop 'matrix_eigenvalues_eigenvectors_sym not yet implemented.'

    end subroutine matrix_eigenvalues_eigenvectors_sym

    ! ---------------------------------------------------------------------------------------------

end module Matrix_utilities

!----------------------------------------------------------------------------------------
