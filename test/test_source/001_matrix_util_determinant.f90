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

program test_matrix_util_determinant

    use precision
    use Matrix_utilities

    implicit none

    real(dl), dimension(2,2) :: A
    real(dl) :: det
    integer  :: error_code

    A(1,1) = 1._dl
    A(1,2) = 2._dl
    A(2,1) = 2._dl
    A(2,2) = 3._dl

    call matrix_determinant(A, det, error_code)

    open (unit = 100, file = "output")

    write(100,*) det

    close(100)

    ! all the successful test should return the following.
    stop 100

end program test_matrix_util_determinant
