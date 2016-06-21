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

program test_Fisher_man_inverr

    use precision
    use Fisher_manipulation

    implicit none

    integer  :: number, i
    real(dl) :: coeff

    open (unit = 100, file = "output")

    number = 100

    do i=1, number

        coeff = 0._dl + 0.99_dl/REAL(number-1)*REAL(i-1)
        write(100,*) confidence_coefficient( coeff )

    end do

    write(100,*)
    write(100,*) confidence_coefficient( 0.68_dl ), confidence_coefficient( 0.95_dl )

    close(100)

    ! all the successful test should return the following.
    stop 100

end program test_Fisher_man_inverr
