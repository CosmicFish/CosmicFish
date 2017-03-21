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

!> @file 000_Utilities.f90
!! This file contains several utilities.

!----------------------------------------------------------------------------------------
!> This module contains several general purpose utilities.

!> @author Marco Raveri

module cosmicfish_utilities

    use precision

    implicit none

    private

    public integer_to_string, CosmicFish_write_header, dierf

contains

    ! ---------------------------------------------------------------------------------------------
    !> This function converts an integer to a string. Usefull for numbered files output.
    function integer_to_string( number )

        implicit none

        integer, intent(in) :: number               !< Input integer number
        character(10)       :: integer_to_string    !< Output string with the number

        write( integer_to_string, '(i10)' ) number

        integer_to_string = trim(adjustl( integer_to_string ))

    end function integer_to_string

    ! ---------------------------------------------------------------------------------------------
    !> This subroutine writes to the terminal the CosmicFish header. To be called at the beginning
    !! of the applications.
    subroutine CosmicFish_write_header( name )

        implicit none

        character(len=*) :: name

        write(*,'(a)') "**************************************************************"
        write(*,'(a)') "   _____               _     _____     __  "
        write(*,'(a)') "  / ___/__  ___ __ _  (_)___/ __(_)__ / /  "
        write(*,'(a)') " / /__/ _ \(_-</  ' \/ / __/ _// (_-</ _ \ "
        write(*,'(a)') " \___/\___/___/_/_/_/_/\__/_/ /_/___/_//_/ "
        write(*,'(a)') " "
        write(*,'(a)') "**************************************************************"
        write(*,'(a)') name
        write(*,'(a)') " This application was developed using the CosmicFish code."
        write(*,'(a)') "**************************************************************"

    end subroutine CosmicFish_write_header

    ! ---------------------------------------------------------------------------------------------
    !> This function computes the inverse error function in double precision.
    !> Taken from: http://www.kurims.kyoto-u.ac.jp/~ooura/gamerf.html and updated to F90.
    function dierf(y)

        implicit none

        real(dl), intent(in) :: y !< Input value
        real(dl) :: dierf

        real(dl) :: s,t,u,w,x,z,y_temp

        real(dl), parameter ::  &
            &    qa = 9.16461398268964d-01, &
            &    qb = 2.31729200323405d-01, &
            &    qc = 4.88826640273108d-01, &
            &    qd = 1.24610454613712d-01, &
            &    q0 = 4.99999303439796d-01, &
            &    q1 = 1.16065025341614d-01, &
            &    q2 = 1.50689047360223d-01, &
            &    q3 = 2.69999308670029d-01, &
            &    q4 = -7.28846765585675d-02

        real(dl), parameter :: &
            &    pa = 3.97886080735226000d+00, &
            &    pb = 1.20782237635245222d-01, &
            &    p0 = 2.44044510593190935d-01, &
            &    p1 = 4.34397492331430115d-01, &
            &    p2 = 6.86265948274097816d-01, &
            &    p3 = 9.56464974744799006d-01, &
            &    p4 = 1.16374581931560831d+00, &
            &    p5 = 1.21448730779995237d+00, &
            &    p6 = 1.05375024970847138d+00, &
            &    p7 = 7.13657635868730364d-01, &
            &    p8 = 3.16847638520135944d-01, &
            &    p9 = 1.47297938331485121d-02, &
            &    p10 = -1.05872177941595488d-01, &
            &    p11 = -7.43424357241784861d-02

        real(dl), parameter :: &
            &    p12 = 2.20995927012179067d-03, &
            &    p13 = 3.46494207789099922d-02, &
            &    p14 = 1.42961988697898018d-02, &
            &    p15 = -1.18598117047771104d-02, &
            &    p16 = -1.12749169332504870d-02, &
            &    p17 = 3.39721910367775861d-03, &
            &    p18 = 6.85649426074558612d-03, &
            &    p19 = -7.71708358954120939d-04, &
            &    p20 = -3.51287146129100025d-03, &
            &    p21 = 1.05739299623423047d-04, &
            &    p22 = 1.12648096188977922d-03

        ! digest the input:

        y_temp = 1._dl - y

        z = y_temp
        if (y_temp .gt. 1) z = 2 - y_temp
        w = qa - log(z)
        u = sqrt(w)
        s = (qc + log(u)) / w
        t = 1 / (u + qb)
        x = u * (1 - s * (0.5d0 + s * qd)) - &
            &    ((((q4 * t + q3) * t + q2) * t + q1) * t + q0) * t
        t = pa / (pa + x)
        u = t - 0.5d0
        s = (((((((((p22 * u + p21) * u + p20) * u + &
            &    p19) * u + p18) * u + p17) * u + p16) * u + &
            &    p15) * u + p14) * u + p13) * u + p12
        s = ((((((((((((s * u + p11) * u + p10) * u + &
            &    p9) * u + p8) * u + p7) * u + p6) * u + p5) * u + &
            &    p4) * u + p3) * u + p2) * u + p1) * u + p0) * t - &
            &    z * exp(x * x - pb)
        x = x + s * (1 + x * s)
        if (y_temp .gt. 1) x = -x
        dierf = x

    end function dierf

end module cosmicfish_utilities
