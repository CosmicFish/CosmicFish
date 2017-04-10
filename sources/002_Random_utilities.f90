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

!> @file 002_Random_utilities.f90
!! This file contains several subroutine and functions that can be used to
!! generate random numbers from different distributions.


!----------------------------------------------------------------------------------------
!> This module contains the subroutine and functions to generate random numbers from
!! different distributions.
!! This file contains part of code from the CosmoMC code with some additions and some modifications.
!! The part coming from the CosmoMC code is copyrighted by its author.

!> @author Antony Lewis (CosmoMC part), Marco Raveri

module Random_utilities

    use precision

    implicit none

    integer :: rand_inst = 0
    integer, parameter :: krand = KIND(1.d0)
    integer :: Rand_Feedback = 1

contains

    ! ---------------------------------------------------------------------------------------------

    subroutine initRandom(i, i2)
        implicit none
        integer, optional, intent(in) :: i
        integer, optional, intent(in) :: i2
        integer seed_in,kl,ij
        character(len=10) :: fred
        real(krand) :: klr

        if (present(i)) then
            seed_in = i
        else
            seed_in = -1
        end if
        if (seed_in /=-1) then
            if (present(i2)) then
                kl=i2
                if (i2 > 30081) stop 'initRandom:second seed too large'
            else
                kl = 9373
            end if
            ij = i
        else
            call system_clock(count=ij)
            ij = mod(ij + rand_inst*100, 31328)
            call date_and_time(time=fred)
            read (fred,'(e10.3)') klr
            kl = mod(int(klr*1000), 30081)
        end if

        if (Rand_Feedback > 0 ) write(*,'(" Random seeds:",1I6,",",1I6," rand_inst:",1I4)') ij,kl,rand_inst
        call rmarin(ij,kl)

    end subroutine initRandom

    ! ---------------------------------------------------------------------------------------------
    !> This function returns a random number generated from the uniform distribution in [a,b(.
    !! To improve performances no check is done in order to ensure that the random number generator
    !! is initialised.
    real(dl) function random_uniform(a,b)

        implicit none

        real(dl), intent(in) :: a !< lower bound of the interval
        real(dl), intent(in) :: b !< upper bound on the interval


        random_uniform = (b-a)*RANMAR() + a

    end function random_uniform

    ! ---------------------------------------------------------------------------------------------
    !> This function returns a random number generated from the gaussian distribution with mean zero
    !! and variance one. It works with the Box-Muller transformation.
    !! To improve performances no check is done in order to ensure that the random number generator
    !! is initialised.
    real(dl) function gaussian_1()

        implicit none

        real(dl) :: R, V1, V2, FAC
        integer , save :: iset = 0
        real(dl), save :: gset

        if ( iset==0 ) then
            R = 2._dl
            do while (R >= 1._dl)
                V1 = 2.d0*ranmar()-1.d0
                V2 = 2.d0*ranmar()-1.d0
                R  = V1**2+V2**2
            end do
            FAC        = sqrt(-2.d0*log(R)/R)
            GSET       = V1*FAC
            gaussian_1 = V2*FAC
            iset = 1
        else
            gaussian_1 = GSET
            iset = 0
        endif

    end function gaussian_1

    ! ---------------------------------------------------------------------------------------------
    !> This function returns a random number generated from the gaussian distribution with mean mu
    !! and variance sigma.
    !! To improve performances no check is done in order to ensure that the random number generator
    !! is initialised.
    real(dl)  function random_gaussian( mu, sigma )

        implicit none

        real(dl), intent(in) :: mu      !< the mean of the distribution
        real(dl), intent(in) :: sigma   !< the variance of the distribution

        random_gaussian = gaussian_1()*sigma + mu

    end function random_gaussian


    ! This random number generator originally appeared in ''Toward a Universal
    ! Random Number Generator'' by George Marsaglia and Arif Zaman.
    ! Florida State University Report: FSU-SCRI-87-50 (1987)
    !
    ! It was later modified by F. James and published in ''A Review of Pseudo-
    ! random Number Generators''
    !
    ! THIS IS THE BEST KNOWN RANDOM NUMBER GENERATOR AVAILABLE.
    !    (However, a newly discovered technique can yield
    !        a period of 10^600. But that is still in the development stage.)
    !
    ! It passes ALL of the tests for random number generators and has a period
    !   of 2^144, is completely portable (gives bit identical results on all
    !   machines with at least 24-bit mantissas in the floating point
    !   representation).
    !
    ! The algorithm is a combination of a Fibonacci sequence (with lags of 97
    !   and 33, and operation "subtraction plus one, modulo one") and an
    !   "arithmetic sequence" (using subtraction).
    !
    ! On a Vax 11/780, this random number generator can produce a number in
    !    13 microseconds.
    !========================================================================
    !
    !      PROGRAM TstRAN
    !     INTEGER IJ, KL, I
    ! Thee are the seeds needed to produce the test case results
    !      IJ = 1802
    !      KL = 9373
    !
    !
    ! Do the initialization
    !      call rmarin(ij,kl)
    !
    ! Generate 20000 random numbers
    !      do 10 I = 1, 20000
    !         x = RANMAR()
    !10    continue
    !
    ! If the random number generator is working properly, the next six random
    !    numbers should be:
    !          6533892.0  14220222.0  7275067.0
    !    6172232.0  8354498.0   10633180.0
    !
    !
    !
    !      write(6,20) (4096.0*4096.0*RANMAR(), I=1,6)
    !20    format (3f12.1)
    !      end
    !
    subroutine RMARIN(IJ,KL)
        ! This is the initialization routine for the random number generator RANMAR()
        ! NOTE: The seed variables can have values between:    0 <= IJ <= 31328
        !                                                      0 <= KL <= 30081
        !The random number sequences created by these two seeds are of sufficient
        ! length to complete an entire calculation with. For example, if sveral
        ! different groups are working on different parts of the same calculation,
        ! each group could be assigned its own IJ seed. This would leave each group
        ! with 30000 choices for the second seed. That is to say, this random
        ! number generator can create 900 million different subsequences -- with
        ! each subsequence having a length of approximately 10^30.
        !
        ! Use IJ = 1802 & KL = 9373 to test the random number generator. The
        ! subroutine RANMAR should be used to generate 20000 random numbers.
        ! Then display the next six random numbers generated multiplied by 4096*4096
        ! If the random number generator is working properly, the random numbers
        !    should be:
        !           6533892.0  14220222.0  7275067.0
        !           6172232.0  8354498.0   10633180.0
        double precision U(97), C, CD, CM, S, T
        integer I97, J97,i,j,k,l,m
        integer ij,kl
        integer ii,jj


        !      INTEGER IRM(103)

        common /RASET1/ U, C, CD, CM, I97, J97
        if( IJ < 0  .or.  IJ > 31328  .or. &
            KL < 0  .or.  KL > 30081 ) then
            print '(A)', ' The first random number seed must have a value  between 0 and 31328'
            print '(A)',' The second seed must have a value between 0 and   30081'
            stop
        endif
        I = mod(IJ/177, 177) + 2
        J = mod(IJ    , 177) + 2
        K = mod(KL/169, 178) + 1
        L = mod(KL,     169)
        do II = 1, 97
            S = 0.0
            T = 0.5
            do JJ = 1, 24
                M = mod(mod(I*J, 179)*K, 179)
                I = J
                J = K
                K = M
                L = mod(53*L+1, 169)
                if (mod(L*M, 64) >= 32) then
                    S = S + T
                endif
                T = 0.5 * T
3           continue
            end do
            U(II) = S
2       continue
        end do
        C = 362436.0 / 16777216.0
        CD = 7654321.0 / 16777216.0
        CM = 16777213.0 /16777216.0
        I97 = 97
        J97 = 33

    end subroutine RMARIN

    ! ---------------------------------------------------------------------------------------------
    !> This function returns a random number generated from the uniform distribution in [0,1(
    !! This is the random number generator proposed by George Marsaglia in
    !! Florida State University Report: FSU-SCRI-87-50
    !! It was slightly modified by F. James to produce an array of pseudorandom
    !! numbers.
    double precision function RANMAR()

        double precision U(97), C, CD, CM
        integer I97, J97
        double precision uni

        common /RASET1/ U, C, CD, CM, I97, J97
        !      INTEGER IVEC
        UNI = U(I97) - U(J97)
        if( UNI < 0.0 ) UNI = UNI + 1.0
        U(I97) = UNI
        I97 = I97 - 1
        if(I97 == 0) I97 = 97
        J97 = J97 - 1
        if(J97 == 0) J97 = 97
        C = C - CD
        if( C < 0.d0 ) C = C + CM
        UNI = UNI - C
        if( UNI < 0.d0 ) UNI = UNI + 1.0
        RANMAR = UNI

    end function RANMAR

    ! ---------------------------------------------------------------------------------------------

end module Random_utilities

!----------------------------------------------------------------------------------------
