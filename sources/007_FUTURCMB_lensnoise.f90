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

!> @file 007_FUTURCMB_lensnoise.f90
!! This file contains the relevant code to compute the noise for CMB lensing.
!! This code was developed by Laurence Perotto and Julien Lesgourgues that retain the copyright
!! for the following code.
!! This source file was modified by Marco Raveri to use it with the CosmicFish code.

!----------------------------------------------------------------------------------------
!> This module contains the code to compute the noise for CMB lensing.

module lensnoise
    !+
    !PURPOSE:
    !    calcul of the deflection field noise power spectrum
    !AUTHOR:
    !    J. Lesgourgues
    !REF:
    !    Okamoto & Hu's algorithm
    !    for the calculation of the quadratic estimator variance
    !VERSION:
    !    june 2006: made an indep. module by LP
    !    june 2005: written by JL
    !-

    use AMLutils
    use precision

    implicit none

    private

    public :: calc_Nldd

contains

  ! ---------------------------------------------------------------------------------------------
  !> This subroutine computes the noise for the Okamoto & Hu's estimator of the CMB lensing
  subroutine calc_Nldd (incls, inlcls, clnoise, nldd, jjmaxTmax, lmaxM)
    ! file_cl_plus(1)='../fidO/fidOa5_scalCls.dat'
    ! les Cl lus dans CAMB sont a convertir en muK^2

    ! declaration des arguments :
    REAL(dl),           INTENT(IN), DIMENSION (:,:) :: incls    !< Input scalar Cls vector:
                                                                !< 1 = Temperature
                                                                !< 2 = E mode polarization
                                                                !< 3 = TE cross correlation
                                                                !< 4 = dd lensing deflection angle
                                                                !< 5 = Td temperature lensing deflection cross correlation

    REAL(dl),           INTENT(IN), DIMENSION (:,:) :: inlcls   !< Input lensed Cls vector:
                                                                !< 1 = Temperature
                                                                !< 2 = E mode polarization
                                                                !< 3 = TE cross correlation
                                                                !< 4 = B mode polarization

    REAL(dl),           INTENT(IN), DIMENSION (:,:) :: clnoise  !< Input noise Cls vector:
                                                                !< 1 = Temperature noise
                                                                !< 2 = Polarization noise

    real(dl),           INTENT(OUT), DIMENSION (:) :: nldd      !< Output CMB lensing noise vector

    integer ,           INTENT(INOUT)              :: jjmaxTmax !< FuturCMB internal parameter. Needs to be a suitable value.
                                                                !< If the value is too small the code will enforce compatibility with lmaxM.

    integer ,           INTENT(IN)                 :: lmaxM     !< Maximum l at which to compute lensing noise.

    ! les variables locales
    integer jjmaxT, jj, ll, l1, l2max, nn, l2, l, err
    integer, dimension(:), allocatable :: llens
    double precision all, gall, al1, an, al2
    double precision aAlTT, aAlEE, aAlTE, aAlTB, aAlEB
    double precision aBlTTEE, aBlTTTE, aBlEETE, aBlTBEB
    double precision Al1L0, Al1L2, wigner0, wigner2
    double precision parity_test, parity, prefactor
    double precision fTT12, fEE12, fTE12, fTE21
    double precision fTB12, fTB21, fEB12, fEB21
    double precision gTE12_den, gTE12_num, gTE21_num
    double precision Wl20_sup, Wl22_sup, Wl20_mid, Wl22_mid, Wl20_inf, Wl22_inf
    double precision aa0, bb0, cc0, aa2, bb2, cc2
    double precision, dimension(5,5) ::  LensingNoise, InvLensingNoise
    double precision  GaussVecLens(5,1)
    double precision, dimension(:), allocatable :: MinVarNoise_L
    double precision, dimension(:), allocatable :: ddMinVarNoise_L
    double precision, dimension(:), allocatable :: alogllens
    integer bonalloc
    CHARACTER(LEN=100), PARAMETER :: code = 'calc_NLDD'


    ! debut !
    ! compute noise for lensing

    ! define values for interpolation
    jjmaxT=int((dble(lmaxM)+1700.d0)/50.d0)

    if (jjmaxT.gt.jjmaxTmax) then
        jjmaxTmax = jjmaxT+1
       write(*,*)'WARNING: jjmaxTmax increased to: ',jjmaxT
    end if

    allocate(MinVarNoise_L(1:jjmaxT), &
    &        alogllens(1:jjmaxT), &
    &        llens(1:jjmaxT), stat = bonalloc)
     if (bonalloc > 0) then
        print*,code,'> pb allocation memoire', bonalloc
        stop
     endif
     !allocate(ddMinVarNoise_L(1::jjmaxT))

    llens(1)=2
    llens(2)=3
    llens(3)=4
    llens(4)=5
    llens(5)=6
    llens(6)=7
    llens(7)=8
    llens(8)=9
    do jj=9,17
       llens(jj)=jj*10-80
    end do
    do jj=18,53
       llens(jj)=jj*25-350
    end do
    do jj=54,jjmaxT
       llens(jj)=jj*50-1700
    end do

    ! sum over L

!    open(unit=37,file="test_Nls.dat",form="formatted", &
!         &                       status='replace', iostat=err)
!    if (err /= 0) print *, "Probleme d'ecriture du fichier test..."

    do jj=1,jjmaxT !! do 1

       ll=llens(jj)
       all=dble(ll)

       !! test !!
       !print *, code, incls(1,ll), incls(2,ll), incls(3,ll)
       !print *, code, inlcls(1,ll),inlcls(2,ll),inlcls(3,ll),inlcls(4,ll)
       !print *, code, clnoise(1,ll),clnoise(2,ll)

       ! determine where to cut the sum over l1 in two parts

       gall=all/2.d0
       if (dint(gall).ne.gall) gall=gall+0.5d0
       if (gall.le.2.d0) gall=2.d0

       ! reset the variable accounting for the various sums over l1,l2

       aAlTT=0.d0
       aAlEE=0.d0
       aAlTE=0.d0
       aAlTB=0.d0
       aAlEB=0.d0
       aBlTTEE=0.d0
       aBlTTTE=0.d0
       aBlEETE=0.d0
       aBlTBEB=0.d0

       ! start loop on l1 (first part)

       do l1=2,int(gall)-1 !! do 2
          al1=dble(l1)

          ! determine maximum value of l2 in the sum over l2

          l2max=ll+l1
          ! if (l2max.gt.lmaxTP) l2max=lmaxTP

          !compute A(l1,L,0) (useful for first values of W_l2)

          Al1L0=1.d0

          do nn=1,ll
             an=dble(nn)
             Al1L0=Al1L0*(an-0.5d0)*(al1+an)/an/(al1+an-0.5d0)
          end do

          Al1L0=dsqrt(Al1L0)

          Al1L2=Al1L0 &
          &           *dsqrt(al1*(al1-1.d0)/(al1+1.d0)/(al1+2.d0)) &
          &           *dsqrt((al1+all+1.d0)*(al1+all+2.d0) &
          &                  /(al1+all)/(al1+all-1.d0))

          ! case l2=l2max

          l2=l2max
          al2=dble(l2)

          Wl20_sup=Al1L0/dsqrt(2.d0*(all+al1)+1.d0)
          Wl22_sup=Al1L2/dsqrt(2.d0*(all+al1)+1.d0)

          parity_test=(all+al1)/2.d0

          if (dint(parity_test).ne.parity_test) then
             Wl20_sup=-Wl20_sup
             Wl22_sup=-Wl22_sup
          end if

          parity=(al1+all+al2)/2.d0

          if (l2.le.lmaxM) then

             wigner0=Wl20_sup
             wigner2=Wl22_sup

             prefactor=dsqrt((2.d0*all+1.d0)*(2.d0*al1+1.d0) &
             &              *(2.d0*al2+1.d0)/(16.d0*pi))

             if (dint(parity).eq.parity) then

                fTT12=prefactor*wigner0*(incls(1,l1) &
                &     *(all*(all+1.d0)+al1*(al1+1.d0)-al2*(al2+1.d0)) &
                &     +incls(1,l2) &
                &     *(all*(all+1.d0)+al2*(al2+1.d0)-al1*(al1+1.d0)))
                fEE12=prefactor*wigner2*(incls(2,l1) &
                &     *(all*(all+1.d0)+al1*(al1+1.d0)-al2*(al2+1.d0)) &
                &     +incls(2,l2) &
                &     *(all*(all+1.d0)+al2*(al2+1.d0)-al1*(al1+1.d0)))
                fTE12=prefactor*(wigner2*incls(3,l1) &
                &     *(all*(all+1.d0)+al1*(al1+1.d0)-al2*(al2+1.d0)) &
                &     +wigner0*incls(3,l2) &
                &     *(all*(all+1.d0)+al2*(al2+1.d0)-al1*(al1+1.d0)))
                fTE21=prefactor*(wigner2*incls(2,l2) &
                &     *(all*(all+1.d0)+al2*(al2+1.d0)-al1*(al1+1.d0)) &
                &     +wigner0*incls(3,l1) &
                &     *(all*(all+1.d0)+al1*(al1+1.d0)-al2*(al2+1.d0)))

                gTE12_den=(inlcls(1,l1)+clnoise(1,l1)) &
                &         *(inlcls(1,l2)+clnoise(1,l2)) &
                &         *(inlcls(2,l1)+clnoise(2,l1)) &
                &         *(inlcls(2,l2)+clnoise(2,l2)) &
                &         -(inlcls(3,l1)*inlcls(3,l2))**2
                gTE12_num=(inlcls(1,l2)+clnoise(1,l2)) &
                &         *(inlcls(2,l1)+clnoise(2,l1)) &
                &         *fTE12-(inlcls(3,l1)) &
                &         *(inlcls(3,l2))*fTE21
                gTE21_num=(inlcls(1,l1)+clnoise(1,l1)) &
                &         *(inlcls(2,l2)+clnoise(2,l2)) &
                &         *fTE21-(inlcls(3,l1)) &
                &         *(inlcls(3,l2))*fTE12

                aAlTT = aAlTT +fTT12**2 &
                &  /(inlcls(1,l1)+clnoise(1,l1))/(inlcls(1,l2)+clnoise(1,l2))
                aAlEE = aAlEE +fEE12**2 &
                &  /(inlcls(2,l1)+clnoise(2,l1))/(inlcls(2,l2)+clnoise(2,l2))
                aAlTE = aAlTE + (fTE12*gTE12_num+fTE21*gTE21_num) &
                &              /gTE12_den

                aBlTTEE= aBlTTEE+fTT12*fEE12 &
                &  *(inlcls(3,l1))*(inlcls(3,l2)) &
                &  /(inlcls(1,l1)+clnoise(1,l1))/(inlcls(1,l2)+clnoise(1,l2)) &
                &  /(inlcls(2,l1)+clnoise(2,l1))/(inlcls(2,l2)+clnoise(2,l2))
                aBlTTTE= aBlTTTE+fTT12/gTE12_den &
                &  /(inlcls(1,l1)+clnoise(1,l1))/(inlcls(1,l2)+clnoise(1,l2)) &
                &  *((inlcls(1,l1)+clnoise(1,l1))*(inlcls(3,l2)) &
                &  *gTE12_num &
                &  +(inlcls(1,l2)+clnoise(1,l2))*(inlcls(3,l1)) &
                &  *gTE21_num)
                aBlEETE= aBlEETE+fEE12/gTE12_den &
                &  /(inlcls(2,l1)+clnoise(2,l1))/(inlcls(2,l2)+clnoise(2,l2)) &
                &  *((inlcls(2,l2)+clnoise(2,l2))*(inlcls(3,l1)) &
                &  *gTE12_num &
                &  +(inlcls(2,l1)+clnoise(2,l1))*(inlcls(3,l2)) &
                &  *gTE21_num)

             else

                fTB12=prefactor*wigner2*incls(3,l1) &
                &   *(all*(all+1.d0)+al1*(al1+1.d0)-al2*(al2+1.d0))
                fTB21=prefactor*wigner2*incls(3,l2) &
                &   *(all*(all+1.d0)+al2*(al2+1.d0)-al1*(al1+1.d0))
                fEB12=prefactor*wigner2*incls(2,l1) &
                &   *(all*(all+1.d0)+al1*(al1+1.d0)-al2*(al2+1.d0))
                fEB21=prefactor*wigner2*incls(2,l2) &
                &   *(all*(all+1.d0)+al2*(al2+1.d0)-al1*(al1+1.d0))

                aAlTB = aAlTB  &
                &     + fTB12**2/(inlcls(1,l1)+clnoise(1,l1)) &
                &               /(inlcls(4,l2)+clnoise(2,l2)) &
                &     + fTB21**2/(inlcls(1,l2)+clnoise(1,l2)) &
                &               /(inlcls(4,l1)+clnoise(2,l1))
                aAlEB = aAlEB &
                &     + fEB12**2/(inlcls(2,l1)+clnoise(2,l1)) &
                &               /(inlcls(4,l2)+clnoise(2,l2)) &
                &     + fEB21**2/(inlcls(2,l2)+clnoise(2,l2)) &
                &               /(inlcls(4,l1)+clnoise(2,l1))
                aBlTBEB = aBlTBEB &
                &       + fTB12*fEB12*(inlcls(3,l1)) &
                &              /(inlcls(4,l2)+clnoise(2,l2)) &
                &   /(inlcls(1,l1)+clnoise(1,l1)) &
                &   /(inlcls(2,l1)+clnoise(2,l1)) &
                &       + fTB21*fEB21*(inlcls(3,l2)) &
                &              /(inlcls(4,l1)+clnoise(2,l1)) &
                &   /(inlcls(1,l2)+clnoise(1,l2))/(inlcls(2,l2)+clnoise(2,l2))

             end if
          end if

          ! case l2=l2max-1

          l2=l2max-1
          al2=dble(l2)

          Wl20_mid=0.
          Wl22_mid=Al1L2*2.d0* &
          &       dsqrt(all/al1/(al1+all-2.d0)/(al1+all+2.d0))

          if (dint(parity_test).ne.parity_test) then
             Wl22_mid=-Wl22_mid
          end if

          parity=(al1+all+al2)/2.d0

          if (l2.le.lmaxM) then

             wigner0=Wl20_mid
             wigner2=Wl22_mid

             prefactor=dsqrt((2.d0*all+1.d0)*(2.d0*al1+1.d0) &
             &         *(2.d0*al2+1.d0)/(16.d0*pi))

             if (dint(parity).eq.parity) then
                fTT12=prefactor*wigner0*(incls(1,l1) &
                &     *(all*(all+1.d0)+al1*(al1+1.d0)-al2*(al2+1.d0)) &
                &    +incls(1,l2) &
                &     *(all*(all+1.d0)+al2*(al2+1.d0)-al1*(al1+1.d0)))
                fEE12=prefactor*wigner2*(incls(2,l1) &
                &     *(all*(all+1.d0)+al1*(al1+1.d0)-al2*(al2+1.d0)) &
                &    +incls(2,l2) &
                &     *(all*(all+1.d0)+al2*(al2+1.d0)-al1*(al1+1.d0)))
                fTE12=prefactor*(wigner2*incls(3,l1) &
                &     *(all*(all+1.d0)+al1*(al1+1.d0)-al2*(al2+1.d0)) &
                &    +wigner0*incls(3,l2) &
                &     *(all*(all+1.d0)+al2*(al2+1.d0)-al1*(al1+1.d0)))
                fTE21=prefactor*(wigner2*incls(3,l2) &
                &     *(all*(all+1.d0)+al2*(al2+1.d0)-al1*(al1+1.d0)) &
                &    +wigner0*incls(3,l1) &
                &     *(all*(all+1.d0)+al1*(al1+1.d0)-al2*(al2+1.d0)))

                gTE12_den=(inlcls(1,l1)+clnoise(1,l1)) &
                &         *(inlcls(1,l2)+clnoise(1,l2)) &
                &         *(inlcls(2,l1)+clnoise(2,l1)) &
                &         *(inlcls(2,l2)+clnoise(2,l2)) &
                &        -((inlcls(3,l1)) &
                &         *(inlcls(3,l2)))**2
                gTE12_num=(inlcls(1,l2)+clnoise(1,l2)) &
                &         *(inlcls(2,l1)+clnoise(2,l1))*fTE12-(inlcls(3,l1)) &
                &         *(inlcls(3,l2))*fTE21
                gTE21_num=(inlcls(1,l1)+clnoise(1,l1)) &
                &         *(inlcls(2,l2)+clnoise(2,l2))*fTE21-(inlcls(3,l1)) &
                &         *(inlcls(3,l2))*fTE12

                aAlTT = aAlTT +fTT12**2 &
                &  /(inlcls(1,l1)+clnoise(1,l1))/(inlcls(1,l2)+clnoise(1,l2))
                aAlEE = aAlEE +fEE12**2 &
                &  /(inlcls(2,l1)+clnoise(2,l1))/(inlcls(2,l2)+clnoise(2,l2))
                aAlTE = aAlTE + (fTE12*gTE12_num+fTE21*gTE21_num) &
                &              /gTE12_den

                aBlTTEE= aBlTTEE+fTT12*fEE12 &
                &       *(inlcls(3,l1))*(inlcls(3,l2)) &
                &  /(inlcls(1,l1)+clnoise(1,l1))/(inlcls(1,l2)+clnoise(1,l2)) &
                &  /(inlcls(2,l1)+clnoise(2,l1))/(inlcls(2,l2)+clnoise(2,l2))
                aBlTTTE= aBlTTTE+fTT12/gTE12_den &
                &  /(inlcls(1,l1)+clnoise(1,l1))/(inlcls(1,l2)+clnoise(1,l2)) &
                &       *((inlcls(1,l1)+clnoise(1,l1))*(inlcls(3,l2)) &
                &       *gTE12_num &
                &      +(inlcls(1,l2)+clnoise(1,l2))*(inlcls(3,l1)) &
                &       *gTE21_num)
                aBlEETE= aBlEETE+fEE12/gTE12_den &
                &  /(inlcls(2,l1)+clnoise(2,l1))/(inlcls(2,l2)+clnoise(2,l2)) &
                &        *((inlcls(2,l2)+clnoise(2,l2))*(inlcls(3,l1)) &
                &        *gTE12_num &
                &      +(inlcls(2,l1)+clnoise(2,l1))*(inlcls(3,l2)) &
                &        *gTE21_num)

             else

                fTB12=prefactor*wigner2*incls(3,l1) &
                &     *(all*(all+1.d0)+al1*(al1+1.d0)-al2*(al2+1.d0))
                fTB21=prefactor*wigner2*incls(3,l2) &
                &     *(all*(all+1.d0)+al2*(al2+1.d0)-al1*(al1+1.d0))
                fEB12=prefactor*wigner2*inlcls(2,l1) &
                &     *(all*(all+1.d0)+al1*(al1+1.d0)-al2*(al2+1.d0))
                fEB21=prefactor*wigner2*inlcls(2,l2) &
                &     *(all*(all+1.d0)+al2*(al2+1.d0)-al1*(al1+1.d0))

                aAlTB = aAlTB &
                &     + fTB12**2/(inlcls(1,l1)+clnoise(1,l1)) &
                &               /(inlcls(4,l2)+clnoise(2,l2)) &
                &     + fTB21**2/(inlcls(1,l2)+clnoise(1,l2)) &
                &               /(inlcls(4,l1)+clnoise(2,l1))
                aAlEB = aAlEB &
                &     + fEB12**2/(inlcls(2,l1)+clnoise(2,l1)) &
                &               /(inlcls(4,l2)+clnoise(2,l2)) &
                &     + fEB21**2/(inlcls(2,l2)+clnoise(2,l2)) &
                &               /(inlcls(4,l1)+clnoise(2,l1))

                aBlTBEB = aBlTBEB &
                &       + fTB12*fEB12*(inlcls(3,l1)) &
                &         /(inlcls(4,l2)+clnoise(2,l2)) &
                &         /(inlcls(1,l1)+clnoise(1,l1)) &
                &         /(inlcls(2,l1)+clnoise(2,l1)) &
                &       + fTB21*fEB21*(inlcls(3,l2)) &
                &         /(inlcls(4,l1)+clnoise(2,l1)) &
                &         /(inlcls(1,l2)+clnoise(1,l2)) &
                &         /(inlcls(2,l2)+clnoise(2,l2))

             end if
          end if

          ! sum over other cases, down to l2=L-l1

          do nn=1,l2max-ll+l1-1 !! do 3

             l2=l2max-1-nn
             al2=dble(l2)

             ! recursion formula giving W_l2

             aa0=(al2+2.d0)*(al2+1.d0)* &
             &              dsqrt((-al2+all+al1)*(al2-all+al1+1.d0) &
             &              *(al2+all-al1+1.d0)*(al2+all+al1+2.d0))
             bb0=0.d0
             cc0=-(al2+1.d0)*(al2+2.d0)* &
             &              dsqrt((-al2+all+al1-1.d0)*(al2-all+al1+2.d0) &
             &              *(al2+all-al1+2.d0)*(al2+all+al1+3.d0))

             Wl20_inf=(bb0*Wl20_mid+cc0*Wl20_sup)/aa0

             aa2=(al2+2.d0)*dsqrt(((al2+1.d0)**2-4.d0) &
             &              *(-al2+all+al1)*(al2-all+al1+1.d0) &
             &              *(al2+all-al1+1.d0)*(al2+all+al1+2.d0))
             bb2=2.d0*(2.d0*al2+3.d0)*((al2+1.d0)*(al2+2.d0) &
             &              +all*(all+1.d0)-al1*(al1+1.d0))
             cc2=-(al2+1.d0)*dsqrt(((al2+2.d0)**2-4.d0) &
             &              *(-al2+all+al1-1.d0)*(al2-all+al1+2.d0) &
             &              *(al2+all-al1+2.d0)*(al2+all+al1+3.d0))

             Wl22_inf=(bb2*Wl22_mid+cc2*Wl22_sup)/aa2

             parity=(al1+all+al2)/2.d0

             if (l2.le.lmaxM) then

                wigner0=Wl20_inf
                wigner2=Wl22_inf

                prefactor=dsqrt((2.d0*all+1.d0)*(2.d0*al1+1.d0) &
                &                 *(2.d0*al2+1.d0)/(16.d0*3.14159d0))

                if (dint(parity).eq.parity) then

                   fTT12=prefactor*wigner0*(incls(1,l1) &
                   &    *(all*(all+1.d0)+al1*(al1+1.d0)-al2*(al2+1.d0)) &
                   &    +incls(1,l2) &
                   &    *(all*(all+1.d0)+al2*(al2+1.d0)-al1*(al1+1.d0)))
                   fEE12=prefactor*wigner2*(incls(2,l1) &
                   &    *(all*(all+1.d0)+al1*(al1+1.d0)-al2*(al2+1.d0)) &
                   &    +incls(2,l2) &
                   &    *(all*(all+1.d0)+al2*(al2+1.d0)-al1*(al1+1.d0)))
                   fTE12=prefactor*(wigner2*incls(3,l1) &
                   &    *(all*(all+1.d0)+al1*(al1+1.d0)-al2*(al2+1.d0)) &
                   &    +wigner0*incls(3,l2) &
                   &    *(all*(all+1.d0)+al2*(al2+1.d0)-al1*(al1+1.d0)))
                   fTE21=prefactor*(wigner2*incls(3,l2) &
                   &    *(all*(all+1.d0)+al2*(al2+1.d0)-al1*(al1+1.d0)) &
                   &    +wigner0*incls(3,l1) &
                   &    *(all*(all+1.d0)+al1*(al1+1.d0)-al2*(al2+1.d0)))

                   gTE12_den=(inlcls(1,l1)+clnoise(1,l1)) &
                   &        *(inlcls(1,l2)+clnoise(1,l2)) &
                   &        *(inlcls(2,l1)+clnoise(2,l1)) &
                   &        *(inlcls(2,l2)+clnoise(2,l2)) &
                   &        -((inlcls(3,l1)) &
                   &        *(inlcls(3,l2)))**2
                   gTE12_num=(inlcls(1,l2)+clnoise(1,l2)) &
                   &        *(inlcls(2,l1)+clnoise(2,l1)) &
                   &        *fTE12-(inlcls(3,l1)) &
                   &        *(inlcls(3,l2))*fTE21
                   gTE21_num=(inlcls(1,l1)+clnoise(1,l1)) &
                   &        *(inlcls(2,l2)+clnoise(2,l2)) &
                   &        *fTE21-(inlcls(3,l1)) &
                   &        *(inlcls(3,l2))*fTE12

                   aAlTT = aAlTT +fTT12**2 &
                   &     /(inlcls(1,l1)+clnoise(1,l1)) &
                   &     /(inlcls(1,l2)+clnoise(1,l2))
                   aAlEE = aAlEE +fEE12**2 &
                   &     /(inlcls(2,l1)+clnoise(2,l1)) &
                   &     /(inlcls(2,l2)+clnoise(2,l2))
                   aAlTE = aAlTE + (fTE12*gTE12_num+fTE21*gTE21_num) &
                   &              /gTE12_den

                   aBlTTEE= aBlTTEE+fTT12*fEE12 &
                   &       *(inlcls(3,l1))*(inlcls(3,l2)) &
                   &       /(inlcls(1,l1)+clnoise(1,l1)) &
                   &       /(inlcls(1,l2)+clnoise(1,l2)) &
                   &       /(inlcls(2,l1)+clnoise(2,l1)) &
                   &       /(inlcls(2,l2)+clnoise(2,l2))
                   aBlTTTE= aBlTTTE+fTT12/gTE12_den &
                   &       /(inlcls(1,l1)+clnoise(1,l1)) &
                   &       /(inlcls(1,l2)+clnoise(1,l2)) &
                   &       *((inlcls(1,l1)+clnoise(1,l1))*(inlcls(3,l2)) &
                   &       *gTE12_num &
                   &       +(inlcls(1,l2)+clnoise(1,l2))*(inlcls(3,l1)) &
                   &       *gTE21_num)
                   aBlEETE= aBlEETE+fEE12/gTE12_den &
                   &       /(inlcls(2,l1)+clnoise(2,l1)) &
                   &       /(inlcls(2,l2)+clnoise(2,l2)) &
                   &       *((inlcls(2,l2)+clnoise(2,l2))*(inlcls(3,l1)) &
                   &       *gTE12_num &
                   &       +(inlcls(2,l1)+clnoise(2,l1))*(inlcls(3,l2)) &
                   &       *gTE21_num)

                else

                   fTB12=prefactor*wigner2*incls(3,l1) &
                   &    *(all*(all+1.d0)+al1*(al1+1.d0)-al2*(al2+1.d0))
                   fTB21=prefactor*wigner2*incls(3,l2) &
                   &    *(all*(all+1.d0)+al2*(al2+1.d0)-al1*(al1+1.d0))
                   fEB12=prefactor*wigner2*incls(2,l1) &
                   &    *(all*(all+1.d0)+al1*(al1+1.d0)-al2*(al2+1.d0))
                   fEB21=prefactor*wigner2*incls(2,l2) &
                   &    *(all*(all+1.d0)+al2*(al2+1.d0)-al1*(al1+1.d0))

                   aAlTB = aAlTB  &
                   &     + fTB12**2/(inlcls(1,l1)+clnoise(1,l1)) &
                   &               /(inlcls(4,l2)+clnoise(2,l2)) &
                   &     + fTB21**2/(inlcls(1,l2)+clnoise(1,l2)) &
                   &               /(inlcls(4,l1)+clnoise(2,l1))
                   aAlEB = aAlEB  &
                   &     + fEB12**2/(inlcls(2,l1)+clnoise(2,l1)) &
                   &               /(inlcls(4,l2)+clnoise(2,l2)) &
                   &     + fEB21**2/(inlcls(2,l2)+clnoise(2,l2)) &
                   &               /(inlcls(4,l1)+clnoise(2,l1))

                   aBlTBEB = aBlTBEB &
                   &       + fTB12*fEB12*(inlcls(3,l1)) &
                   &                    /(inlcls(4,l2)+clnoise(2,l2)) &
                   &                    /(inlcls(1,l1)+clnoise(1,l1)) &
                   &                    /(inlcls(2,l1)+clnoise(2,l1)) &
                   &       + fTB21*fEB21*(inlcls(3,l2)) &
                   &                    /(inlcls(4,l1)+clnoise(2,l1)) &
                   &                    /(inlcls(1,l2)+clnoise(1,l2)) &
                   &                    /(inlcls(2,l2)+clnoise(2,l2))

                end if
             end if

             ! shift for the next recursion

             Wl20_sup=Wl20_mid
             Wl20_mid=Wl20_inf

             Wl22_sup=Wl22_mid
             Wl22_mid=Wl22_inf

          end do !! do 3

       end do !! do 2

       ! start loop on l1 (second part)

       do l1=int(gall),lmaxM !! do 2
          al1=dble(l1)

          ! determine maximum value of l2 in the sum over l2

          l2max=ll+l1
          ! if (l2max.gt.lmaxTP) l2max=lmaxTP

          ! compute A(l1,L,0) (useful for first values of W_l2)

          Al1L0=1.d0

          do nn=1,ll
             an=dble(nn)
             Al1L0=Al1L0*(an-0.5d0)*(al1+an)/an/(al1+an-0.5d0)
          end do

          Al1L0=dsqrt(Al1L0)

          Al1L2=Al1L0 &
          &           *dsqrt(al1*(al1-1.d0)/(al1+1.d0)/(al1+2.d0)) &
          &           *dsqrt((al1+all+1.d0)*(al1+all+2.d0) &
          &                  /(al1+all)/(al1+all-1.d0))

          ! case l2=l2max

          l2=l2max
          al2=dble(l2)

          Wl20_sup=Al1L0/dsqrt(2.d0*(all+al1)+1.d0)
          Wl22_sup=Al1L2/dsqrt(2.d0*(all+al1)+1.d0)

          parity_test=(all+al1)/2.d0

          if (dint(parity_test).ne.parity_test) then
             Wl20_sup=-Wl20_sup
             Wl22_sup=-Wl22_sup
          end if

          parity=(al1+all+al2)/2.d0

          if (l2.le.lmaxM) then

             wigner0=Wl20_sup
             wigner2=Wl22_sup

             if (dint(parity).eq.parity) then

                prefactor=dsqrt((2.d0*all+1.d0)*(2.d0*al1+1.d0) &
                &                *(2.d0*al2+1.d0)/(16.d0*pi))

                fTT12=prefactor*wigner0*(incls(1,l1) &
                &    *(all*(all+1.d0)+al1*(al1+1.d0)-al2*(al2+1.d0)) &
                &    +incls(1,l2) &
                &    *(all*(all+1.d0)+al2*(al2+1.d0)-al1*(al1+1.d0)))
                fEE12=prefactor*wigner2*(incls(2,l1) &
                &    *(all*(all+1.d0)+al1*(al1+1.d0)-al2*(al2+1.d0)) &
                &    +incls(2,l2) &
                &    *(all*(all+1.d0)+al2*(al2+1.d0)-al1*(al1+1.d0)))
                fTE12=prefactor*(wigner2*incls(3,l1) &
                &    *(all*(all+1.d0)+al1*(al1+1.d0)-al2*(al2+1.d0)) &
                &    +wigner0*incls(3,l2) &
                &    *(all*(all+1.d0)+al2*(al2+1.d0)-al1*(al1+1.d0)))
                fTE21=prefactor*(wigner2*incls(3,l2) &
                &    *(all*(all+1.d0)+al2*(al2+1.d0)-al1*(al1+1.d0)) &
                &    +wigner0*incls(3,l1) &
                &    *(all*(all+1.d0)+al1*(al1+1.d0)-al2*(al2+1.d0)))

                gTE12_den=(inlcls(1,l1)+clnoise(1,l1)) &
                &         *(inlcls(1,l2)+clnoise(1,l2)) &
                &         *(inlcls(2,l1)+clnoise(2,l1)) &
                &         *(inlcls(2,l2)+clnoise(2,l2)) &
                &        -((inlcls(3,l1)) &
                &         *(inlcls(3,l2)))**2
                gTE12_num=(inlcls(1,l2)+clnoise(1,l2)) &
                &         *(inlcls(2,l1)+clnoise(2,l1))*fTE12-(inlcls(3,l1)) &
                &         *(inlcls(3,l2))*fTE21
                gTE21_num=(inlcls(1,l1)+clnoise(1,l1)) &
                &         *(inlcls(2,l2)+clnoise(2,l2))*fTE21-(inlcls(3,l1)) &
                &         *(inlcls(3,l2))*fTE12

                aAlTT = aAlTT +fTT12**2 &
                &              /(inlcls(1,l1)+clnoise(1,l1)) &
                &              /(inlcls(1,l2)+clnoise(1,l2))
                aAlEE = aAlEE +fEE12**2 &
                &              /(inlcls(2,l1)+clnoise(2,l1)) &
                &              /(inlcls(2,l2)+clnoise(2,l2))
                aAlTE = aAlTE + (fTE12*gTE12_num+fTE21*gTE21_num) &
                &              /gTE12_den

                aBlTTEE= aBlTTEE+fTT12*fEE12 &
                &               *(inlcls(3,l1))*(inlcls(3,l2)) &
                &               /(inlcls(1,l1)+clnoise(1,l1)) &
                &               /(inlcls(1,l2)+clnoise(1,l2)) &
                &               /(inlcls(2,l1)+clnoise(2,l1)) &
                &               /(inlcls(2,l2)+clnoise(2,l2))
                aBlTTTE= aBlTTTE+fTT12/gTE12_den &
                &               /(inlcls(1,l1)+clnoise(1,l1)) &
                &               /(inlcls(1,l2)+clnoise(1,l2)) &
                &               *((inlcls(1,l1)+clnoise(1,l1))*(inlcls(3,l2)) &
                &               *gTE12_num &
                &               +(inlcls(1,l2)+clnoise(1,l2))*(inlcls(3,l1)) &
                &               *gTE21_num)
                aBlEETE= aBlEETE+fEE12/gTE12_den &
                &               /(inlcls(2,l1)+clnoise(2,l1)) &
                &               /(inlcls(2,l2)+clnoise(2,l2)) &
                &               *((inlcls(2,l2)+clnoise(2,l2))*(inlcls(3,l1)) &
                &               *gTE12_num &
                &               +(inlcls(2,l1)+clnoise(2,l1))*(inlcls(3,l2)) &
                &               *gTE21_num)

             else

                fTB12=prefactor*wigner2*incls(3,l1) &
                &     *(all*(all+1.d0)+al1*(al1+1.d0)-al2*(al2+1.d0))
                fTB21=prefactor*wigner2*incls(3,l2) &
                &     *(all*(all+1.d0)+al2*(al2+1.d0)-al1*(al1+1.d0))
                fEB12=prefactor*wigner2*incls(2,l1) &
                &     *(all*(all+1.d0)+al1*(al1+1.d0)-al2*(al2+1.d0))
                fEB21=prefactor*wigner2*incls(2,l2) &
                &     *(all*(all+1.d0)+al2*(al2+1.d0)-al1*(al1+1.d0))

                aAlTB = aAlTB  &
                &     + fTB12**2/(inlcls(1,l1)+clnoise(1,l1)) &
                &               /(inlcls(4,l2)+clnoise(2,l2)) &
                &     + fTB21**2/(inlcls(1,l2)+clnoise(1,l2)) &
                &               /(inlcls(4,l1)+clnoise(2,l1))
                aAlEB = aAlEB  &
                &     + fEB12**2/(inlcls(2,l1)+clnoise(2,l1)) &
                &               /(inlcls(4,l2)+clnoise(2,l2)) &
                &     + fEB21**2/(inlcls(2,l2)+clnoise(2,l2)) &
                &               /(inlcls(4,l1)+clnoise(2,l1))

                aBlTBEB = aBlTBEB &
                &       + fTB12*fEB12*(inlcls(3,l1)) &
                &                    /(inlcls(4,l2)+clnoise(2,l2)) &
                &                    /(inlcls(1,l1)+clnoise(1,l1)) &
                &                    /(inlcls(2,l1)+clnoise(2,l1)) &
                &       + fTB21*fEB21*(inlcls(3,l2)) &
                &                    /(inlcls(4,l1)+clnoise(2,l1)) &
                &                    /(inlcls(1,l2)+clnoise(1,l2)) &
                &                    /(inlcls(2,l2)+clnoise(2,l2))

             end if
          end if

          ! case l2=l2max-1

          l2=l2max-1
          al2=dble(l2)

          Wl20_mid=0.
          Wl22_mid=Al1L2*2.d0* &
          &               dsqrt(all/al1/(al1+all-2.d0)/(al1+all+2.d0))

          if (dint(parity_test).ne.parity_test) then
             Wl22_mid=-Wl22_mid
          end if

          parity=(al1+all+al2)/2.d0

          if (l2.le.lmaxM) then

             wigner0=Wl20_mid
             wigner2=Wl22_mid

             prefactor=dsqrt((2.d0*all+1.d0)*(2.d0*al1+1.d0) &
             &              *(2.d0*al2+1.d0)/(16.d0*3.14159d0))

             if (dint(parity).eq.parity) then

                fTT12=prefactor*wigner0*(incls(1,l1) &
                &    *(all*(all+1.d0)+al1*(al1+1.d0)-al2*(al2+1.d0)) &
                &    +incls(1,l2) &
                &    *(all*(all+1.d0)+al2*(al2+1.d0)-al1*(al1+1.d0)))
                fEE12=prefactor*wigner2*(incls(2,l1) &
                &    *(all*(all+1.d0)+al1*(al1+1.d0)-al2*(al2+1.d0)) &
                &    +incls(2,l2) &
                &    *(all*(all+1.d0)+al2*(al2+1.d0)-al1*(al1+1.d0)))
                fTE12=prefactor*(wigner2*incls(3,l1) &
                &    *(all*(all+1.d0)+al1*(al1+1.d0)-al2*(al2+1.d0)) &
                &    +wigner0*incls(3,l2) &
                &    *(all*(all+1.d0)+al2*(al2+1.d0)-al1*(al1+1.d0)))
                fTE21=prefactor*(wigner2*incls(3,l2) &
                &    *(all*(all+1.d0)+al2*(al2+1.d0)-al1*(al1+1.d0)) &
                &    +wigner0*incls(3,l1) &
                &    *(all*(all+1.d0)+al1*(al1+1.d0)-al2*(al2+1.d0)))

                gTE12_den=(inlcls(1,l1)+clnoise(1,l1)) &
                &        *(inlcls(1,l2)+clnoise(1,l2)) &
                &        *(inlcls(2,l1)+clnoise(2,l1)) &
                &        *(inlcls(2,l2)+clnoise(2,l2)) &
                &        -((inlcls(3,l1)) &
                &        *(inlcls(3,l2)))**2
                gTE12_num=(inlcls(1,l2)+clnoise(1,l2)) &
                &        *(inlcls(2,l1)+clnoise(2,l1)) &
                &        *fTE12-(inlcls(3,l1)) &
                &        *(inlcls(3,l2))*fTE21
                gTE21_num=(inlcls(1,l1)+clnoise(1,l1)) &
                &        *(inlcls(2,l2)+clnoise(2,l2)) &
                &        *fTE21-(inlcls(3,l1)) &
                &        *(inlcls(3,l2))*fTE12

                aAlTT = aAlTT +fTT12**2 &
                &              /(inlcls(1,l1)+clnoise(1,l1)) &
                &              /(inlcls(1,l2)+clnoise(1,l2))
                aAlEE = aAlEE +fEE12**2 &
                &              /(inlcls(2,l1)+clnoise(2,l1)) &
                &              /(inlcls(2,l2)+clnoise(2,l2))
                aAlTE = aAlTE + (fTE12*gTE12_num+fTE21*gTE21_num) &
                &              /gTE12_den

                aBlTTEE= aBlTTEE+fTT12*fEE12 &
                &                *(inlcls(3,l1))*(inlcls(3,l2)) &
                &                /(inlcls(1,l1)+clnoise(1,l1)) &
                &                /(inlcls(1,l2)+clnoise(1,l2)) &
                &                /(inlcls(2,l1)+clnoise(2,l1)) &
                &                /(inlcls(2,l2)+clnoise(2,l2))
                aBlTTTE= aBlTTTE+fTT12/gTE12_den &
                &                /(inlcls(1,l1)+clnoise(1,l1)) &
                &                /(inlcls(1,l2)+clnoise(1,l2)) &
                &                *((inlcls(1,l1)+clnoise(1,l1)) &
                &                *(inlcls(3,l2))*gTE12_num &
                &               +(inlcls(1,l2)+clnoise(1,l2)) &
                &                 *(inlcls(3,l1))*gTE21_num)
                aBlEETE= aBlEETE+fEE12/gTE12_den &
                &                /(inlcls(2,l1)+clnoise(2,l1)) &
                &                /(inlcls(2,l2)+clnoise(2,l2)) &
                &                *((inlcls(2,l2)+clnoise(2,l2)) &
                &                *(inlcls(3,l1))*gTE12_num &
                &               +(inlcls(2,l1)+clnoise(2,l1)) &
                &                *(inlcls(3,l2))*gTE21_num)

             else

                fTB12=prefactor*wigner2*incls(3,l1) &
                &     *(all*(all+1.d0)+al1*(al1+1.d0)-al2*(al2+1.d0))
                fTB21=prefactor*wigner2*incls(3,l2) &
                &     *(all*(all+1.d0)+al2*(al2+1.d0)-al1*(al1+1.d0))
                fEB12=prefactor*wigner2*incls(2,l1) &
                &     *(all*(all+1.d0)+al1*(al1+1.d0)-al2*(al2+1.d0))
                fEB21=prefactor*wigner2*incls(2,l2) &
                &     *(all*(all+1.d0)+al2*(al2+1.d0)-al1*(al1+1.d0))

                aAlTB = aAlTB  &
                &     + fTB12**2/(inlcls(1,l1)+clnoise(1,l1)) &
                &               /(inlcls(4,l2)+clnoise(2,l2)) &
                &     + fTB21**2/(inlcls(1,l2)+clnoise(1,l2)) &
                &               /(inlcls(4,l1)+clnoise(2,l1))
                aAlEB = aAlEB  &
                &     + fEB12**2/(inlcls(2,l1)+clnoise(2,l1)) &
                &               /(inlcls(4,l2)+clnoise(2,l2)) &
                &     + fEB21**2/(inlcls(2,l2)+clnoise(2,l2)) &
                &               /(inlcls(4,l1)+clnoise(2,l1))

                aBlTBEB = aBlTBEB &
                &       + fTB12*fEB12*(inlcls(3,l1)) &
                &                    /(inlcls(4,l2)+clnoise(2,l2)) &
                &                    /(inlcls(1,l1)+clnoise(1,l1)) &
                &                    /(inlcls(2,l1)+clnoise(2,l1)) &
                &       + fTB21*fEB21*(inlcls(3,l2)) &
                &                    /(inlcls(4,l1)+clnoise(2,l1)) &
                &                    /(inlcls(1,l2)+clnoise(1,l2)) &
                &                    /(inlcls(2,l2)+clnoise(2,l2))

             end if
          end if

          ! sum over other cases, down to l2=l1+1

          do nn=1,l2max-1-l1-1 !! do 3

             l2=l2max-1-nn
             al2=dble(l2)

             ! recursion formula giving W_l2

             aa0=(al2+2.d0)*(al2+1.d0) &
             &             *dsqrt((-al2+all+al1)*(al2-all+al1+1.d0) &
             &             *(al2+all-al1+1.d0)*(al2+all+al1+2.d0))
             bb0=0.d0
             cc0=-(al2+1.d0)*(al2+2.d0) &
             &              *dsqrt((-al2+all+al1-1.d0)*(al2-all+al1+2.d0) &
             &              *(al2+all-al1+2.d0)*(al2+all+al1+3.d0))

             Wl20_inf=(bb0*Wl20_mid+cc0*Wl20_sup)/aa0

             aa2=(al2+2.d0)*dsqrt(((al2+1.d0)**2-4.d0) &
             &              *(-al2+all+al1)*(al2-all+al1+1.d0) &
             &              *(al2+all-al1+1.d0)*(al2+all+al1+2.d0))
             bb2=2.d0*(2.d0*al2+3.d0)*((al2+1.d0)*(al2+2.d0) &
             &              +all*(all+1.d0)-al1*(al1+1.d0))
             cc2=-(al2+1.d0)*dsqrt(((al2+2.d0)**2-4.d0) &
             &              *(-al2+all+al1-1.d0)*(al2-all+al1+2.d0) &
             &              *(al2+all-al1+2.d0)*(al2+all+al1+3.d0))

             Wl22_inf=(bb2*Wl22_mid+cc2*Wl22_sup)/aa2

             parity=(al1+all+al2)/2.d0

             if (l2.le.lmaxM) then

                wigner0=Wl20_inf
                wigner2=Wl22_inf

                prefactor=dsqrt((2.d0*all+1.d0)*(2.d0*al1+1.d0) &
                &                *(2.d0*al2+1.d0)/(16.d0*3.14159d0))

                if (dint(parity).eq.parity) then

                   fTT12=prefactor*wigner0*(incls(1,l1) &
                   &    *(all*(all+1.d0)+al1*(al1+1.d0)-al2*(al2+1.d0)) &
                   &    +incls(1,l2) &
                   &    *(all*(all+1.d0)+al2*(al2+1.d0)-al1*(al1+1.d0)))
                   fEE12=prefactor*wigner2*(incls(2,l1) &
                   &    *(all*(all+1.d0)+al1*(al1+1.d0)-al2*(al2+1.d0)) &
                   &    +incls(2,l2) &
                   &    *(all*(all+1.d0)+al2*(al2+1.d0)-al1*(al1+1.d0)))
                   fTE12=prefactor*(wigner2*incls(3,l1) &
                   &    *(all*(all+1.d0)+al1*(al1+1.d0)-al2*(al2+1.d0)) &
                   &    +wigner0*incls(3,l2) &
                   &    *(all*(all+1.d0)+al2*(al2+1.d0)-al1*(al1+1.d0)))
                   fTE21=prefactor*(wigner2*incls(3,l2) &
                   &    *(all*(all+1.d0)+al2*(al2+1.d0)-al1*(al1+1.d0)) &
                   &    +wigner0*incls(3,l1) &
                   &    *(all*(all+1.d0)+al1*(al1+1.d0)-al2*(al2+1.d0)))

                   gTE12_den=(inlcls(1,l1)+clnoise(1,l1)) &
                   &        *(inlcls(1,l2)+clnoise(1,l2)) &
                   &        *(inlcls(2,l1)+clnoise(2,l1)) &
                   &        *(inlcls(2,l2)+clnoise(2,l2)) &
                   &        -((inlcls(3,l1)) &
                   &        *(inlcls(3,l2)))**2
                   gTE12_num=(inlcls(1,l2)+clnoise(1,l2)) &
                   &         *(inlcls(2,l1)+clnoise(2,l1))*fTE12 &
                   &        -(inlcls(3,l1))*(inlcls(3,l2))*fTE21
                   gTE21_num=(inlcls(1,l1)+clnoise(1,l1)) &
                   &         *(inlcls(2,l2)+clnoise(2,l2))*fTE21 &
                   &         -(inlcls(3,l1))*(inlcls(3,l2))*fTE12

                   aAlTT = aAlTT +fTT12**2 &
                   &             /(inlcls(1,l1)+clnoise(1,l1)) &
                   &             /(inlcls(1,l2)+clnoise(1,l2))
                   aAlEE = aAlEE +fEE12**2 &
                   &             /(inlcls(2,l1)+clnoise(2,l1)) &
                   &             /(inlcls(2,l2)+clnoise(2,l2))
                   aAlTE = aAlTE + (fTE12*gTE12_num+fTE21*gTE21_num) &
                   &              /gTE12_den

                   aBlTTEE= aBlTTEE+fTT12*fEE12 &
                   &                *(inlcls(3,l1))*(inlcls(3,l2)) &
                   &                /(inlcls(1,l1)+clnoise(1,l1)) &
                   &                /(inlcls(1,l2)+clnoise(1,l2)) &
                   &                /(inlcls(2,l1)+clnoise(2,l1)) &
                   &                /(inlcls(2,l2)+clnoise(2,l2))
                   aBlTTTE= aBlTTTE+fTT12/gTE12_den &
                   &                /(inlcls(1,l1)+clnoise(1,l1)) &
                   &                /(inlcls(1,l2)+clnoise(1,l2)) &
                   &                *((inlcls(1,l1)+clnoise(1,l1)) &
                   &                *(inlcls(3,l2))*gTE12_num &
                   &               +(inlcls(1,l2)+clnoise(1,l2)) &
                   &                *(inlcls(3,l1))*gTE21_num)
                   aBlEETE= aBlEETE+fEE12/gTE12_den &
                   &                /(inlcls(2,l1)+clnoise(2,l1)) &
                   &                /(inlcls(2,l2)+clnoise(2,l2)) &
                   &                *((inlcls(2,l2)+clnoise(2,l2)) &
                   &                *(inlcls(3,l1))*gTE12_num &
                   &               +(inlcls(2,l1)+clnoise(2,l1)) &
                   &                 *(inlcls(3,l2))*gTE21_num)

                else

                   fTB12=prefactor*wigner2*incls(3,l1) &
                   &    *(all*(all+1.d0)+al1*(al1+1.d0)-al2*(al2+1.d0))
                   fTB21=prefactor*wigner2*incls(3,l2) &
                   &    *(all*(all+1.d0)+al2*(al2+1.d0)-al1*(al1+1.d0))
                   fEB12=prefactor*wigner2*incls(2,l1) &
                   &    *(all*(all+1.d0)+al1*(al1+1.d0)-al2*(al2+1.d0))
                   fEB21=prefactor*wigner2*incls(2,l2) &
                   &    *(all*(all+1.d0)+al2*(al2+1.d0)-al1*(al1+1.d0))

                   aAlTB = aAlTB  &
                   &     + fTB12**2/(inlcls(1,l1)+clnoise(1,l1)) &
                   &               /(inlcls(4,l2)+clnoise(2,l2)) &
                   &     + fTB21**2/(inlcls(1,l2)+clnoise(1,l2)) &
                   &               /(inlcls(4,l1)+clnoise(2,l1))
                   aAlEB = aAlEB  &
                   &     + fEB12**2/(inlcls(2,l1)+clnoise(2,l1)) &
                   &               /(inlcls(4,l2)+clnoise(2,l2)) &
                   &     + fEB21**2/(inlcls(2,l2)+clnoise(2,l2)) &
                   &               /(inlcls(4,l1)+clnoise(2,l1))

                   aBlTBEB = aBlTBEB &
                   &       + fTB12*fEB12*(inlcls(3,l1)) &
                   &               /(inlcls(4,l2)+clnoise(2,l2)) &
                   &               /(inlcls(1,l1)+clnoise(1,l1)) &
                   &               /(inlcls(2,l1)+clnoise(2,l1)) &
                   &       + fTB21*fEB21*inlcls(3,l2) &
                   &               /(inlcls(4,l1)+clnoise(2,l1)) &
                   &               /(inlcls(1,l2)+clnoise(1,l2)) &
                   &               /(inlcls(2,l2)+clnoise(2,l2))

                end if
             end if

             ! shift for the next recursion

             Wl20_sup=Wl20_mid
             Wl20_mid=Wl20_inf

             Wl22_sup=Wl22_mid
             Wl22_mid=Wl22_inf

          end do !! do 3

          ! last case, for l2=l1

          l2=l1
          al2=al1

          ! recursion formula giving W_l2

          aa0=(al2+2.d0)*(al2+1.d0)* &
          &              dsqrt((-al2+all+al1)*(al2-all+al1+1.d0) &
          &              *(al2+all-al1+1.d0)*(al2+all+al1+2.d0))
          bb0=0.d0
          cc0=-(al2+1.d0)*(al2+2.d0)* &
          &              dsqrt((-al2+all+al1-1.d0)*(al2-all+al1+2.d0) &
          &              *(al2+all-al1+2.d0)*(al2+all+al1+3.d0))

          Wl20_inf=(bb0*Wl20_mid+cc0*Wl20_sup)/aa0

          aa2=(al2+2.d0)*dsqrt(((al2+1.d0)**2-4.d0) &
          &              *(-al2+all+al1)*(al2-all+al1+1.d0) &
          &              *(al2+all-al1+1.d0)*(al2+all+al1+2.d0))
          bb2=2.d0*(2.d0*al2+3.d0)*((al2+1.d0)*(al2+2.d0) &
          &              +all*(all+1.d0)-al1*(al1+1.d0))
          cc2=-(al2+1.d0)*dsqrt(((al2+2.d0)**2-4.d0) &
          &              *(-al2+all+al1-1.d0)*(al2-all+al1+2.d0) &
          &              *(al2+all-al1+2.d0)*(al2+all+al1+3.d0))

          Wl22_inf=(bb2*Wl22_mid+cc2*Wl22_sup)/aa2

          parity=(al1+all+al2)/2.d0

          if (l2.le.lmaxM) then

             wigner0=Wl20_inf
             wigner2=Wl22_inf

             prefactor=dsqrt((2.d0*all+1.d0)*(2.d0*al1+1.d0) &
             &                *(2.d0*al2+1.d0)/(16.d0*pi))

             if (dint(parity).eq.parity) then

                fTT12=prefactor*wigner0*2.d0*incls(1,l1)*all*(all+1.d0)
                fEE12=prefactor*wigner2*2.d0*incls(2,l1)*all*(all+1.d0)
                fTE12=prefactor*(wigner2+wigner0)*incls(3,l1) &
                &                *all*(all+1.d0)

                gTE12_den=(inlcls(1,l1)+clnoise(1,l1))**2 &
                &        *(inlcls(2,l1)+clnoise(2,l1))**2 &
                &        -(inlcls(3,l1))**4
                gTE12_num=((inlcls(1,l1)+clnoise(1,l1)) &
                &         *(inlcls(2,l1)+clnoise(2,l1)) &
                &         -(inlcls(3,l1))**2)*fTE21

                aAlTT = aAlTT + fTT12**2 &
                &               /(2.d0*(inlcls(1,l1)+clnoise(1,l1))**2)
                aAlEE = aAlEE + fEE12**2 &
                &               /(2.d0*(inlcls(2,l1)+clnoise(2,l1))**2)
                aAlTE = aAlTE + fTE12*gTE12_num &
                &               /gTE12_den

                aBlTTEE= aBlTTEE+fTT12*fEE12 &
                &               /2.d0*(inlcls(3,l1))**2 &
                &               /(inlcls(1,l1)+clnoise(1,l1))**2 &
                &               /(inlcls(2,l1)+clnoise(2,l1))**2
                aBlTTTE= aBlTTTE+fTT12/gTE12_den &
                &              /(inlcls(1,l1)+clnoise(1,l1)) &
                &              *(inlcls(3,l1))*gTE12_num
                aBlEETE= aBlEETE+fEE12/gTE12_den &
                &              /(inlcls(2,l1)+clnoise(2,l1)) &
                &              *(inlcls(3,l1))*gTE12_num

             else

                fTB12=prefactor*wigner2*incls(3,l1)*all*(all+1.d0)
                fEB12=prefactor*wigner2*incls(2,l1)*all*(all+1.d0)

                aAlTB = aAlTB  &
                &     + fTB12**2/(inlcls(1,l1)+clnoise(1,l1)) &
                &               /(inlcls(4,l1)+clnoise(2,l1))
                aAlEB = aAlEB  &
                &     + fEB12**2/(inlcls(2,l1)+clnoise(2,l1)) &
                &               /(inlcls(4,l1)+clnoise(2,l1))

                aBlTBEB = aBlTBEB &
                &       + fTB12*fEB12*(inlcls(3,l1)) &
                &                    /(inlcls(4,l1)+clnoise(2,l1)) &
                &                    /(inlcls(1,l1)+clnoise(1,l1)) &
                &                    /(inlcls(2,l1)+clnoise(2,l1))

             end if
          end if

       end do !! do 2

       ! finally computing the noise

       aAlTT = all*(all+1.d0)*(2.d0*all+1.d0)/aAlTT
       aAlEE = all*(all+1.d0)*(2.d0*all+1.d0)/aAlEE
       aAlTE = all*(all+1.d0)*(2.d0*all+1.d0)/aAlTE
       aAlTB = all*(all+1.d0)*(2.d0*all+1.d0)/aAlTB
       aAlEB = all*(all+1.d0)*(2.d0*all+1.d0)/aAlEB
       aBlTTEE= all*(all+1.d0)*(2.d0*all+1.d0)/aBlTTEE
       aBlTTTE= all*(all+1.d0)*(2.d0*all+1.d0)/aBlTTTE
       aBlEETE= all*(all+1.d0)*(2.d0*all+1.d0)/aBlEETE
       aBlTBEB= all*(all+1.d0)*(2.d0*all+1.d0)/aBlTBEB

       LensingNoise(1,1) = aAlTT
       LensingNoise(2,2) = aAlEE
       LensingNoise(3,3) = aAlTE
       LensingNoise(4,4) = aAlTB
       LensingNoise(5,5) = aAlEB

       LensingNoise(1,2) = aAlTT*aAlEE/aBlTTEE
       LensingNoise(1,3) = aAlTT*aAlTE/aBlTTTE
       LensingNoise(1,4) = 0.
       LensingNoise(1,5) = 0.
       LensingNoise(2,3) = aAlEE*aAlTE/aBlEETE
       LensingNoise(2,4) = 0.
       LensingNoise(2,5) = 0.
       LensingNoise(3,4) = 0.
       LensingNoise(3,5) = 0.
       LensingNoise(4,5) = aAlTB*aAlEB/aBlTBEB

       LensingNoise(2,1) = LensingNoise(1,2)
       LensingNoise(3,1) = LensingNoise(1,3)
       LensingNoise(4,1) = LensingNoise(1,4)
       LensingNoise(5,1) = LensingNoise(1,5)
       LensingNoise(3,2) = LensingNoise(2,3)
       LensingNoise(4,2) = LensingNoise(2,4)
       LensingNoise(5,2) = LensingNoise(2,5)
       LensingNoise(4,3) = LensingNoise(3,4)
       LensingNoise(5,3) = LensingNoise(3,5)
       LensingNoise(5,4) = LensingNoise(4,5)

       ! compute minimum variance estimator based on TT,EE,TE,TB,EB

       do l1=1,5
          GaussVecLens(l1,1)=1.
          do l2=1,5
             InvLensingNoise(l1,l2)=LensingNoise(l1,l2)
          end do
       end do

       call gaussj(InvLensingNoise,5,5,GaussVecLens,1,1)

       MinVarNoise_L(jj)=0.

       do l1=1,5
          do l2=1,5
             MinVarNoise_L(jj)=MinVarNoise_L(jj)+InvLensingNoise(l1,l2)
          end do
       end do

       MinVarNoise_L(jj)=1.d0/MinVarNoise_L(jj)

       ! compute minimum variance estimator based only on EE,EB
       !MinVarNoisePol_L(jj)=1.d0/
       !&       (1.d0/LensingNoise(2,2)+1.d0/LensingNoise(5,5))
       ! compute minimum variance estimator based only on TT
       !MinVarNoiseTemp_L(jj)=LensingNoise(1,1)

       ! print the results

!       write(37,'(1I5,8E13.4)')ll, &
!       &        LensingNoise(1,1)*all*(all+1.d0)/2.d0/pi, &
!       &        LensingNoise(2,2)*all*(all+1.d0)/2.d0/pi, &
!       &        LensingNoise(3,3)*all*(all+1.d0)/2.d0/pi, &
!       &        LensingNoise(4,4)*all*(all+1.d0)/2.d0/pi, &
!       &        LensingNoise(5,5)*all*(all+1.d0)/2.d0/pi, &
!       &        MinVarNoise_L(jj)*all*(all+1.d0)/2.d0/pi, &
!       &        incls(4,int(all))*all*(all+1.d0)/2.d0/pi

       alogllens(jj)=log(all)
    end do

    !close(37)
    !     interpolate to any value of l
    !     (interpolation performed in log(l) space)
    allocate(ddMinVarNoise_L(1:jjmaxT))
    call spline(alogllens,MinVarNoise_L,jjmaxT, &
    &     0.d0,0.d0,ddMinVarNoise_L)
    !call spline(alogllens,MinVarNoisePol_L,jjmaxT,
    !&     0.d0,0.d0,ddMinVarNoisePol_L)
    !call spline(alogllens,MinVarNoiseTemp_L,jjmaxT,
    !&     0.d0,0.d0,ddMinVarNoiseTemp_L)
!    open(unit=34,file="test_Nldd.dat",form="formatted", &
!         &                       status='replace', iostat=err)
!    if (err /= 0) print *, "Probleme d'ecriture du fichier test..."
    do l=2,lmaxM
       call splint(alogllens,MinVarNoise_L,ddMinVarNoise_L, &
       &        jjmaxT,log(dble(l)),nldd(l))
       !call splint(alogllens,MinVarNoisePol_L,ddMinVarNoisePol_L,
       !&        jjmaxT,log(dble(l)),omBPol_D(l))
       !call splint(alogllens,MinVarNoiseTemp_L,ddMinVarNoiseTemp_L,
       !&        jjmaxT,log(dble(l)),omBTemp_D(l))
       !write(34,'(1I5,1E13.4)') l, nldd(l)*l*(l+1.d0)/2.d0/pi
    end do
    !close(34)

    !print*,code,"FIN"
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

  end subroutine calc_Nldd


 !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
 SUBROUTINE gaussj(a,n,np,b,m,mp)
   !c      implicit double precision (a-h,o-z)
   INTEGER m,mp,n,np,NMAX
   REAL*8 a(np,np),b(np,mp)
   PARAMETER (NMAX=50)
   !c      Linear equation solution by Gauss-Jordan elimination,
   !c      equation (2.1.1) above. a(1:n,1:n) is an input matrix stored in an
   !c      array of physical dimensions np by np.
   !c      b(1:n,1:m) is an in-put matrix containing the m right-hand side
   !c      vectors, stored in an array of physical dimensions np by mp.
   !c      On output, a(1:n,1:n) is replaced by its matrix inverse,
   !c      and b(1:n,1:m) is replaced by the corresponding set of solution
   !c      vectors. Parameter: NMAX is the largest anticipated value of n.
   INTEGER i,icol,irow,j,k,l,ll,indxc(NMAX),indxr(NMAX),ipiv(NMAX)
   !c      The integer arrays ipiv, indxr,andindxc are used for bookkeeping
   !c      on the pivoting.
   REAL*8 big,dum,pivinv

   do 11 j=1,n
      ipiv(j)=0
11 enddo
   do 22 i=1,n
      !c      This is the main loop over the columns to be re-duced.
      big=0.
      do 13 j=1,n
         !c      This is the outer loop of the search for a pivot ele-ment.
            if(ipiv(j).ne.1)then
               do 12 k=1,n
                  if (ipiv(k).eq.0) then
                     if (abs(a(j,k)).ge.big)then
                        big=abs(a(j,k))
                        irow=j
                        icol=k
                     endif
                  endif
12             enddo
            endif
13       enddo
         ipiv(icol)=ipiv(icol)+1
         if (irow.ne.icol) then
            do 14 l=1,n
               dum=a(irow,l)
               a(irow,l)=a(icol,l)
               a(icol,l)=dum
14          enddo
            do 15 l=1,m
               dum=b(irow,l)
               b(irow,l)=b(icol,l)
               b(icol,l)=dum
15          enddo
         endif
         indxr(i)=irow
         !c We are now ready to divide the pivot row by the pivot element,
         !c located at irow and icol.
         indxc(i)=icol
         if (a(icol,icol).eq.0.) stop 'singular matrix in gaussj'
         !c    singular matrix in gaussj
         pivinv=1./a(icol,icol)
         a(icol,icol)=1.
         do 16 l=1,n
            a(icol,l)=a(icol,l)*pivinv
16       enddo
         do 17 l=1,m
            b(icol,l)=b(icol,l)*pivinv
17       enddo
         do 21 ll=1,n
            !c   Next, we reduce the rows...
            if(ll.ne.icol)then
               !c   ...except for the pivot one, of course.
               dum=a(ll,icol)
               a(ll,icol)=0.
               do 18 l=1,n
                  a(ll,l)=a(ll,l)-a(icol,l)*dum
18             enddo
               do 19 l=1,m
                  b(ll,l)=b(ll,l)-b(icol,l)*dum
19             enddo
            endif
21       enddo
22    enddo
      !c      This is the end of the main loop over columns of the reduction.
      do 24 l=n,1,-1
         !c    It only remains to unscramble the solution in view of the
         !c    column interchanges. We do this by in-terchanging pairs of
         !c    columns in the reverse order that the permutation was built up.
         if(indxr(l).ne.indxc(l))then
            do 23 k=1,n
               dum=a(k,indxr(l))
               a(k,indxr(l))=a(k,indxc(l))
               a(k,indxc(l))=dum
23          enddo
         endif
24    enddo
      return
      !c  And we are done.
    END SUBROUTINE gaussj
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

    SUBROUTINE spline(x,y,n,yp1,ypn,y2)
      INTEGER n,NMAX
      double precision yp1,ypn,x(n),y(n),y2(n)
      PARAMETER (NMAX=500)
      INTEGER i,k
      double precision p,qn,sig,un,u(NMAX)
      if (yp1.gt..99d30) then
         y2(1)=0.d0
         u(1)=0.d0
      else
         y2(1)=-0.5d0
         u(1)=(3.d0/(x(2)-x(1)))*((y(2)-y(1))/(x(2)-x(1))-yp1)
      endif
      do 11 i=2,n-1
         sig=(x(i)-x(i-1))/(x(i+1)-x(i-1))
         p=sig*y2(i-1)+2.d0
         y2(i)=(sig-1.d0)/p
         u(i)=(6.d0*((y(i+1)-y(i))/(x(i+1)-x(i))-(y(i)-y(i-1)) &
         &       /(x(i)-x(i-1)))/(x(i+1)-x(i-1))-sig*u(i-1))/p
11       continue
         if (ypn.gt..99d30) then
            qn=0.d0
            un=0.d0
         else
            qn=0.5d0
            un=(3.d0/(x(n)-x(n-1)))*(ypn-(y(n)-y(n-1))/(x(n)-x(n-1)))
         endif
         y2(n)=(un-qn*u(n-1))/(qn*y2(n-1)+1.d0)
         do 12 k=n-1,1,-1
            y2(k)=y2(k)*y2(k+1)+u(k)
12          continue
            return
            !C  (C) Copr. 1986-92 Numerical Recipes Software =$j*m,).
     END subroutine spline

     !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

     subroutine splint(xa,ya,y2a,n,x,y)
       !c
       integer n
       double precision x,y,xa(n),y2a(n),ya(n)
       !c
       integer k,khi,klo
       double precision a,b,h
       !c
       klo=1
       khi=n
1      if (khi-klo.gt.1) then
          k=(khi+klo)/2
          if(xa(k).gt.x)then
             khi=k
          else
             klo=k
          end if
          goto 1
       end if
       h=xa(khi)-xa(klo)
       if (h.eq.0.) stop 'bad xa input in splint'
       a=(xa(khi)-x)/h
       b=(x-xa(klo))/h
       y=a*ya(klo)+b*ya(khi) + &
       &       ((a**3-a)*y2a(klo)+(b**3-b)*y2a(khi))*(h**2)/6.
       return
     end subroutine splint
     !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

end module lensnoise

!----------------------------------------------------------------------------------------
