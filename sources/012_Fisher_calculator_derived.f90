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

!> @file 012_Fisher_calculator_derived.f90
!! This file contains the code that computes the matrix that is needed to obtain the
!! Fisher matrix for derived parameters.

!----------------------------------------------------------------------------------------
!> This module contains the code that computes the matrix that can be used to compute the
!! the Fisher matrix for derived parameters.

!> @author Marco Raveri

module Fisher_calculator_Derived

    use constants, only: k_B, eV, mass_ratio_He_H
    use precision
    use Transfer, only: Transfer_SortAndIndexRedshifts, Transfer_output_sig8, MT
    use ModelParams
    use cosmicfish_types
    use omp_lib
    use cosmicfish_utilities
    use cosmicfish_camb_interface
    use CAMB

#ifdef COSMICFISH_EFTCAMB
    use EFTinitialization
#endif

    implicit none

    private

    public  dimension_derived_parameters, Compute_Derived_Parameters, Fisher_derived_param_names, &
        & Derived_Parameter_Derivative, Fisher_derived, save_derived_Fisher_to_file, save_derived_paramnames_to_file

contains

    ! ---------------------------------------------------------------------------------------------
    !> This subroutine checks the compatibility of input parameters. It will take the FP parameters
    !> and adapt the P parameters so that you get what you want.
    subroutine derived_fisher_parameter_compatibility( P, FP )

        implicit none

        Type(CAMBparams)        :: P         !< Input and output CAMBparams
        Type(cosmicfish_params)  :: FP        !< Input and output Cosmicfish params

        integer :: i

        integer :: j
        character(LEN=20) :: str, str2

        ! if we want tomography for the derived parameters fix the transfer functions parameters accorgingly:
        if ( FP%fisher_der%want_sigma8 ) then
            ! write the relevant quantities in Transfer:
            P%PK_WantTransfer = .True.
            P%WantTransfer    = .True.
            P%Transfer%high_precision = .True.
            P%Transfer%PK_num_redshifts = FP%fisher_der%FD_num_redshift
            do i=1, FP%fisher_der%FD_num_redshift
                P%transfer%PK_redshifts(i) = FP%fisher_der%FD_redshift(i)
            end do
            ! ensure ordering in redshift:
            call Transfer_SortAndIndexRedshifts(P%Transfer)
            ! copy back the redshift:
            do i=1, FP%fisher_der%FD_num_redshift
                FP%fisher_der%FD_redshift(i) = P%transfer%PK_redshifts(i)
            end do
        end if

        ! we do not want other things to be computed:
        num_redshiftwindows = 0

    end subroutine derived_fisher_parameter_compatibility

    ! ---------------------------------------------------------------------------------------------
    !> This subroutine returns the number of derived parameters
    subroutine dimension_derived_parameters( P, FP, num )

        implicit none

        Type(CAMBparams)      , intent(in)  :: P         !< Input CAMBparams
        Type(cosmicfish_params), intent(in)  :: FP        !< Input Cosmicfish params
        integer               , intent(out) :: num       !< Output size of derived parameters

        ! start with the base derived parameters:
        num = 0
        if ( FP%fisher_der%want_omegab ) then
            num = num + 1
        end if
        if ( FP%fisher_der%want_omegac ) then
            num = num + 1
        end if
        if ( FP%fisher_der%want_omegan ) then
            num = num + 1
        end if
        if ( FP%fisher_der%want_omegav ) then
            num = num + 1
        end if
        if ( FP%fisher_der%want_omegak ) then
            num = num + 1
        end if
        if ( FP%fisher_der%want_omegam ) then
            num = num + 1
        end if
        if ( FP%fisher_der%want_theta  ) then
            num = num + 1
        end if
        if ( FP%fisher_der%want_mnu    ) then
            num = num + 1
        end if
        if ( FP%fisher_der%want_zre    ) then
            num = num + 1
        end if
        if ( FP%fisher_der%want_neff   ) then
            num = num + 1
        end if

#ifdef COSMICFISH_MGCAMB
        if ( FP%fisher_der%want_mu0 ) then
            num = num + 1
        end if

        if ( FP%fisher_der%want_gam0 ) then
            num = num + 1
        end if
#endif

        ! add sigma8 if wanted:
        if ( FP%fisher_der%want_sigma8 ) then
            num = num + FP%fisher_der%FD_num_redshift
        end if
        ! add Hubble if wanted:
        if ( FP%fisher_der%want_loghubble ) then
            num = num + FP%fisher_der%FD_num_redshift
        end if
        ! add DA if wanted:
        if ( FP%fisher_der%want_logDA ) then
            num = num + FP%fisher_der%FD_num_redshift
        end if


    end subroutine dimension_derived_parameters

    ! ---------------------------------------------------------------------------------------------
    !> This subroutine returns an array of derived parameters given CAMB parameters
    subroutine Compute_Derived_Parameters( P, FP, derived, num, err )

        implicit none

        Type(CAMBparams)      , intent(in)  :: P            !< Input CAMBparams
        Type(cosmicfish_params), intent(in)  :: FP           !< Input Cosmicfish params
        real(dl), dimension(:), intent(out) :: derived      !< Output array with the supernovae magnitude
        integer               , intent(in)  :: num          !< number of derived parameters.
        integer, intent(out)    :: err                      !< Output error code:
                                                            !< 0 = all fine
                                                            !< 1 = error in input
                                                            !< 2 = error in computation
                                                            !< 3 = error in quality
                                                            !< 4 = routine specific

        ! other definitions:
        integer  :: i, j, nu_i, ind, error
        logical  :: EFTsuccess
        real(dl) :: conv, z
        type (CAMBdata)  :: OutData
        ! initialization:
        derived = 0._dl
        err     = 0

        ! set the camb parameters:
        call CAMBParams_Set(P)

#ifdef COSMICFISH_EFTCAMB
        ! check stability if EFTCAMB is used:
        if (P%EFTflag /= 0) then
            call EFTCAMB_initialization(EFTsuccess)
            if (.not. EFTsuccess) then
                write(*,*) 'EFTCAMB: unstable model'
                err = 4
                return
            end if
        end if
#endif

        ! get the transfer functions for sigma_8
        if ( FP%fisher_der%want_sigma8 ) then
            error = 0
            call CAMB_GetTransfers(P, OutData, error)
            ! test for errors in the computation
            if (error/=0) then
                write (*,*) 'Error result '//trim(global_error_message)
                err = 2
                return
            endif
        end if

        ! write the parameter array:
        ind = 0
        if ( FP%fisher_der%want_omegab ) then
            ind = ind + 1
            derived(ind) = CP%omegab ! Baryons density
        end if
        if ( FP%fisher_der%want_omegac ) then
            ind = ind + 1
            derived(ind) = CP%omegac ! CDM density
        end if
        if ( FP%fisher_der%want_omegan ) then
            ind = ind + 1
            derived(ind) = CP%omegan ! Neutrino density
        end if
        if ( FP%fisher_der%want_omegav ) then
            ind = ind + 1
            derived(ind) = CP%omegav ! DE density
        end if
        if ( FP%fisher_der%want_omegak ) then
            ind = ind + 1
            derived(ind) = CP%omegak ! Curvature density
        end if
        if ( FP%fisher_der%want_omegam ) then
            ind = ind + 1
            derived(ind) = 1-CP%omegak-CP%omegav ! Matter density
        end if
        if ( FP%fisher_der%want_theta  ) then
            ind = ind + 1
            derived(ind) = 100*CosmomcTheta()    ! Theta CosmoMC
        end if
        if ( FP%fisher_der%want_mnu    ) then
            ind = ind + 1
            if (CP%Num_Nu_Massive == 1) then
                do nu_i=1, CP%Nu_mass_eigenstates
                    conv = k_B*(8*grhor/grhog/7)**0.25*CP%tcmb/eV * &
                        (CP%nu_mass_degeneracies(nu_i)/CP%nu_mass_numbers(nu_i))**0.25 !approx 1.68e-4
                    derived(ind) = conv*nu_masses(nu_i)
                end do
            else
                derived(ind) = 0._dl ! m_nu in eV
            end if
        end if
        if ( FP%fisher_der%want_zre    ) then
            ind = ind + 1
            derived(ind)  = CP%Reion%redshift ! reionization redshift
        end if
        if ( FP%fisher_der%want_neff   ) then
            ind = ind + 1
            derived(ind) = CP%Num_Nu_massless + CP%Num_Nu_massive
        end if

#ifdef COSMICFISH_MGCAMB
        if ( FP%fisher_der%want_mu0 ) then
            ind = ind + 1
            if (CP%MGC_model ==11) derived(ind) = 1 + CP%E11_mg * CP%omegav
            if (CP%MGC_model ==12) derived(ind) = 1 + CP%E11_mg
        end if

        if ( FP%fisher_der%want_gam0 ) then
            ind = ind + 1
            if (CP%MGC_model ==11) derived(ind) = 1 + CP%E11_mg * CP%omegav
            if (CP%MGC_model ==12) derived(ind) = 1 + CP%E11_mg
        end if
#endif

        ! add sigma8 if wanted:
        if ( FP%fisher_der%want_sigma8 ) then
            do i=1, FP%fisher_der%FD_num_redshift
                ind = ind + 1
                j = P%Transfer%PK_redshifts_index(i)
                derived(ind) = OutData%MTrans%sigma_8(i,1)
            end do
        end if
        ! add Hubble if wanted:
        if ( FP%fisher_der%want_loghubble ) then
            do i=1, FP%fisher_der%FD_num_redshift
                ind = ind + 1
                z = FP%fisher_der%FD_redshift(i)
                derived(ind) = log10(Hofz(z))
            end do
        end if
        ! add DA if wanted:
        if ( FP%fisher_der%want_logDA ) then
            do i=1, FP%fisher_der%FD_num_redshift
                ind = ind + 1
                z = FP%fisher_der%FD_redshift(i)
                if (AngularDiameterDistance(z)==0._dl) then
                    derived(ind) = -16._dl ! protection against redshift zero...
                else
                    derived(ind) = log10(AngularDiameterDistance(z))
                end if
            end do
        end if


        if ( num/=ind ) then
            write(*,*) 'ERROR: Compute_Derived_Parameters num is not: ', num
            write(*,*) 'If you added a parameter change also dimension_derived_parameters'
            stop
        end if

    end subroutine Compute_Derived_Parameters

    ! ---------------------------------------------------------------------------------------------
    !> This subroutine returns the name corresponding to a derived parameter number.
    subroutine Fisher_derived_param_names(P, FP, param_number, param_name, param_name_latex)

        implicit none

        Type(CAMBparams)      , intent(in)  :: P            !< Input CAMBparams
        Type(cosmicfish_params), intent(in)  :: FP           !< Input Cosmicfish params
        integer               , intent(in)  :: param_number !< Input number of parameter for which the name is wanted
        character(*)          , intent(out) :: param_name   !< Output name of the parameter
        character(*)          , intent(out), optional :: param_name_latex   !< Output name of the parameter in LaTeX format

        integer  :: i, j, ind, num_param
        character(len=20) str
        real(dl) :: z

        num_param = 0

        ! standard parameters:

        if ( FP%fisher_der%want_omegab          )  num_param  = num_param + 1 ! omegab
        if ( param_number == num_param          )  then; param_name = 'omegab'; if ( present(param_name_latex) ) param_name_latex = '\Omega_b'; return ;
        end if
        if ( FP%fisher_der%want_omegac          )  num_param  = num_param + 1 ! omegac
        if ( param_number == num_param          )  then; param_name = 'omegac'; if ( present(param_name_latex) ) param_name_latex = '\Omega_c'; return ;
        end if
        if ( FP%fisher_der%want_omegan          )  num_param  = num_param + 1 ! omegan
        if ( param_number == num_param          )  then; param_name = 'omeganu'; if ( present(param_name_latex) ) param_name_latex = '\Omega_{\nu}'; return ;
        end if
        if ( FP%fisher_der%want_omegav          )  num_param  = num_param + 1 ! omegav
        if ( param_number == num_param          )  then; param_name = 'omegav'; if ( present(param_name_latex) ) param_name_latex = '\Omega_{\Lambda}'; return ;
        end if
        if ( FP%fisher_der%want_omegak          )  num_param  = num_param + 1 ! omegak
        if ( param_number == num_param          )  then; param_name = 'omegak'; if ( present(param_name_latex) ) param_name_latex = '\Omega_{K}'; return ;
        end if
        if ( FP%fisher_der%want_omegam          )  num_param  = num_param + 1 ! omegam
        if ( param_number == num_param          )  then; param_name = 'omegam'; if ( present(param_name_latex) ) param_name_latex = '\Omega_{m}'; return ;
        end if
        if ( FP%fisher_der%want_theta           )  num_param  = num_param + 1 ! theta
        if ( param_number == num_param          )  then; param_name = 'theta'; if ( present(param_name_latex) ) param_name_latex = '100\theta_{MC}'; return ;
        end if
        if ( FP%fisher_der%want_mnu             )  num_param  = num_param + 1 ! mnu
        if ( param_number == num_param          )  then; param_name = 'mnu'; if ( present(param_name_latex) ) param_name_latex = '\Sigma m_\nu'; return ;
        end if
        if ( FP%fisher_der%want_zre             )  num_param  = num_param + 1 ! zre
        if ( param_number == num_param          )  then; param_name = 'z_re'; if ( present(param_name_latex) ) param_name_latex = 'z_{\rm re}'; return ;
        end if
        if ( FP%fisher_der%want_neff            )  num_param  = num_param + 1 ! neff
        if ( param_number == num_param          )  then; param_name = 'Neff'; if ( present(param_name_latex) ) param_name_latex = 'N_{\rm eff}'; return ;
        end if


#ifdef COSMICFISH_MGCAMB
        if ( FP%fisher_der%want_mu0            )  num_param  = num_param + 1 !mu(z=0)
        if ( param_number == num_param          )  then; param_name = 'mu0'; if( present(param_name_latex) ) param_name_latex='\mu_0'; return ;
        end if
        if ( FP%fisher_der%want_gam0            )  num_param  = num_param + 1 !gamma(z=0)
        if ( param_number == num_param          )  then; param_name = 'gamma0';if ( present(param_name_latex) ) param_name_latex='\gamma_0'; return ;
        end if
#endif

        ind = num_param
        ! add sigma8 if wanted:
        if ( FP%fisher_der%want_sigma8 ) then
            do i=1, FP%fisher_der%FD_num_redshift
                ind = ind + 1
                if ( param_number == ind )  then
                    j = P%Transfer%PK_redshifts_index(i)
                    write (str, *) i
                    str = adjustl(str)
                    param_name = 'sigma8_'//trim(str)
                    write (str,'(F20.2)') P%Transfer%redshifts(j)
                    str = adjustl(str)
                    if ( present(param_name_latex) ) param_name_latex = '\sigma_{8}('//trim(str)//')'
                    return
                end if
            end do
        end if
        ! add Hubble if wanted:
        if ( FP%fisher_der%want_loghubble ) then
            do i=1, FP%fisher_der%FD_num_redshift
                ind = ind + 1
                if ( param_number == ind )  then
                    z = FP%fisher_der%FD_redshift(i)
                    write (str, *) i
                    str = adjustl(str)
                    param_name = 'logHoz_'//trim(str)
                    write (str, '(F20.2)') z
                    str = adjustl(str)
                    if ( present(param_name_latex) ) param_name_latex = '\log_{10}H('//trim(str)//')'
                    return
                end if
            end do
        end if
        ! add DA if wanted:
        if ( FP%fisher_der%want_logDA ) then
            do i=1, FP%fisher_der%FD_num_redshift
                ind = ind + 1
                if ( param_number == ind )  then
                    z = FP%fisher_der%FD_redshift(i)
                    write (str, *) i
                    str = adjustl(str)
                    param_name = 'logDAz_'//trim(str)
                    write (str, '(F20.2)') z
                    str = adjustl(str)
                    if ( present(param_name_latex) ) param_name_latex = '\log_{10}D_{\rm A}('//trim(str)//')'
                    return
                end if
            end do
        end if


    end subroutine Fisher_derived_param_names

    ! ---------------------------------------------------------------------------------------------
    !> This subroutine computes the derivative of derived parameters.
    !! The calculation is carried out with a fixed stepsize finite difference formula or by
    !! the Richardson’s deferred approach to the limit.
    subroutine Derived_Parameter_Derivative( P, FP, param_num, &
        & derived_initial, derived_derivative, derived_der_error, &
        & initial_step, &
        & num_param, &
        & error_code  )

        implicit none

        Type(CAMBparams)       , intent(in) :: P       !< Input CAMBparams
        Type(cosmicfish_params) , intent(in) :: FP      !< Input Cosmicfish params

        integer , intent(in)  :: param_num  !< Input number of the parameter wrt the derivative is computed
        integer , intent(in)  :: num_param  !< Input number of parameters involved in the calculation

        real(dl), dimension(:), intent(in)  :: derived_initial     !< Derived parameters computed for h=0
        real(dl), dimension(:), intent(out) :: derived_derivative  !< Output returns the derivative of the derived parameter vector
        real(dl), dimension(:), intent(out) :: derived_der_error   !< Output returns the error on the derived parameter derivative

        real(dl), intent(in)  :: initial_step  !< Input initial value of h. Does not need to be small for the adaptive algorithm.
        integer , intent(out) :: error_code    !< Output error code. When run the code will not stop on error but report it
                                               !<  0 = everything was fine computation exited correctly
                                               !<  1 = error on input parameters
                                               !<  2 = error on computation
                                               !<  3 = error on quality of the results

        ! definitions:
        real(dl), dimension(:), allocatable :: der_temp_1, der_temp_2, errt
        integer , dimension(:), allocatable :: accept_mask
        real(dl), dimension(:,:,:), allocatable :: a
        real(dl), dimension(num_param) :: param_vector, fiducial_param_vector, temp_param_vector
        ! parameters of the algorithm:
        real(dl), parameter :: con  = 1.4_dl  ! parameter that decides the decrease in the stepsize h
        real(dl), parameter :: big  = 1.d+30  ! a big number used to fill the error array initially
        real(dl), parameter :: safe = 2._dl   ! parameter used to define the stopping cryterium
        integer , parameter :: ntab = 30      ! maximum number of Cls computation
        real(dl) :: con2 = con*con
        integer  :: i, j, k, l, m, side, err, AllocateStatus
        real(dl) :: fac, hh
        real(dl) :: temp1, temp2, temp3
        Type(CAMBparams) :: P_temp
        Type(cosmicfish_params) :: FP_temp
        integer  :: accepted
        real(dl) :: acceptance
        integer  :: derived_param_size

        ! initialize output variables:
        error_code         = 0
        derived_derivative = 0._dl
        derived_der_error  = 0._dl

        ! check input values are correct
        if ( initial_step == 0._dl ) then
            stop 'Derived_Parameter_Derivative initial stepsize is zero'
        end if

        call dimension_derived_parameters( P, FP, derived_param_size )

        ! allocate the arrays:
        if ( FP%adaptivity ) then
            allocate( der_temp_1(derived_param_size),  &
                & der_temp_2(derived_param_size),  &
                & errt(derived_param_size),       &
                & accept_mask(derived_param_size),&
                & stat = AllocateStatus )
            if (AllocateStatus /= 0) stop "Derived_Parameter_Derivative: Allocation failed. Not enough memory"
            allocate( a(ntab, ntab, derived_param_size) , stat = AllocateStatus )
            if (AllocateStatus /= 0) stop "Derived_Parameter_Derivative: Allocation failed. Not enough memory for the derivative array"
        else
            allocate( der_temp_1(derived_param_size),  &
                & der_temp_2(derived_param_size),  &
                & stat = AllocateStatus )
            if (AllocateStatus /= 0) stop "Derived_Parameter_Derivative: Allocation failed. Not enough memory"
            allocate( a(1,1, derived_param_size) , stat = AllocateStatus )
            if (AllocateStatus /= 0) stop "Derived_Parameter_Derivative: Allocation failed. Not enough memory for the derivative array"
        end if

        ! initialize the parameter array:
        P_temp  = P
        FP_temp = FP

        call CAMB_params_to_params_array( P_temp, FP, num_param, param_vector )

        fiducial_param_vector = param_vector
        temp_param_vector     = param_vector

        ! initialize the computation:
        hh                = initial_step
        derived_der_error = big
        side              = 0 ! both sides
        if ( FP%adaptivity ) accept_mask = 0

        ! do the initial step on the right:
        param_vector(param_num)      = fiducial_param_vector(param_num) + hh
        call params_array_to_CAMB_param( P_temp, FP_temp, num_param, param_vector )
        call Compute_Derived_Parameters( P_temp, FP_temp, der_temp_1, derived_param_size, err )
        ! if on the right we cannot go we go to the left:
        if ( err /= 0 ) side = 1
        ! do the initial step on the left:
        temp_param_vector(param_num) = fiducial_param_vector(param_num) - hh
        call params_array_to_CAMB_param( P_temp, FP_temp, num_param, temp_param_vector )
        call Compute_Derived_Parameters( P_temp, FP_temp, der_temp_2, derived_param_size, err )

        ! if on the left we cannot go decide what to do:
        if ( err /= 0 ) then
            if ( side == 1 ) then ! we cannot go to the left nor to the right. Return derivative zero with error code.
                error_code           = 4
                derived_der_error    = 0.0_dl
                derived_derivative   = 0.0_dl
                return
            else if ( side == 0 ) then ! we cannot go to the left but we can go to the right:
                side = 2
            end if
        end if

        ! store the results based on the side we decided to go:
        if (side == 0) then ! both sides
            a(1,1,:) = (der_temp_1(:) - der_temp_2(:))/(2.0_dl*hh)
        else if (side == 1) then ! left side: increment negative
            a(1,1,:) = (der_temp_2(:) - derived_initial(:))/(-hh)
        else if (side == 2) then ! right side: increment positive
            a(1,1,:) = (der_temp_1(:) - derived_initial(:))/(hh)
        end if

        ! fixed stepsize derivative if adaptivity is not wanted:
        if ( .not. FP%adaptivity ) then
            derived_derivative(:) = a(1,1,:)
            derived_der_error     = 0.0_dl
            deallocate( der_temp_1, der_temp_2, a )
            return
        end if

        ! start the derivative computation:
        do i = 2, ntab

            ! feedback for the adaptive derivative:
            if ( FP%cosmicfish_feedback >= 2 ) then
                accepted   = sum(accept_mask)
                acceptance = REAL(accepted)/REAL(derived_param_size)
                print*, 'iteration:', i
                print*, '  stepsize:', hh
                print*, '  accepted elements', accepted, 'accepted ratio', acceptance
            end if

            ! decrease the stepsize:
            hh = hh/con

            ! compute the finite difference:
            if (side == 0) then
                param_vector(param_num)      = fiducial_param_vector(param_num) + hh
                temp_param_vector(param_num) = fiducial_param_vector(param_num) - hh
                ! compute for + hh
                call params_array_to_CAMB_param( P_temp, FP_temp, num_param, param_vector )
                call Compute_Derived_Parameters( P_temp, FP_temp, der_temp_1, derived_param_size, err )
                if ( err /= 0 ) then
                    error_code     = 4
                    derived_der_error  = 0.0_dl
                    derived_derivative = 0.0_dl
                    return
                end if
                ! compute for - hh
                call params_array_to_CAMB_param( P_temp, FP_temp, num_param, temp_param_vector )
                call Compute_Derived_Parameters( P_temp, FP_temp, der_temp_2, derived_param_size, err )
                if ( err /= 0 ) then
                    error_code     = 4
                    derived_der_error  = 0.0_dl
                    derived_derivative = 0.0_dl
                    return
                end if
                ! store the result
                a(1,i,:) = (der_temp_1(:) - der_temp_2(:))/(2.0_dl*hh)
            else if (side==1) then
                param_vector(param_num)      = fiducial_param_vector(param_num) - hh
                call params_array_to_CAMB_param( P_temp, FP_temp, num_param, param_vector )
                call Compute_Derived_Parameters( P_temp, FP_temp, der_temp_1, derived_param_size, err )
                if ( err /= 0 ) then
                    error_code     = 4
                    derived_der_error  = 0.0_dl
                    derived_derivative = 0.0_dl
                    return
                end if
                ! store the result
                a(1,i,:) = (der_temp_1(:) - derived_initial(:))/(-hh)
            else if (side==2) then
                param_vector(param_num)      = fiducial_param_vector(param_num) + hh
                call params_array_to_CAMB_param( P_temp, FP_temp, num_param, param_vector )
                call Compute_Derived_Parameters( P_temp, FP_temp, der_temp_1, derived_param_size, err )
                if ( err /= 0 ) then
                    error_code     = 4
                    derived_der_error  = 0.0_dl
                    derived_derivative = 0.0_dl
                    return
                end if
                ! store the result
                a(1,i,:) = (der_temp_1(:) - derived_initial(:))/(hh)
            end if

            ! now do the extrapolation table:
            fac=CON2
            do j = 2, i

                forall ( k = 1:derived_param_size, accept_mask(k) == 0 )
                    ! compute the interpolation table:
                    a(j,i,k) = (a(j-1,i,k)*fac-a(j-1,i-1,k))/(fac-1._dl)
                    ! estimate the error:
                    errt(k)  = max( abs(a(j,i,k)-a(j-1,i,k)), &
                        & abs(a(j,i,k)-a(j-1,i-1,k)))
                end forall

                fac=CON2*fac

                ! if there is an improvement save the result:

                forall ( k = 1:derived_param_size, &
                    & errt(k) .le. derived_der_error(k) .and. accept_mask(k) == 0 )

                    derived_der_error(k)  = errt(k)
                    derived_derivative(k) = a(j,i,k)
                end forall

            end do
            ! accept a value of the derivative if the value of the extrapolation is stable
            ! Store this information into the mask matrix
            ! so that the value contained in cl_error and cl_derivative is not overwritten.
            forall ( k = 1:derived_param_size, &
                & accept_mask(k) == 0 .and. &
                & (a(i,i,k) - a(i-1,i-1,k)) .ge. safe*derived_der_error(k) )

                accept_mask(k) = 1
            end forall

            ! if all the entries match the required accuracy exit the loop
            if ( sum(accept_mask) == derived_param_size ) exit

        end do

        ! quality check on the results:
        if ( sum(accept_mask) /= derived_param_size ) then
            error_code = 3
        end if

        ! deallocate the arrays:
        if ( FP%adaptivity ) then
            deallocate( der_temp_1, der_temp_2, errt, accept_mask, a )
        end if

        return

    end subroutine Derived_Parameter_Derivative

    ! ---------------------------------------------------------------------------------------------
    !> This subroutine computes the Fisher matrix for derived parameters.
    subroutine Fisher_derived( P_in, FP_in, num_param, num_derived, Derived_Matrix, outroot )

        implicit none

        Type(CAMBparams)       , intent(in) :: P_in       !< Input CAMBparams
        Type(cosmicfish_params) , intent(in) :: FP_in      !< Input Cosmicfish params

        integer , intent(in)  :: num_param                                         !< Input number of parameters
        integer , intent(in)  :: num_derived                                       !< Input number of derived parameters
        real(dl), dimension(num_param,num_derived), intent(out) :: Derived_Matrix  !< Output Fisher matrix
        character(len=*), intent(in), optional :: outroot                          !< Optional input: filename that will be used to dump to file optional output if feedback is greater than 1.

        Type(CAMBparams)       :: P
        Type(cosmicfish_params) :: FP

        real(dl) :: time_1, time_2
        real(dl) :: initial_step, temp, temp1, temp2
        integer  :: AllocateStatus, err, ind, ind2
        integer  :: i, j, k, l, number
        real(dl), allocatable :: param_array(:)
        real(dl), allocatable, dimension(:)   :: derived_fiducial, derived_derivative, derived_temp
        real(dl), allocatable, dimension(:,:) :: derived_derivative_array

        ! copy input parameters:
        P  = P_in
        FP = FP_in

        ! enforce parameter compatibility:
        call derived_fisher_parameter_compatibility( P, FP )

        ! warning against 1 parameter derived:
        if ( num_derived==1 ) then
            write(*,*) 'WARNING: for better compatibility with CosmicFish PyLib use at least two derived parameters.'
        end if

        ! allocate a parameter array:
        allocate( param_array(num_param) )
        ! write the parameters in the parameter array:
        call CAMB_params_to_params_array( P, FP, num_param, param_array )

        allocate( derived_fiducial(num_derived), derived_derivative(num_derived), derived_temp(num_derived), stat = AllocateStatus )
        if (AllocateStatus /= 0) stop "Allocation failed. Not enough memory to allocate derived fiducial."
        allocate( derived_derivative_array( num_param, num_derived ), stat = AllocateStatus)
        if (AllocateStatus /= 0) stop "Allocation failed. Not enough memory to allocate derived_derivative_array"

        ! compute the fiducial:
        time_1 = omp_get_wtime()
        call Compute_Derived_Parameters( P, FP, derived_fiducial, num_derived, err )
        time_2 = omp_get_wtime() - time_1

        if ( err == 0 ) then
        else if ( err == 4 ) then
            write(*,*) 'Fiducial model unstable. Calculation cannot proceed.'
            stop
        else
            write(*,*) 'Error in computing the fiducial. Calculation cannot proceed.'
            stop
        end if

        ! print some feedback
        if ( FP%cosmicfish_feedback >= 1 ) then
            write(*,'(a)')         '**************************************************************'
            write(*,'(a,i10)')         'Number of derived parameters : ', num_derived
            write(*,'(a,i10)')         'Number of parameters         : ', num_param
            write(*,'(a,f10.3,a)')     'Time for derived calculation : ', time_2, ' (s)'
            write(*,'(a,i10)')         'Number of omp processes      : ', OMP_GET_MAX_THREADS()
            if ( FP%adaptivity ) then
                write(*,'(a,f10.3,a)') 'Estimated completion time    : ', time_2*20._dl*num_param, ' (s)'
            else
                write(*,'(a,f10.3,a)') 'Estimated completion time    : ', time_2*2._dl*num_param, ' (s)'
            end if
            write(*,'(a)')         '**************************************************************'
        end if

        ! compute the derivatives of the Cls:
        if ( FP%cosmicfish_feedback >= 1 ) then
            write(*,'(a)') 'Computation of the derived parameters derivatives'
        end if

        time_1 = omp_get_wtime()
        do ind = 1, num_param

            if ( FP%cosmicfish_feedback >= 1 ) then
                write(*,'(a,i5,a,E15.5)') 'Doing parameter number: ', ind, ' value: ', param_array(ind)
            end if

            if ( param_array(ind) /= 0._dl ) then
                initial_step = param_array(ind)*3.0_dl/100._dl
            else
                initial_step = 3.0_dl/100._dl
            end if

            call Derived_Parameter_Derivative( P, FP, ind, &
                & derived_fiducial, derived_derivative, derived_temp, &
                & initial_step, &
                & num_param, &
                & err  )

            if (err /= 0) print*, 'WARNING: something went wrong with the Derived derivative of parameter:', ind, 'Error code:', err

            derived_derivative_array(ind,:) = derived_derivative(:)

        end do
        time_2 = omp_get_wtime() - time_1

        if ( FP%cosmicfish_feedback >= 1 ) then
            write(*,*)
            write(*,'(a,f8.3,a)') 'Time to compute all the derivatives:', time_2, ' (s)'
        end if

        ! compute the Fisher:
        do ind = 1, num_param
            do ind2 = 1, num_derived

                Derived_Matrix(ind, ind2) = derived_derivative_array(ind, ind2)

            end do
        end do
        time_2 = omp_get_wtime() - time_1

        if ( FP%cosmicfish_feedback >= 1 ) then
            write(*,'(a,f8.3,a)') 'Time to compute the derived matrix:', time_2, ' (s)'
        end if

    end subroutine Fisher_derived

    ! ---------------------------------------------------------------------------------------------
    !> This subroutine saves the Fisher derived matrix to file taking care of putting into the file all the
    !> relevant informations.
    subroutine save_derived_Fisher_to_file( P_in, FP_in, num_params, num_derived, matrix, filename )

        implicit none

        Type(CAMBparams)                          :: P_in                   !< Input CAMBparams
        Type(cosmicfish_params)                    :: FP_in                  !< Input Cosmicfish params
        integer , intent(in)                      :: num_params             !< Number of parameters
        integer , intent(in)                      :: num_derived            !< Number of derived parameters
        real(dl), intent(in), dimension(num_params, num_derived) :: matrix  !< Input matrix
        character(len=*), intent(in)              :: filename               !< Name of the file where the Fisher matrix will be saved

        Type(CAMBparams)        :: P
        Type(cosmicfish_params)  :: FP
        integer :: ind, err
        character(LEN=500) :: param_name
        character(LEN=500) :: param_name_latex

        real(dl), allocatable, dimension(:) :: parameters, parameters_derived

        ! copy input parameters:
        P  = P_in
        FP = FP_in
        ! enforce parameter compatibility:
        call derived_fisher_parameter_compatibility( P, FP )

        allocate( parameters(num_params), parameters_derived(num_derived) )
        call CAMB_params_to_params_array( P, FP, num_params, parameters )
        call Compute_Derived_Parameters( P, FP, parameters_derived, num_derived, err )
        if ( err /= 0 ) stop 'Error: Compute_Derived_Parameters failed when save_derived_Fisher_to_file'

        open(unit=666, FILE=filename, ACTION="write", STATUS="replace")

        ! write the header:
        write(666,'(a)') '#'
        write(666,'(a)') '# This file contains a Fisher matrix created with the CosmicFish code.'
        ! write the parameter names:
        write(666,'(a)') '#'
        write(666,'(a)') '# The parameters of this derived parameters Fisher matrix are:'
        write(666,'(a)') '#'

        do ind = 1, num_params
            call Fisher_param_names(P, FP, ind, param_name, param_name_latex)
            write(666,'(a,I5,a,a,a,a,a,E25.16)') '# ', ind, '    ', trim(param_name), '    ', trim(param_name_latex), '    ', parameters(ind)
        end do
        ! write the parameter names:
        write(666,'(a)') '#'
        write(666,'(a)') '# The derived parameters of this matrix are:'
        write(666,'(a)') '#'

        do ind = 1, num_derived
            call Fisher_derived_param_names(P, FP, ind, param_name, param_name_latex)
            write(666,'(a,I5,a,a,a,a,a,E25.16)') '# ', ind, '    ', trim(param_name), '    ', trim(param_name_latex), '    ', parameters_derived(ind)
        end do
        write(666,'(a)') '#'
        write(666,'(a)')

        ! write the Fisher matrix:
        do ind = 1, num_params
            write(666,'(2x,*(E25.16,2x))') matrix(ind, :)
        end do

        close(666)

    end subroutine save_derived_Fisher_to_file

    ! ---------------------------------------------------------------------------------------------
    !> This subroutine saves a file containing the parameter names for the derived Fisher matrix run.
    subroutine save_derived_paramnames_to_file( P_in, FP_in, num_params, num_derived, filename )

        implicit none

        Type(CAMBparams)                          :: P_in                     !< Input CAMBparams
        Type(cosmicfish_params)                    :: FP_in                    !< Input Cosmicfish params
        integer , intent(in)                      :: num_params               !< Number of parameters of the Fisher matrix
        integer , intent(in)                      :: num_derived              !< Number of derived parameters
        character(len=*), intent(in)              :: filename                 !< Name of the file where the Fisher matrix paramnames file will be saved

        Type(CAMBparams)        :: P
        Type(cosmicfish_params)  :: FP
        integer :: ind, err
        character(LEN=500) :: param_name
        character(LEN=500) :: param_name_latex

        real(dl), allocatable, dimension(:) :: parameters, parameters_derived

        ! copy input parameters:
        P  = P_in
        FP = FP_in
        ! enforce parameter compatibility:
        call derived_fisher_parameter_compatibility( P, FP )
        allocate( parameters(num_params), parameters_derived(num_derived) )
        call CAMB_params_to_params_array( P, FP, num_params, parameters )

        call Compute_Derived_Parameters( P, FP, parameters_derived, num_derived, err )

        if ( err /= 0 ) stop 'Error: Compute_Derived_Parameters failed when save_derived_Fisher_to_file'

        open(unit=666, FILE=filename, ACTION="write", STATUS="replace")

        ! write the header:

        write(666,'(a)') '#'
        write(666,'(a)') '# This file contains the parameter names for a derived Fisher matrix.'
        write(666,'(a)') '# Derived parameters are indicated with a * at the end of the name.'
        write(666,'(a)') '#'
        ! write the parameter names:
        do ind = 1, num_params
            call Fisher_param_names(P, FP, ind, param_name, param_name_latex)
            write(666,'(a,a,a,a,E25.16)') trim(param_name), '    ', trim(param_name_latex), '    ', parameters(ind)
        end do
        do ind = 1, num_derived
            call Fisher_derived_param_names(P, FP, ind, param_name, param_name_latex)
            write(666,'(a,a,a,a,E25.16)') trim(param_name)//'*', '    ', trim(param_name_latex), '    ', parameters_derived(ind)
        end do

        close(666)

    end subroutine save_derived_paramnames_to_file

end module Fisher_calculator_Derived
