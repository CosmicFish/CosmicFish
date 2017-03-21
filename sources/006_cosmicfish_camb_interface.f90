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

!> @file 006_cosmicfish_camb_interface.f90
!! This file contains the interface layer between CosmicFish and CAMB. It decide how to handle
!! Fisher parameters and how to call camb appropriately.

!----------------------------------------------------------------------------------------
!> This module contains the interface layer between CosmicFish and CAMB.

!> @author Marco Raveri

module cosmicfish_camb_interface

    use precision
    use cosmicfish_types
    use ModelParams
    use CAMB
    use LambdaGeneral

    implicit none

    private

    public num_param_Fisher, Fisher_param_names, params_array_to_CAMB_param, &
        & CAMB_params_to_params_array, save_Fisher_to_file, save_paramnames_to_file

contains

    ! ---------------------------------------------------------------------------------------------
    !> This subroutine returns the number of parameters in the Fisher matrix based on the input
    !! parameters.
    subroutine num_param_Fisher(P, FP, num_param)

        implicit none

        Type(CAMBparams)      , intent(in)  :: P         !< Input CAMBparams
        Type(cosmicfish_params), intent(in)  :: FP        !< Input Cosmicfish params
        integer               , intent(out) :: num_param !< Output number of parameters

        integer :: i

        num_param = 0

        ! standard parameters:
        if ( FP%fisher_par%want_ombh2           )  num_param = num_param + 1 ! ombh2
        if ( FP%fisher_par%want_omch2           )  num_param = num_param + 1 ! omch2
        if ( FP%fisher_par%want_omnuh2          )  num_param = num_param + 1 ! omnuh2
        if ( FP%fisher_par%want_hubble          )  num_param = num_param + 1 ! small hubble
        if ( FP%fisher_par%want_helium_fraction )  num_param = num_param + 1 ! helium abundance
        if ( FP%fisher_par%want_massless        )  num_param = num_param + 1 ! massless neutrinos
        ! primordial spectrum parameters:
        if ( FP%fisher_par%want_scalar_amp            ) num_param = num_param + 1 ! scalar_amp
        if ( FP%fisher_par%want_scalar_spectral_index ) num_param = num_param + 1 ! scalar_spectral_index
        if ( FP%fisher_par%want_scalar_nrun           ) num_param = num_param + 1 ! scalar_nrun
        if ( FP%fisher_par%want_tensor_spectral_index ) num_param = num_param + 1 ! tensor_spectral_index
        if ( FP%fisher_par%want_initial_ratio         ) num_param = num_param + 1 ! initial_ratio
        ! reionization:
        if ( FP%fisher_par%want_re_optical_depth      ) num_param = num_param + 1 ! reionization optical depth

        ! Nuisance parameters:
        ! bias in every counts redshift windows:
        if ( FP%fisher_par%want_bias ) then
            do i = 1, num_redshiftwindows
                if ( Redshift_w(i)%kind == window_counts ) then
                    num_param = num_param + 1
                end if
            end do
        end if
        ! nuissance parameters for supernovae
        if ( FP%fisher_par%want_alpha_SN ) num_param = num_param + 1 ! alpha_SN
        if ( FP%fisher_par%want_beta_SN  ) num_param = num_param + 1 ! beta_SN
        if ( FP%fisher_par%want_M0_SN    ) num_param = num_param + 1 ! M0_SN

#ifdef COSMICFISH_CAMB
        ! ppf parameters
        if ( FP%fisher_par%want_w0_ppf   ) num_param = num_param + 1 ! w0_ppf
        if ( FP%fisher_par%want_wa_ppf   ) num_param = num_param + 1 ! wa_ppf
        if ( FP%fisher_par%want_cs_ppf   ) num_param = num_param + 1 ! cs_ppf
        if ( FP%fisher_par%want_cT       ) num_param = num_param + 1 ! cT
#endif

#ifdef COSMICFISH_EFTCAMB
        ! EFTCAMB parameters
        ! w_DE parameters:
        if ( ( P%EFTflag == 1 .or. P%EFTflag == 2 .or. P%EFTflag == 3 ) .and. P%EFTwDE /= 0 ) then
            if ( P%EFTwDE == 1 ) then ! DE with constant Eos
                num_param = num_param + 1 ! EFTw0
            else if ( P%EFTwDE == 2 ) then ! CPL parametrization
                num_param = num_param + 2 ! EFTw0, EFTwa
            else if ( P%EFTwDE == 3 ) then ! JBP parametrization
                num_param = num_param + 3 ! EFTw0, EFTwa, EFTwn
            else if ( P%EFTwDE == 4 ) then ! turning point parametrization
                num_param = num_param + 3 ! EFTw0, EFTwa, EFTwat
            else if ( P%EFTwDE == 5 ) then ! Taylor expansion
                num_param = num_param + 4 ! EFTw0, EFTwa, EFTw2, EFTw3
            else if ( P%EFTwDE == 6 ) then ! User defined
                num_param = num_param + 1
            end if
        end if
        ! EFT functions parameters:
        if ( P%EFTflag == 1 ) then ! Pure EFT models
            ! EFTOmega:
            if ( P%PureEFTmodelOmega == 1 ) then
                num_param = num_param + 1
            else if ( P%PureEFTmodelOmega == 2 ) then
                num_param = num_param + 1
            else if ( P%PureEFTmodelOmega == 3 ) then
                num_param = num_param + 2
            else if ( P%PureEFTmodelOmega == 4 ) then
                num_param = num_param + 2
            else if ( P%PureEFTmodelOmega == 5 ) then
                num_param = num_param + 1
            end if
            ! EFTGamma1
            if ( P%PureEFTmodelGamma1 == 1 ) then
                num_param = num_param + 1
            else if ( P%PureEFTmodelGamma1 == 2 ) then
                num_param = num_param + 1
            else if ( P%PureEFTmodelGamma1 == 3 ) then
                num_param = num_param + 2
            else if ( P%PureEFTmodelGamma1 == 4 ) then
                num_param = num_param + 2
            else if ( P%PureEFTmodelGamma1 == 5 ) then
                num_param = num_param + 1
            end if
            ! EFTGamma2
            if ( P%PureEFTmodelGamma2 == 1 ) then
                num_param = num_param + 1
            else if ( P%PureEFTmodelGamma2 == 2 ) then
                num_param = num_param + 1
            else if ( P%PureEFTmodelGamma2 == 3 ) then
                num_param = num_param + 2
            else if ( P%PureEFTmodelGamma2 == 4 ) then
                num_param = num_param + 2
            else if ( P%PureEFTmodelGamma2 == 5 ) then
                num_param = num_param + 1
            end if
            ! EFTGamma3
            if ( P%PureEFTmodelGamma3 == 1 ) then
                num_param = num_param + 1
            else if ( P%PureEFTmodelGamma3 == 2 ) then
                num_param = num_param + 1
            else if ( P%PureEFTmodelGamma3 == 3 ) then
                num_param = num_param + 2
            else if ( P%PureEFTmodelGamma3 == 4 ) then
                num_param = num_param + 2
            else if ( P%PureEFTmodelGamma3 == 5 ) then
                num_param = num_param + 1
            end if
            ! EFTGamma4
            if ( P%PureEFTmodelGamma4 == 1 ) then
                num_param = num_param + 1
            else if ( P%PureEFTmodelGamma4 == 2 ) then
                num_param = num_param + 1
            else if ( P%PureEFTmodelGamma4 == 3 ) then
                num_param = num_param + 2
            else if ( P%PureEFTmodelGamma4 == 4 ) then
                num_param = num_param + 2
            else if ( P%PureEFTmodelGamma4 == 5 ) then
                num_param = num_param + 1
            end if
            ! EFTGamma5
            if ( P%PureEFTmodelGamma5 == 1 ) then
                num_param = num_param + 1
            else if ( P%PureEFTmodelGamma5 == 2 ) then
                num_param = num_param + 1
            else if ( P%PureEFTmodelGamma5 == 3 ) then
                num_param = num_param + 2
            else if ( P%PureEFTmodelGamma5 == 4 ) then
                num_param = num_param + 2
            else if ( P%PureEFTmodelGamma5 == 5 ) then
                num_param = num_param + 1
            end if
            ! EFTGamma6
            if ( P%PureEFTmodelGamma6 == 1 ) then
                num_param = num_param + 1
            else if ( P%PureEFTmodelGamma6 == 2 ) then
                num_param = num_param + 1
            else if ( P%PureEFTmodelGamma6 == 3 ) then
                num_param = num_param + 2
            else if ( P%PureEFTmodelGamma6 == 4 ) then
                num_param = num_param + 2
            else if ( P%PureEFTmodelGamma6 == 5 ) then
                num_param = num_param + 1
            end if
        else if ( P%EFTflag == 2 ) then ! designer EFT models
            if ( P%DesignerEFTmodel == 1 ) then
                num_param = num_param + 1 ! B0
            else if ( P%DesignerEFTmodel == 2 ) then
                num_param = num_param + 0 ! mc 5e models
            end if
        else if ( P%EFTflag == 3 ) then ! EFT alternative parametrizations
            if ( P%AltParEFTmodel == 1 ) then
                if ( P%RPHmassPmodel == 1 ) then
                    num_param = num_param + 1
                else if ( P%RPHmassPmodel == 2 ) then
                    num_param = num_param + 2
                end if
                if ( P%RPHkineticitymodel == 1 ) then
                    num_param = num_param + 1
                else if ( P%RPHkineticitymodel == 2 ) then
                    num_param = num_param + 2
                end if
                if ( P%RPHbraidingmodel == 1 ) then
                    num_param = num_param + 1
                else if ( P%RPHbraidingmodel == 2 ) then
                    num_param = num_param + 2
                end if
                if ( P%RPHtensormodel == 1 ) then
                    num_param = num_param + 1
                else if ( P%RPHtensormodel == 2 ) then
                    num_param = num_param + 2
                end if
            end if
        else if ( P%EFTflag == 4 ) then ! Full mapping EFT model
            if ( P%FullMappingEFTmodel == 1 ) then
                if ( P%HoravaSolarSystem .eqv. .True. ) then
                    num_param = num_param + 2
                else
                    num_param = num_param + 3
                end if
            end if
        end if
#endif

#ifdef COSMICFISH_MGCAMB
        if ( P%MGC_model == 1 ) then
            num_param = num_param + 5
        else if ( P%MGC_model == 2 ) then
            num_param = num_param + 2
        else if ( P%MGC_model == 3 ) then
            num_param = num_param + 3
        else if ( P%MGC_model == 4 ) then
            num_param = num_param + 1
        else if ( P%MGC_model ==5 ) then
            num_param = num_param + 3
        else if ( P%MGC_model ==6 ) then
            num_param = num_param + 1
        else if ( P%MGC_model == 7 ) then
            num_param = num_param + 3
        else if ( P%MGC_model == 8 ) then
            num_param = num_param + 4
        else if ( P%MGC_model == 9 ) then
            num_param = num_param + 2
        else if ( P%MGC_model ==10 ) then
            num_param = num_param + 2
        else if ( P%MGC_model ==11 ) then
            num_param = num_param + 2
            if ( FP%fisher_par%want_c1 ) num_param = num_param + 1 ! c1
            if ( FP%fisher_par%want_c2 ) num_param = num_param + 1 ! c2
            if ( FP%fisher_par%want_lambda ) num_param = num_param + 1 ! lambda
        else if ( P%MGC_model ==12 ) then
            num_param = num_param + 4
            if ( FP%fisher_par%want_c1 ) num_param = num_param + 1 ! c1
            if ( FP%fisher_par%want_c2 ) num_param = num_param + 1 ! c2
            if ( FP%fisher_par%want_lambda ) num_param = num_param + 1 ! lambda
        end if
#endif

    end subroutine num_param_Fisher

    ! ---------------------------------------------------------------------------------------------
    !> This subroutine returns the name corresponding to a parameter number.
    subroutine Fisher_param_names(P, FP, param_number, param_name, param_name_latex)

        implicit none

        Type(CAMBparams)      , intent(in)  :: P            !< Input CAMBparams
        Type(cosmicfish_params), intent(in)  :: FP           !< Input Cosmicfish params
        integer               , intent(in)  :: param_number !< Input number of parameter for which the name is wanted
        character(*)          , intent(out) :: param_name   !< Output name of the parameter
        character(*)          , intent(out), optional :: param_name_latex   !< Output name of the parameter in LaTeX format

        integer :: i, num_param
        character(len=20) str

        num_param = 0

        ! standard parameters:
        if ( FP%fisher_par%want_ombh2           )  num_param  = num_param + 1 ! ombh2
        if ( num_param == param_number          )  then; param_name = 'omegabh2'; if ( present(param_name_latex) ) param_name_latex = '\Omega_b h^2'; return ;
        end if

        if ( FP%fisher_par%want_omch2           )  num_param  = num_param + 1 ! omch2
        if ( num_param == param_number          )  then; param_name = 'omegach2'; if ( present(param_name_latex) ) param_name_latex = '\Omega_c h^2'; return ;
        end if

        if ( FP%fisher_par%want_omnuh2          )  num_param  = num_param + 1 ! omnuh2
        if ( num_param == param_number          )  then; param_name = 'omeganuh2'; if ( present(param_name_latex) ) param_name_latex = '\Omega_\nu h^2'; return ;
        end if

        if ( FP%fisher_par%want_hubble          )  num_param  = num_param + 1 ! omegav or small hubble
        if ( num_param == param_number          )  then; param_name = 'h'; if ( present(param_name_latex) ) param_name_latex = 'h'; return ;
        end if

        if ( FP%fisher_par%want_helium_fraction )  num_param  = num_param + 1 ! helium abundance
        if ( num_param == param_number          )  then; param_name = 'yhe'; if ( present(param_name_latex) ) param_name_latex = 'Y_{He}'; return ;
        end if

        if ( FP%fisher_par%want_massless        )  num_param  = num_param + 1 ! massless neutrinos
        if ( num_param == param_number          )  then; param_name = 'numless'; if ( present(param_name_latex) ) param_name_latex = '\nu_{mless}'; return ;
        end if

        ! primordial spectrum parameters:
        if ( FP%fisher_par%want_scalar_amp            )  num_param  = num_param + 1 ! scalar_amp
        if ( num_param == param_number                )  then; param_name = 'logA'; if ( present(param_name_latex) ) param_name_latex = '{\rm{ln}}(10^{10} A_s)'; return ;
        end if

        if ( FP%fisher_par%want_scalar_spectral_index )  num_param  = num_param + 1 ! scalar_spectral_index
        if ( num_param == param_number                )  then; param_name = 'ns'; if ( present(param_name_latex) ) param_name_latex = 'n_s'; return ;
        end if

        if ( FP%fisher_par%want_scalar_nrun           )  num_param  = num_param + 1 ! scalar_nrun
        if ( num_param == param_number                )  then; param_name = 'nrun'; if ( present(param_name_latex) ) param_name_latex = 'n_{\rm run}'; return ;
        end if

        if ( FP%fisher_par%want_tensor_spectral_index )  num_param  = num_param + 1 ! tensor_spectral_index
        if ( num_param == param_number                )  then; param_name = 'nt'; if ( present(param_name_latex) ) param_name_latex = 'n_t'; return ;
        end if

        if ( FP%fisher_par%want_initial_ratio         )  num_param  = num_param + 1 ! initial_ratio
        if ( num_param == param_number                )  then; param_name = 'r'; if ( present(param_name_latex) ) param_name_latex = 'r'; return ;
        end if

        ! reionization:
        if ( FP%fisher_par%want_re_optical_depth      )  num_param  = num_param + 1 ! reionization optical depth
        if ( num_param == param_number                )  then; param_name = 'tau'; if ( present(param_name_latex) ) param_name_latex = '\tau'; return ;
        end if

        ! Nuisance parameters:
        ! bias in every counts redshift windows:
        if ( FP%fisher_par%want_bias ) then
            do i = 1, num_redshiftwindows
                if ( Redshift_w(i)%kind == window_counts ) then
                    num_param = num_param + 1
                    write (str, *) i
                    str = adjustl(str)
                    if ( num_param == param_number ) then; param_name = 'Bias_W_'//trim(str) ; if ( present(param_name_latex) ) param_name_latex = 'b_{'//trim(str)//'}'; return ;
                    end if
                end if
            end do
        end if

        ! nuissance parameters for supernovae
        if ( FP%fisher_par%want_alpha_SN )  num_param  = num_param + 1 ! alpha_SN
        if ( num_param == param_number   )  then; param_name = 'alpha_SN'; if ( present(param_name_latex) ) param_name_latex = '\alpha_{\rm SN}'; return ;
        end if
        if ( FP%fisher_par%want_beta_SN  )  num_param  = num_param + 1 ! beta_SN
        if ( num_param == param_number   )  then; param_name = 'beta_SN'; if ( present(param_name_latex) ) param_name_latex = '\beta_{\rm SN}'; return ;
        end if
        if ( FP%fisher_par%want_M0_SN    )  num_param  = num_param + 1 ! M0_SN
        if ( num_param == param_number   )  then; param_name = 'M0_SN'; if ( present(param_name_latex) ) param_name_latex = 'M_0^{\rm SN}'; return ;
        end if

#ifdef COSMICFISH_CAMB
        ! ppf parameters
        if ( FP%fisher_par%want_w0_ppf   ) num_param = num_param + 1 ! w0_ppf
        if ( num_param == param_number   )  then; param_name = 'w0_ppf'; if ( present(param_name_latex) ) param_name_latex = 'w_{0}^{\rm ppf}'; return ;
        end if
        if ( FP%fisher_par%want_wa_ppf   ) num_param = num_param + 1 ! wa_ppf
        if ( num_param == param_number   )  then; param_name = 'wa_ppf'; if ( present(param_name_latex) ) param_name_latex = 'w_{a}^{\rm ppf}'; return ;
        end if
        if ( FP%fisher_par%want_cs_ppf   ) num_param = num_param + 1 ! cs_ppf
        if ( num_param == param_number   )  then; param_name = 'cs_ppf'; if ( present(param_name_latex) ) param_name_latex = 'c_{s}^{\rm ppf}'; return ;
        end if
        ! speed of gravitational waves:
        if ( FP%fisher_par%want_cT       ) num_param = num_param + 1 ! cT
        if ( num_param == param_number   )  then; param_name = 'cT'; if ( present(param_name_latex) ) param_name_latex = 'c_{\rm GW}^{2}'; return ;
        end if
#endif

#ifdef COSMICFISH_EFTCAMB
        ! EFTCAMB parameters
        ! w_DE parameters:
        if ( ( P%EFTflag == 1 .or. P%EFTflag == 2 .or. P%EFTflag == 3 ) .and. P%EFTwDE /= 0 ) then
            if ( P%EFTwDE == 1 ) then ! DE with constant Eos
                num_param = num_param + 1 ! EFTw0
                if ( num_param == param_number ) then; param_name = 'EFTw0' ; if ( present(param_name_latex) ) param_name_latex = 'w_0'; return ;
                end if
            else if ( P%EFTwDE == 2 ) then ! CPL parametrization
                num_param = num_param + 1 ! EFTw0
                if ( num_param == param_number ) then; param_name = 'EFTw0' ; if ( present(param_name_latex) ) param_name_latex = 'w_0'; return ;
                end if
                num_param = num_param + 1 ! EFTwa
                if ( num_param == param_number ) then; param_name = 'EFTwa' ; if ( present(param_name_latex) ) param_name_latex = 'w_a'; return ;
                end if
            else if ( P%EFTwDE == 3 ) then ! JBP parametrization
                num_param = num_param + 1 ! EFTw0
                if ( num_param == param_number ) then; param_name = 'EFTw0' ; if ( present(param_name_latex) ) param_name_latex = 'w_0'; return ;
                end if
                num_param = num_param + 1 ! EFTwa
                if ( num_param == param_number ) then; param_name = 'EFTwa' ; if ( present(param_name_latex) ) param_name_latex = 'w_a'; return ;
                end if
                num_param = num_param + 1 ! EFTwn
                if ( num_param == param_number ) then; param_name = 'EFTwn' ; if ( present(param_name_latex) ) param_name_latex = 'w_n'; return ;
                end if
            else if ( P%EFTwDE == 4 ) then ! turning point parametrization
                num_param = num_param + 1 ! EFTw0
                if ( num_param == param_number ) then; param_name = 'EFTw0' ; if ( present(param_name_latex) ) param_name_latex = 'w_0'; return ;
                end if
                num_param = num_param + 1 ! EFTwa
                if ( num_param == param_number ) then; param_name = 'EFTwa' ; if ( present(param_name_latex) ) param_name_latex = 'w_a'; return ;
                end if
                num_param = num_param + 1 ! EFTwat
                if ( num_param == param_number ) then; param_name = 'EFTwat' ; if ( present(param_name_latex) ) param_name_latex = 'wa_t'; return ;
                end if
            else if ( P%EFTwDE == 5 ) then ! Taylor expansion
                num_param = num_param + 1 ! EFTw0
                if ( num_param == param_number ) then; param_name = 'EFTw0' ; if ( present(param_name_latex) ) param_name_latex = 'w_0'; return ;
                end if
                num_param = num_param + 1 ! EFTwa
                if ( num_param == param_number ) then; param_name = 'EFTwa' ; if ( present(param_name_latex) ) param_name_latex = 'w_a'; return ;
                end if
                num_param = num_param + 1 ! EFTw2
                if ( num_param == param_number ) then; param_name = 'EFTw2' ; if ( present(param_name_latex) ) param_name_latex = 'w_2'; return ;
                end if
                num_param = num_param + 1 ! EFTw3
                if ( num_param == param_number ) then; param_name = 'EFTw3' ; if ( present(param_name_latex) ) param_name_latex = 'w_3'; return ;
                end if
            else if ( P%EFTwDE == 6 ) then ! User defined

            end if
        end if

        ! EFT functions parameters:
        if ( P%EFTflag == 1 ) then ! Pure EFT models
            ! EFTOmega:
            if ( P%PureEFTmodelOmega == 1 ) then
                num_param = num_param + 1
                if ( num_param == param_number ) then; param_name = 'EFTOmega0' ; if ( present(param_name_latex) ) param_name_latex = '\Omega_0^{\rm EFT}'; return ;
                end if
            else if ( P%PureEFTmodelOmega == 2 ) then
                num_param = num_param + 1
                if ( num_param == param_number ) then; param_name = 'EFTOmega0' ; if ( present(param_name_latex) ) param_name_latex = '\Omega_0^{\rm EFT}'; return ;
                end if
            else if ( P%PureEFTmodelOmega == 3 ) then
                num_param = num_param + 1
                if ( num_param == param_number ) then; param_name = 'EFTOmega0' ; if ( present(param_name_latex) ) param_name_latex = '\Omega_0^{\rm EFT}'; return ;
                end if
                num_param = num_param + 1
                if ( num_param == param_number ) then; param_name = 'EFTOmegaExp' ; if ( present(param_name_latex) ) param_name_latex = 'n^{\rm EFT}'; return ;
                end if
            else if ( P%PureEFTmodelOmega == 4 ) then
                num_param = num_param + 1
                if ( num_param == param_number ) then; param_name = 'EFTOmega0' ; if ( present(param_name_latex) ) param_name_latex = '\Omega_0^{\rm EFT}'; return ;
                end if
                num_param = num_param + 1
                if ( num_param == param_number ) then; param_name = 'EFTOmegaExp' ; if ( present(param_name_latex) ) param_name_latex = 'n^{\rm EFT}'; return ;
                end if
            else if ( P%PureEFTmodelOmega == 5 ) then

            end if
            ! EFTGamma1
            if ( P%PureEFTmodelGamma1 == 1 ) then
                num_param = num_param + 1
                if ( num_param == param_number ) then; param_name = 'EFTGamma10' ; if ( present(param_name_latex) ) param_name_latex = '\gamma_0^{(1) {\rm EFT}}'; return ;
                end if
            else if ( P%PureEFTmodelGamma1 == 2 ) then
                num_param = num_param + 1
                if ( num_param == param_number ) then; param_name = 'EFTGamma10' ; if ( present(param_name_latex) ) param_name_latex = '\gamma_0^{(1) {\rm EFT}}'; return ;
                end if
            else if ( P%PureEFTmodelGamma1 == 3 ) then
                num_param = num_param + 1
                if ( num_param == param_number ) then; param_name = 'EFTGamma10' ; if ( present(param_name_latex) ) param_name_latex = '\gamma_0^{(1) {\rm EFT}}'; return ;
                end if
                num_param = num_param + 1
                if ( num_param == param_number ) then; param_name = 'EFTGamma1Exp' ; if ( present(param_name_latex) ) param_name_latex = '\gamma_{\rm exp}^{(1) {\rm EFT}}'; return ;
                end if
            else if ( P%PureEFTmodelGamma1 == 4 ) then
                num_param = num_param + 1
                if ( num_param == param_number ) then; param_name = 'EFTGamma10' ; if ( present(param_name_latex) ) param_name_latex = '\gamma_0^{(1) {\rm EFT}}'; return ;
                end if
                num_param = num_param + 1
                if ( num_param == param_number ) then; param_name = 'EFTGamma1Exp' ; if ( present(param_name_latex) ) param_name_latex = '\gamma_{\rm exp}^{(1) {\rm EFT}}'; return ;
                end if
            else if ( P%PureEFTmodelGamma1 == 5 ) then

            end if
            ! EFTGamma2
            if ( P%PureEFTmodelGamma2 == 1 ) then
                num_param = num_param + 1
                if ( num_param == param_number ) then; param_name = 'EFTGamma20' ; if ( present(param_name_latex) ) param_name_latex = '\gamma_0^{(2) {\rm EFT}}'; return ;
                end if
            else if ( P%PureEFTmodelGamma2 == 2 ) then
                num_param = num_param + 1
                if ( num_param == param_number ) then; param_name = 'EFTGamma20' ; if ( present(param_name_latex) ) param_name_latex = '\gamma_0^{(2) {\rm EFT}}'; return ;
                end if
            else if ( P%PureEFTmodelGamma2 == 3 ) then
                num_param = num_param + 1
                if ( num_param == param_number ) then; param_name = 'EFTGamma20' ; if ( present(param_name_latex) ) param_name_latex = '\gamma_0^{(2) {\rm EFT}}'; return ;
                end if
                num_param = num_param + 1
                if ( num_param == param_number ) then; param_name = 'EFTGamma2Exp' ; if ( present(param_name_latex) ) param_name_latex = '\gamma_{\rm exp}^{(2) {\rm EFT}}'; return ;
                end if
            else if ( P%PureEFTmodelGamma2 == 4 ) then
                num_param = num_param + 1
                if ( num_param == param_number ) then; param_name = 'EFTGamma20' ; if ( present(param_name_latex) ) param_name_latex = '\gamma_0^{(2) {\rm EFT}}'; return ;
                end if
                num_param = num_param + 1
                if ( num_param == param_number ) then; param_name = 'EFTGamma2Exp' ; if ( present(param_name_latex) ) param_name_latex = '\gamma_{\rm exp}^{(2) {\rm EFT}}'; return ;
                end if
            else if ( P%PureEFTmodelGamma2 == 5 ) then

            end if
            ! EFTGamma3
            if ( P%PureEFTmodelGamma3 == 1 ) then
                num_param = num_param + 1
                if ( num_param == param_number ) then; param_name = 'EFTGamma30' ; if ( present(param_name_latex) ) param_name_latex = '\gamma_0^{(3) {\rm EFT}}'; return ;
                end if
            else if ( P%PureEFTmodelGamma3 == 2 ) then
                num_param = num_param + 1
                if ( num_param == param_number ) then; param_name = 'EFTGamma30' ; if ( present(param_name_latex) ) param_name_latex = '\gamma_0^{(3) {\rm EFT}}'; return ;
                end if
            else if ( P%PureEFTmodelGamma3 == 3 ) then
                num_param = num_param + 1
                if ( num_param == param_number ) then; param_name = 'EFTGamma30' ; if ( present(param_name_latex) ) param_name_latex = '\gamma_0^{(3) {\rm EFT}}'; return ;
                end if
                num_param = num_param + 1
                if ( num_param == param_number ) then; param_name = 'EFTGamma3Exp' ; if ( present(param_name_latex) ) param_name_latex = '\gamma_{\rm exp}^{(3) {\rm EFT}}'; return ;
                end if
            else if ( P%PureEFTmodelGamma3 == 4 ) then
                num_param = num_param + 1
                if ( num_param == param_number ) then; param_name = 'EFTGamma30' ; if ( present(param_name_latex) ) param_name_latex = '\gamma_0^{(3) {\rm EFT}}'; return ;
                end if
                num_param = num_param + 1
                if ( num_param == param_number ) then; param_name = 'EFTGamma3Exp' ; if ( present(param_name_latex) ) param_name_latex = '\gamma_{\rm exp}^{(3) {\rm EFT}}'; return ;
                end if
            else if ( P%PureEFTmodelGamma3 == 5 ) then

            end if
            ! EFTGamma4
            if ( P%PureEFTmodelGamma4 == 1 ) then
                num_param = num_param + 1
                if ( num_param == param_number ) then; param_name = 'EFTGamma40' ; if ( present(param_name_latex) ) param_name_latex = '\gamma_0^{(4) {\rm EFT}}'; return ;
                end if
            else if ( P%PureEFTmodelGamma4 == 2 ) then
                num_param = num_param + 1
                if ( num_param == param_number ) then; param_name = 'EFTGamma40' ; if ( present(param_name_latex) ) param_name_latex = '\gamma_0^{(4) {\rm EFT}}'; return ;
                end if
            else if ( P%PureEFTmodelGamma4 == 3 ) then
                num_param = num_param + 1
                if ( num_param == param_number ) then; param_name = 'EFTGamma40' ; if ( present(param_name_latex) ) param_name_latex = '\gamma_0^{(4) {\rm EFT}}'; return ;
                end if
                num_param = num_param + 1
                if ( num_param == param_number ) then; param_name = 'EFTGamma4Exp' ; if ( present(param_name_latex) ) param_name_latex = '\gamma_{\rm exp}^{(4) {\rm EFT}}'; return ;
                end if
            else if ( P%PureEFTmodelGamma4 == 4 ) then
                num_param = num_param + 1
                if ( num_param == param_number ) then; param_name = 'EFTGamma40' ; if ( present(param_name_latex) ) param_name_latex = '\gamma_0^{(4) {\rm EFT}}'; return ;
                end if
                num_param = num_param + 1
                if ( num_param == param_number ) then; param_name = 'EFTGamma4Exp' ; if ( present(param_name_latex) ) param_name_latex = '\gamma_{\rm exp}^{(4) {\rm EFT}}'; return ;
                end if
            else if ( P%PureEFTmodelGamma4 == 5 ) then

            end if
            ! EFTGamma5
            if ( P%PureEFTmodelGamma5 == 1 ) then
                num_param = num_param + 1
                if ( num_param == param_number ) then; param_name = 'EFTGamma50' ; if ( present(param_name_latex) ) param_name_latex = '\gamma_0^{(5) {\rm EFT}}'; return ;
                end if
            else if ( P%PureEFTmodelGamma5 == 2 ) then
                num_param = num_param + 1
                if ( num_param == param_number ) then; param_name = 'EFTGamma50' ; if ( present(param_name_latex) ) param_name_latex = '\gamma_0^{(5) {\rm EFT}}'; return ;
                end if
            else if ( P%PureEFTmodelGamma5 == 3 ) then
                num_param = num_param + 1
                if ( num_param == param_number ) then; param_name = 'EFTGamma50' ; if ( present(param_name_latex) ) param_name_latex = '\gamma_0^{(5) {\rm EFT}}'; return ;
                end if
                num_param = num_param + 1
                if ( num_param == param_number ) then; param_name = 'EFTGamma5Exp' ; if ( present(param_name_latex) ) param_name_latex = '\gamma_{\rm exp}^{(5) {\rm EFT}}'; return ;
                end if
            else if ( P%PureEFTmodelGamma5 == 4 ) then
                num_param = num_param + 1
                if ( num_param == param_number ) then; param_name = 'EFTGamma50' ; if ( present(param_name_latex) ) param_name_latex = '\gamma_0^{(5) {\rm EFT}}'; return ;
                end if
                num_param = num_param + 1
                if ( num_param == param_number ) then; param_name = 'EFTGamma5Exp' ; if ( present(param_name_latex) ) param_name_latex = '\gamma_{\rm exp}^{(5) {\rm EFT}}'; return ;
                end if
            else if ( P%PureEFTmodelGamma5 == 5 ) then

            end if
            ! EFTGamma6
            if ( P%PureEFTmodelGamma6 == 1 ) then
                num_param = num_param + 1
                if ( num_param == param_number ) then; param_name = 'EFTGamma60' ; if ( present(param_name_latex) ) param_name_latex = '\gamma_0^{(6) {\rm EFT}}'; return ;
                end if
            else if ( P%PureEFTmodelGamma6 == 2 ) then
                num_param = num_param + 1
                if ( num_param == param_number ) then; param_name = 'EFTGamma60' ; if ( present(param_name_latex) ) param_name_latex = '\gamma_0^{(6) {\rm EFT}}'; return ;
                end if
            else if ( P%PureEFTmodelGamma6 == 3 ) then
                num_param = num_param + 1
                if ( num_param == param_number ) then; param_name = 'EFTGamma60' ; if ( present(param_name_latex) ) param_name_latex = '\gamma_0^{(6) {\rm EFT}}'; return ;
                end if
                num_param = num_param + 1
                if ( num_param == param_number ) then; param_name = 'EFTGamma6Exp' ; if ( present(param_name_latex) ) param_name_latex = '\gamma_{\rm exp}^{(6) {\rm EFT}}'; return ;
                end if
            else if ( P%PureEFTmodelGamma6 == 4 ) then
                num_param = num_param + 1
                if ( num_param == param_number ) then; param_name = 'EFTGamma60' ; if ( present(param_name_latex) ) param_name_latex = '\gamma_0^{(6) {\rm EFT}}'; return ;
                end if
                num_param = num_param + 1
                if ( num_param == param_number ) then; param_name = 'EFTGamma6Exp' ; if ( present(param_name_latex) ) param_name_latex = '\gamma_{\rm exp}^{(6) {\rm EFT}}'; return ;
                end if
            else if ( P%PureEFTmodelGamma6 == 5 ) then

            end if
        else if ( P%EFTflag == 2 ) then ! designer EFT models
            if ( P%DesignerEFTmodel == 1 ) then
                num_param = num_param + 1 ! B0
                if ( num_param == param_number ) then; param_name = 'B0' ; if ( present(param_name_latex) ) param_name_latex = 'B_0'; return ;
                end if
            else if ( P%DesignerEFTmodel == 2 ) then
                num_param = num_param + 0 ! mc 5e models
            end if
        else if ( P%EFTflag == 3 ) then ! EFT alternative parametrizations
            if ( P%AltParEFTmodel == 1 ) then
                if ( P%RPHmassPmodel == 1 ) then
                    num_param = num_param + 1 ! RPHmassP0
                    if ( num_param == param_number ) then; param_name = 'RPHmassP0' ; if ( present(param_name_latex) ) param_name_latex = '\tilde{M_0}'; return ;
                    end if
                else if ( P%RPHmassPmodel == 2 ) then
                    num_param = num_param + 1 ! RPHmassP0
                    if ( num_param == param_number ) then; param_name = 'RPHmassP0' ; if ( present(param_name_latex) ) param_name_latex = '\tilde{M_0}'; return ;
                    end if
                    num_param = num_param + 1 ! RPHmassPexp
                    if ( num_param == param_number ) then; param_name = 'RPHmassPexp' ; if ( present(param_name_latex) ) param_name_latex = '\tilde{M}_{\rm exp}^{\rm RPH}'; return ;
                    end if
                end if
                if ( P%RPHkineticitymodel == 1 ) then
                    num_param = num_param + 1 ! RPHkineticity0
                    if ( num_param == param_number ) then; param_name = 'RPHkineticity0' ; if ( present(param_name_latex) ) param_name_latex = '\alpha_0^{\rm K}'; return ;
                    end if
                else if ( P%RPHkineticitymodel == 2 ) then
                    num_param = num_param + 1 ! RPHkineticity0
                    if ( num_param == param_number ) then; param_name = 'RPHkineticity0' ; if ( present(param_name_latex) ) param_name_latex = '\alpha_0^{\rm K}'; return ;
                    end if
                    num_param = num_param + 1 ! RPHkineticityexp
                    if ( num_param == param_number ) then; param_name = 'RPHkineticityexp' ; if ( present(param_name_latex) ) param_name_latex = '\alpha_{\rm exp}^{\rm K}'; return ;
                    end if
                end if
                if ( P%RPHbraidingmodel == 1 ) then
                    num_param = num_param + 1 ! RPHbraiding0
                    if ( num_param == param_number ) then; param_name = 'RPHbraiding0' ; if ( present(param_name_latex) ) param_name_latex = '\alpha_0^{\rm B}'; return ;
                    end if
                else if ( P%RPHbraidingmodel == 2 ) then
                    num_param = num_param + 1 ! RPHbraiding0
                    if ( num_param == param_number ) then; param_name = 'RPHbraiding0' ; if ( present(param_name_latex) ) param_name_latex = '\alpha_0^{\rm B}'; return ;
                    end if
                    num_param = num_param + 1 ! RPHbraidingexp
                    if ( num_param == param_number ) then; param_name = 'RPHbraidingexp' ; if ( present(param_name_latex) ) param_name_latex = '\alpha_{\rm exp}^{\rm B}'; return ;
                    end if
                end if
                if ( P%RPHtensormodel == 1 ) then
                    num_param = num_param + 1 ! RPHtensor0
                    if ( num_param == param_number ) then; param_name = 'RPHtensor0' ; if ( present(param_name_latex) ) param_name_latex = '\alpha_0^{\rm T}'; return ;
                    end if
                else if ( P%RPHtensormodel == 2 ) then
                    num_param = num_param + 1 ! RPHtensor0
                    if ( num_param == param_number ) then; param_name = 'RPHtensor0' ; if ( present(param_name_latex) ) param_name_latex = '\alpha_0^{\rm T}'; return ;
                    end if
                    num_param = num_param + 1 ! RPHtensorexp
                    if ( num_param == param_number ) then; param_name = 'RPHtensorexp' ; if ( present(param_name_latex) ) param_name_latex = '\alpha_{\rm exp}^{\rm T}'; return ;
                    end if
                end if
            end if
        else if ( P%EFTflag == 4 ) then ! Full mapping EFT model
            if ( P%FullMappingEFTmodel == 1 ) then
                if ( P%HoravaSolarSystem .eqv. .True. ) then
                    num_param = num_param + 1 ! Horava_lambda
                    if ( num_param == param_number ) then; param_name = 'Horava_lambda' ; if ( present(param_name_latex) ) param_name_latex = '{\rm Horava}_{\lambda}'; return ;
                    end if
                    num_param = num_param + 1 ! Horava_eta
                    if ( num_param == param_number ) then; param_name = 'Horava_eta' ; if ( present(param_name_latex) ) param_name_latex = '{\rm Horava}_{\eta}'; return ;
                    end if
                else
                    num_param = num_param + 1 ! Horava_lambda
                    if ( num_param == param_number ) then; param_name = 'Horava_lambda' ; if ( present(param_name_latex) ) param_name_latex = '{\rm Horava}_{\lambda}'; return ;
                    end if
                    num_param = num_param + 1 ! Horava_eta
                    if ( num_param == param_number ) then; param_name = 'Horava_eta' ; if ( present(param_name_latex) ) param_name_latex = '{\rm Horava}_{\eta}'; return ;
                    end if
                    num_param = num_param + 1 ! Horava_xi
                    if ( num_param == param_number ) then; param_name = 'Horava_xi' ; if ( present(param_name_latex) ) param_name_latex = '{\rm Horava}_{\xi}'; return ;
                    end if
                end if
            end if
        end if
#endif

#ifdef COSMICFISH_MGCAMB
        if ( P%MGC_model == 1 ) then
            num_param = num_param + 1 ! B1
            if ( num_param == param_number ) then; param_name = 'B1' ; if ( present(param_name_latex) ) param_name_latex = 'B_{1}'; return ;
            end if
            num_param = num_param + 1 ! B2
            if ( num_param == param_number ) then; param_name = 'B2' ; if ( present(param_name_latex) ) param_name_latex = 'B_{2}'; return ;
            end if
            num_param = num_param + 1 ! lambda1_2
            if ( num_param == param_number ) then; param_name = 'lambda1_2' ; if ( present(param_name_latex) ) param_name_latex = '\lambda_{1}^{2}'; return ;
            end if
            num_param = num_param + 1 ! lambda2_2
            if ( num_param == param_number ) then; param_name = 'lambda2_2' ; if ( present(param_name_latex) ) param_name_latex = '\lambda_{2}^{2}'; return ;
            end if
            num_param = num_param + 1 ! ss
            if ( num_param == param_number ) then; param_name = 'ss' ; if ( present(param_name_latex) ) param_name_latex = 's'; return ;
            end if
        else if ( P%MGC_model == 2 ) then
            num_param = num_param + 1 ! MGQfix
            if ( num_param == param_number ) then; param_name = 'MGQfix' ; if ( present(param_name_latex) ) param_name_latex = 'Q_{0}'; return ;
            end if
            num_param = num_param + 1 ! MGRfix
            if ( num_param == param_number ) then; param_name = 'MGRfix' ; if ( present(param_name_latex) ) param_name_latex = 'R_{0}'; return ;
            end if
        else if ( P%MGC_model == 3 ) then
            num_param = num_param + 1 ! Qnot
            if ( num_param == param_number ) then; param_name = 'Qnot' ; if ( present(param_name_latex) ) param_name_latex = 'Q_{0}'; return ;
            end if
            num_param = num_param + 1 ! Rnot
            if ( num_param == param_number ) then; param_name = 'Rnot' ; if ( present(param_name_latex) ) param_name_latex = 'R_{0}'; return ;
            end if
            num_param = num_param + 1 ! sss
            if ( num_param == param_number ) then; param_name = 'sss' ; if ( present(param_name_latex) ) param_name_latex = 's'; return ;
            end if
        else if ( P%MGC_model == 4 ) then
            num_param = num_param + 1 ! B0
            if ( num_param == param_number ) then; param_name = 'B0' ; if ( present(param_name_latex) ) param_name_latex = 'B_{0}'; return ;
            end if
        else if ( P%MGC_model ==5 ) then
            num_param = num_param + 1 ! B0
            if ( num_param == param_number ) then; param_name = 'B0' ; if ( present(param_name_latex) ) param_name_latex = 'B_{0}'; return ;
            end if
            num_param = num_param + 1 ! beta1
            if ( num_param == param_number ) then; param_name = 'beta1' ; if ( present(param_name_latex) ) param_name_latex = '\beta_{1}'; return ;
            end if
            num_param = num_param + 1 ! s
            if ( num_param == param_number ) then; param_name = 's' ; if ( present(param_name_latex) ) param_name_latex = 's'; return ;
            end if
        else if ( P%MGC_model ==6 ) then
            num_param = num_param + 1 ! Linder_gamma
            if ( num_param == param_number ) then; param_name = 'Linder_gamma' ; if ( present(param_name_latex) ) param_name_latex = '\gamma_{\rm L}'; return ;
            end if
        else if ( P%MGC_model == 7 ) then
            num_param = num_param + 1 ! beta_star
            if ( num_param == param_number ) then; param_name = 'beta_star' ; if ( present(param_name_latex) ) param_name_latex = '\beta_{\star}'; return ;
            end if
            num_param = num_param + 1 ! xi_star
            if ( num_param == param_number ) then; param_name = 'xi_star' ; if ( present(param_name_latex) ) param_name_latex = '\xi_{\star}'; return ;
            end if
            num_param = num_param + 1 ! a_star
            if ( num_param == param_number ) then; param_name = 'a_star' ; if ( present(param_name_latex) ) param_name_latex = 'a_{\star}'; return ;
            end if
        else if ( P%MGC_model == 8 ) then
            num_param = num_param + 1 ! beta0
            if ( num_param == param_number ) then; param_name = 'beta0' ; if ( present(param_name_latex) ) param_name_latex = '\beta_{0}'; return ;
            end if
            num_param = num_param + 1 ! xi0
            if ( num_param == param_number ) then; param_name = 'xi0' ; if ( present(param_name_latex) ) param_name_latex = '\xi_{0}'; return ;
            end if
            num_param = num_param + 1 ! DilR
            if ( num_param == param_number ) then; param_name = 'DilR' ; if ( present(param_name_latex) ) param_name_latex = 'R_{\rm dil}'; return ;
            end if
            num_param = num_param + 1 ! DilS
            if ( num_param == param_number ) then; param_name = 'DilS' ; if ( present(param_name_latex) ) param_name_latex = 'S_{\rm dil}'; return ;
            end if
        else if ( P%MGC_model == 9 ) then
            num_param = num_param + 1 ! F_R0
            if ( num_param == param_number ) then; param_name = 'F_R0' ; if ( present(param_name_latex) ) param_name_latex = 'f_{R,0}'; return ;
            end if
            num_param = num_param + 1 ! FRn
            if ( num_param == param_number ) then; param_name = 'FRn' ; if ( present(param_name_latex) ) param_name_latex = 'n'; return ;
            end if
        else if ( P%MGC_model ==10 ) then
            num_param = num_param + 1 ! beta0
            if ( num_param == param_number ) then; param_name = 'beta0' ; if ( present(param_name_latex) ) param_name_latex = '\beta_0'; return ;
            end if
            num_param = num_param + 1 ! A2
            if ( num_param == param_number ) then; param_name = 'A2' ; if ( present(param_name_latex) ) param_name_latex = 'A_2'; return ;
            end if
        else if ( P%MGC_model ==11 ) then
            num_param = num_param + 1 ! E11
            if ( num_param == param_number ) then; param_name = 'E11' ; if ( present(param_name_latex) ) param_name_latex = 'E_{11}'; return ;
            end if
            num_param = num_param + 1 ! E22
            if ( num_param == param_number ) then; param_name = 'E22' ; if ( present(param_name_latex) ) param_name_latex = 'E_{22}'; return ;
            end if
            if ( FP%fisher_par%want_c1 ) num_param = num_param + 1 ! c1
            if ( num_param == param_number ) then; param_name = 'c1' ; if ( present(param_name_latex) )  param_name_latex = 'c_{1}'; return ;
            end if
            if ( FP%fisher_par%want_c2 ) num_param = num_param + 1 ! c2
            if ( num_param == param_number ) then; param_name = 'c2' ; if ( present(param_name_latex) )  param_name_latex = 'c_{2}'; return ;
            end if
            if ( FP%fisher_par%want_lambda ) num_param = num_param + 1 ! lambda
            if ( num_param == param_number ) then; param_name = 'lambda' ; if ( present(param_name_latex) ) param_name_latex = '\lambda'; return ;
            end if
         else if ( P%MGC_model ==12 ) then
            num_param = num_param + 1 ! E11
            if ( num_param == param_number ) then; param_name = 'E11' ; if ( present(param_name_latex) ) param_name_latex = 'E_{11}'; return ;
            end if
            num_param = num_param + 1 ! E22
            if ( num_param == param_number ) then; param_name = 'E22' ; if ( present(param_name_latex) ) param_name_latex = 'E_{22}'; return ;
            end if
            num_param = num_param + 1 ! E12
            if ( num_param == param_number ) then; param_name = 'E12' ; if ( present(param_name_latex) ) param_name_latex = 'E_{12}'; return ;
            end if
            num_param = num_param + 1 ! E21
            if ( num_param == param_number ) then; param_name = 'E21' ; if ( present(param_name_latex) ) param_name_latex = 'E_{21}'; return ;
            end if
            if ( FP%fisher_par%want_c1 ) num_param = num_param + 1 ! c1
            if ( num_param == param_number ) then; param_name = 'c1' ; if ( present(param_name_latex) )  param_name_latex = 'c_{1}';  return ;
            end if
            if ( FP%fisher_par%want_c2 ) num_param = num_param + 1 ! c2
            if ( num_param == param_number ) then; param_name = 'c2' ; if ( present(param_name_latex) )  param_name_latex = 'c_{2}';  return ;
            end if
            if ( FP%fisher_par%want_lambda ) num_param = num_param + 1 ! lambda
            if ( num_param == param_number ) then; param_name = 'lambda' ; if ( present(param_name_latex) ) param_name_latex = '\lambda'; return ;
            end if
        end if
#endif

    end subroutine Fisher_param_names

    ! ---------------------------------------------------------------------------------------------
    !> This subroutine updates the values of the parameters array based on the values of P and FP.
    subroutine CAMB_params_to_params_array( P, FP, num_param, params_array )

        implicit none

        Type(CAMBparams)               , intent(in)  :: P            !< Input CAMBparams
        Type(cosmicfish_params)         , intent(in)  :: FP           !< Input Cosmicfish params
        integer                        , intent(in)  :: num_param    !< Input length of the parameter array
        real(dl), dimension(num_param) , intent(out) :: params_array !< Output updated parameter array

        integer :: ind, ind2, i, j, k

        ind = 1

        ! standard parameters:
        ! omegab
        if ( FP%fisher_par%want_ombh2 ) then
            params_array(ind) = P%omegab*(P%h0/100._dl)**2
            ind = ind + 1
        end if
        ! omegac
        if ( FP%fisher_par%want_omch2 ) then
            params_array(ind) = P%omegac*(P%h0/100._dl)**2
            ind = ind + 1
        end if
        ! omegan
        if ( FP%fisher_par%want_omnuh2 ) then
            params_array(ind) = P%omegan*(P%h0/100._dl)**2
            ind = ind + 1
        end if
        ! hubble
        if ( FP%fisher_par%want_hubble ) then
            params_array(ind) = P%h0/100._dl
            ind = ind + 1
        end if
        ! helium abundance
        if ( FP%fisher_par%want_helium_fraction ) then
            params_array(ind) = P%yhe
            ind = ind + 1
        end if
        ! massless neutrinos
        if ( FP%fisher_par%want_massless ) then
            params_array(ind) = P%Num_Nu_massless
            ind = ind + 1
        end if
        ! primordial spectrum parameters:
        ! scalar_amp
        if ( FP%fisher_par%want_scalar_amp            ) then
            params_array(ind) = log(10._dl**10) + log(P%InitPower%ScalarPowerAmp(1))
            ind = ind + 1
        end if
        ! scalar_spectral_index
        if ( FP%fisher_par%want_scalar_spectral_index ) then
            params_array(ind) = P%InitPower%an(1)
            ind = ind + 1
        end if
        ! scalar_nrun
        if ( FP%fisher_par%want_scalar_nrun           ) then
            params_array(ind) = P%InitPower%n_run(1)
            ind = ind + 1
        end if
        ! tensor_spectral_index
        if ( FP%fisher_par%want_tensor_spectral_index ) then
            params_array(ind) = P%InitPower%ant(1)
            ind = ind + 1
        end if
        ! initial_ratio
        if ( FP%fisher_par%want_initial_ratio         ) then
            params_array(ind) = P%InitPower%rat(1)
            ind = ind + 1
        end if
        ! reionization: reionization optical depth
        if ( FP%fisher_par%want_re_optical_depth      ) then
            params_array(ind) = P%Reion%optical_depth
            ind = ind + 1
        end if

        ! Nuisance parameters for LSS:
        ! bias in every counts redshift windows:
        if ( FP%fisher_par%want_bias ) then
            ind2 = 0
            do i = 1, num_redshiftwindows
                if ( Redshift_w(i)%kind == window_counts ) then
                    ind2 = ind2 +1
                    params_array(ind) = FP%fisher_cls%bias(ind2)
                    ind = ind +1
                end if
            end do
        end if

        ! nuissance parameters for supernovae:
        ! alpha_SN
        if ( FP%fisher_par%want_alpha_SN ) then
            params_array(ind) = FP%fisher_SN%alpha_SN
            ind = ind +1
        end if
        ! beta_SN
        if ( FP%fisher_par%want_beta_SN  ) then
            params_array(ind) = FP%fisher_SN%beta_SN
            ind = ind +1
        end if
        ! M0_SN
        if ( FP%fisher_par%want_M0_SN    ) then
            params_array(ind) = FP%fisher_SN%M0_SN
            ind = ind +1
        end if

#ifdef COSMICFISH_CAMB
        ! ppf parameters:
        ! w0_ppf
        if ( FP%fisher_par%want_w0_ppf   ) then
            params_array(ind) = P%w_lam
            ind = ind +1
        end if
        ! wa_ppf
        if ( FP%fisher_par%want_wa_ppf   ) then
            params_array(ind) = P%wa_ppf
            ind = ind +1
        end if
        ! cs_ppf
        if ( FP%fisher_par%want_cs_ppf   ) then
            params_array(ind) = P%cs2_lam
            ind = ind +1
        end if
        ! cT
        if ( FP%fisher_par%want_cT       ) then
            params_array(ind) = P%csT2
            ind = ind +1
        end if
#endif

#ifdef COSMICFISH_EFTCAMB
        ! EFTCAMB parameters
        ! w_DE parameters:
        if ( P%EFTflag /= 0 .and. P%EFTwDE /= 0 ) then
            if ( P%EFTwDE == 1 ) then ! DE with constant Eos
                params_array(ind) = P%EFTw0
                ind = ind +1
            else if ( P%EFTwDE == 2 ) then ! CPL parametrization
                params_array(ind)   = P%EFTw0
                params_array(ind+1) = P%EFTwa
                ind = ind +2
            else if ( P%EFTwDE == 3 ) then ! JBP parametrization
                params_array(ind)   = P%EFTw0
                params_array(ind+1) = P%EFTwa
                params_array(ind+2) = P%EFTwn
                ind = ind +3
            else if ( P%EFTwDE == 4 ) then ! turning point parametrization
                params_array(ind)   = P%EFTw0
                params_array(ind+1) = P%EFTwa
                params_array(ind+2) = P%EFTwat
                ind = ind +3
            else if ( P%EFTwDE == 5 ) then ! Taylor expansion
                params_array(ind)   = P%EFTw0
                params_array(ind+1) = P%EFTwa
                params_array(ind+2) = P%EFTw2
                params_array(ind+3) = P%EFTw3
                ind = ind +4
            else if ( P%EFTwDE == 6 ) then ! User defined
                stop 'You have to define the parametrization first'
            end if
        end if

        ! EFT functions parameters:
        if ( P%EFTflag == 1 ) then ! Pure EFT models
            ! EFTOmega:
            if ( P%PureEFTmodelOmega == 1 ) then
                params_array(ind) = P%EFTOmega0
                ind = ind +1
            else if ( P%PureEFTmodelOmega == 2 ) then
                params_array(ind) = P%EFTOmega0
                ind = ind +1
            else if ( P%PureEFTmodelOmega == 3 ) then
                params_array(ind)   = P%EFTOmega0
                params_array(ind+1) = P%EFTOmegaExp
                ind = ind +2
            else if ( P%PureEFTmodelOmega == 4 ) then
                params_array(ind)   = P%EFTOmega0
                params_array(ind+1) = P%EFTOmegaExp
                ind = ind +2
            else if ( P%PureEFTmodelOmega == 5 ) then
                stop 'You have to define the parametrization first'
            end if
            ! EFTGamma1
            if ( P%PureEFTmodelGamma1 == 1 ) then
                params_array(ind) = P%EFTGamma10
                ind = ind +1
            else if ( P%PureEFTmodelGamma1 == 2 ) then
                params_array(ind) = P%EFTGamma10
                ind = ind +1
            else if ( P%PureEFTmodelGamma1 == 3 ) then
                params_array(ind)   = P%EFTGamma10
                params_array(ind+1) = P%EFTGamma1Exp
                ind = ind +2
            else if ( P%PureEFTmodelGamma1 == 4 ) then
                params_array(ind)   = P%EFTGamma10
                params_array(ind+1) = P%EFTGamma1Exp
                ind = ind +2
            else if ( P%PureEFTmodelGamma1 == 5 ) then
                stop 'You have to define the parametrization first'
            end if
            ! EFTGamma2
            if ( P%PureEFTmodelGamma2 == 1 ) then
                params_array(ind) = P%EFTGamma20
                ind = ind +1
            else if ( P%PureEFTmodelGamma2 == 2 ) then
                params_array(ind) = P%EFTGamma20
                ind = ind +1
            else if ( P%PureEFTmodelGamma2 == 3 ) then
                params_array(ind)   = P%EFTGamma20
                params_array(ind+1) = P%EFTGamma2Exp
                ind = ind +2
            else if ( P%PureEFTmodelGamma2 == 4 ) then
                params_array(ind)   = P%EFTGamma20
                params_array(ind+1) = P%EFTGamma2Exp
                ind = ind +2
            else if ( P%PureEFTmodelGamma2 == 5 ) then
                stop 'You have to define the parametrization first'
            end if
            ! EFTGamma3
            if ( P%PureEFTmodelGamma3 == 1 ) then
                params_array(ind) = P%EFTGamma30
                ind = ind +1
            else if ( P%PureEFTmodelGamma3 == 2 ) then
                params_array(ind) = P%EFTGamma30
                ind = ind +1
            else if ( P%PureEFTmodelGamma3 == 3 ) then
                params_array(ind)   = P%EFTGamma30
                params_array(ind+1) = P%EFTGamma3Exp
                ind = ind +2
            else if ( P%PureEFTmodelGamma3 == 4 ) then
                params_array(ind)   = P%EFTGamma30
                params_array(ind+1) = P%EFTGamma3Exp
                ind = ind +2
            else if ( P%PureEFTmodelGamma3 == 5 ) then
                stop 'You have to define the parametrization first'
            end if
            ! EFTGamma4
            if ( P%PureEFTmodelGamma4 == 1 ) then
                params_array(ind) = P%EFTGamma40
                ind = ind +1
            else if ( P%PureEFTmodelGamma4 == 2 ) then
                params_array(ind) = P%EFTGamma40
                ind = ind +1
            else if ( P%PureEFTmodelGamma4 == 3 ) then
                params_array(ind)   = P%EFTGamma40
                params_array(ind+1) = P%EFTGamma4Exp
                ind = ind +2
            else if ( P%PureEFTmodelGamma4 == 4 ) then
                params_array(ind)   = P%EFTGamma40
                params_array(ind+1) = P%EFTGamma4Exp
                ind = ind +2
            else if ( P%PureEFTmodelGamma4 == 5 ) then
                stop 'You have to define the parametrization first'
            end if
            ! EFTGamma5
            if ( P%PureEFTmodelGamma5 == 1 ) then
                params_array(ind) = P%EFTGamma50
                ind = ind +1
            else if ( P%PureEFTmodelGamma5 == 2 ) then
                params_array(ind) = P%EFTGamma50
                ind = ind +1
            else if ( P%PureEFTmodelGamma5 == 3 ) then
                params_array(ind)   = P%EFTGamma50
                params_array(ind+1) = P%EFTGamma5Exp
                ind = ind +2
            else if ( P%PureEFTmodelGamma5 == 4 ) then
                params_array(ind)   = P%EFTGamma50
                params_array(ind+1) = P%EFTGamma5Exp
                ind = ind +2
            else if ( P%PureEFTmodelGamma5 == 5 ) then
                stop 'You have to define the parametrization first'
            end if
            ! EFTGamma6
            if ( P%PureEFTmodelGamma6 == 1 ) then
                params_array(ind) = P%EFTGamma60
                ind = ind +1
            else if ( P%PureEFTmodelGamma6 == 2 ) then
                params_array(ind) = P%EFTGamma60
                ind = ind +1
            else if ( P%PureEFTmodelGamma6 == 3 ) then
                params_array(ind)   = P%EFTGamma60
                params_array(ind+1) = P%EFTGamma6Exp
                ind = ind +2
            else if ( P%PureEFTmodelGamma6 == 4 ) then
                params_array(ind)   = P%EFTGamma60
                params_array(ind+1) = P%EFTGamma6Exp
                ind = ind +2
            else if ( P%PureEFTmodelGamma6 == 5 ) then
                stop 'You have to define the parametrization first'
            end if
        else if ( P%EFTflag == 2 ) then ! designer EFT models
            if ( P%DesignerEFTmodel == 1 ) then
                params_array(ind) = P%EFTB0
                ind = ind +1
            else if ( P%DesignerEFTmodel == 2 ) then
            end if
        else if ( P%EFTflag == 3 ) then ! EFT alternative parametrizations
            if ( P%AltParEFTmodel == 1 ) then
                if ( P%RPHmassPmodel == 1 ) then
                    params_array(ind) = P%RPHmassP0
                    ind = ind +1
                else if ( P%RPHmassPmodel == 2 ) then
                    params_array(ind)   = P%RPHmassP0
                    params_array(ind+1) = P%RPHmassPexp
                    ind = ind +2
                end if
                if ( P%RPHkineticitymodel == 1 ) then
                    params_array(ind) = P%RPHkineticity0
                    ind = ind +1
                else if ( P%RPHkineticitymodel == 2 ) then
                    params_array(ind)   = P%RPHkineticity0
                    params_array(ind+1) = P%RPHkineticityexp
                    ind = ind +2
                end if
                if ( P%RPHbraidingmodel == 1 ) then
                    params_array(ind) = P%RPHbraiding0
                    ind = ind +1
                else if ( P%RPHbraidingmodel == 2 ) then
                    params_array(ind)   = P%RPHbraiding0
                    params_array(ind+1) = P%RPHbraidingexp
                    ind = ind +2
                end if
                if ( P%RPHtensormodel == 1 ) then
                    params_array(ind) = P%RPHtensor0
                    ind = ind +1
                else if ( P%RPHtensormodel == 2 ) then
                    params_array(ind)   = P%RPHtensor0
                    params_array(ind+1) = P%RPHtensorexp
                    ind = ind +2
                end if
            end if
        else if ( P%EFTflag == 4 ) then ! Full mapping EFT model
            if ( P%FullMappingEFTmodel == 1 ) then
                if ( P%HoravaSolarSystem .eqv. .True. ) then
                    params_array(ind)   = P%Horava_lambda
                    params_array(ind+1) = P%Horava_eta
                    ind = ind +2
                else
                    params_array(ind)   = P%Horava_lambda
                    params_array(ind+1) = P%Horava_eta
                    params_array(ind+2) = P%Horava_xi
                    ind = ind +3
                end if
            end if
        end if
#endif

#ifdef COSMICFISH_MGCAMB
        if ( P%MGC_model == 1 ) then
            params_array(ind)   = P%B1
            params_array(ind+1) = P%B2
            params_array(ind+2) = P%lambda1_2
            params_array(ind+3) = P%lambda2_2
            params_array(ind+4) = P%ss
            ind = ind +5
        else if ( P%MGC_model == 2 ) then
            params_array(ind)   = P%MGQfix
            params_array(ind+1) = P%MGRfix
            ind = ind +2
        else if ( P%MGC_model == 3 ) then
            params_array(ind)   = P%Qnot
            params_array(ind+1) = P%Rnot
            params_array(ind+2) = P%sss
            ind = ind +3
        else if ( P%MGC_model == 4 ) then
            params_array(ind)   = P%B0
            ind = ind +1
        else if ( P%MGC_model ==5 ) then
            params_array(ind)   = P%B0
            params_array(ind+1) = P%B1
            params_array(ind+2) = P%ss
            ind = ind +3
        else if ( P%MGC_model ==6 ) then
            params_array(ind)   = P%Linder_gamma
            ind = ind +1
        else if ( P%MGC_model == 7 ) then
            params_array(ind)   = P%beta_star
            params_array(ind+1) = P%xi_star
            params_array(ind+2) = P%a_star
            ind = ind +3
        else if ( P%MGC_model == 8 ) then
            params_array(ind)   = P%beta0
            params_array(ind+1) = P%xi0
            params_array(ind+2) = P%DilR
            params_array(ind+3) = P%DilS
            ind = ind +4
        else if ( P%MGC_model == 9 ) then
            params_array(ind)   = P%F_R0
            params_array(ind+1) = P%FRn
            ind = ind +2
        else if ( P%MGC_model ==10 ) then
            params_array(ind)   = P%beta0
            params_array(ind+1) = P%A_2
            ind = ind +2
        else if ( P%MGC_model ==11 ) then
            params_array(ind)   = P%E11_mg
            params_array(ind+1) = P%E22_mg
            ind = ind +2
            if ( FP%fisher_par%want_c1   ) then
               params_array(ind) = P%c1_mg
               ind = ind +1
            end if
            if ( FP%fisher_par%want_c2   ) then
               params_array(ind) = P%c2_mg
               ind = ind +1
            end if
            if ( FP%fisher_par%want_lambda   ) then
               params_array(ind) = P%lam_mg
               ind = ind +1
            end if
        else if ( P%MGC_model ==12 ) then
            params_array(ind)   = P%E11_mg
            params_array(ind+1) = P%E22_mg
            params_array(ind+2) = P%E12_mg
            params_array(ind+3) = P%E21_mg
            ind = ind +4
            if ( FP%fisher_par%want_c1   ) then
               params_array(ind) = P%c1_mg
               ind = ind +1
            end if
            if ( FP%fisher_par%want_c2   ) then
               params_array(ind) = P%c2_mg
               ind = ind +1
            end if
            if ( FP%fisher_par%want_lambda   ) then
               params_array(ind) = P%lam_mg
               ind = ind +1
            end if
        end if
#endif

        ! Quality check:
        if ( ind-1 /= num_param ) stop 'CAMB_params_to_params_array: Wrong assignment of parameters'

    end subroutine CAMB_params_to_params_array

    ! ---------------------------------------------------------------------------------------------
    !> This subroutine updates the values of the parameters P and FP based on the values contained
    !! in the parameter array.
    subroutine params_array_to_CAMB_param( P, FP, num_param, params_array )

        implicit none

        Type(CAMBparams)                             :: P            !< Output CAMBparams
        Type(cosmicfish_params)                       :: FP           !< Output Cosmicfish params
        integer                        , intent(in)  :: num_param    !< Input length of the parameter array
        real(dl), dimension(num_param) , intent(in)  :: params_array !< Input parameter array

        integer  :: ind, ind2, i, j, k
        real(dl) :: temp

        ind = 1
        i   = 1

        if ( FP%fisher_par%want_hubble ) then
            if ( FP%fisher_par%want_ombh2  ) i = i + 1
            if ( FP%fisher_par%want_omch2  ) i = i + 1
            if ( FP%fisher_par%want_omnuh2 ) i = i + 1
            P%h0 = params_array(i)*100._dl
        end if

        if ( FP%fisher_par%want_ombh2  ) then
            P%omegab = params_array(ind)/(P%h0/100._dl)**2
            ind = ind + 1
        end if

        if ( FP%fisher_par%want_omch2  ) then
            P%omegac = params_array(ind)/(P%h0/100._dl)**2
            ind = ind + 1
        end if

        if ( FP%fisher_par%want_omnuh2  ) then
            P%omegan = params_array(ind)/(P%h0/100._dl)**2
            ind = ind + 1
        end if

        if ( FP%fisher_par%want_hubble ) then
            ind = ind + 1
        end if

        P%omegav = 1- P%omegab- P%omegac - P%omegan

        ! helium abundance
        if ( FP%fisher_par%want_helium_fraction ) then
            P%yhe = params_array(ind)
            ind = ind + 1
        end if
        ! massless neutrinos
        if ( FP%fisher_par%want_massless ) then
            P%Num_Nu_massless = params_array(ind)
            ind = ind + 1
        end if
        ! primordial spectrum parameters:
        ! scalar_amp
        if ( FP%fisher_par%want_scalar_amp            ) then
            P%InitPower%ScalarPowerAmp(1) = exp(params_array(ind))/(10._dl**10)
            ind = ind + 1
        end if
        ! scalar_spectral_index
        if ( FP%fisher_par%want_scalar_spectral_index ) then
            P%InitPower%an(1) = params_array(ind)
            ind = ind + 1
        end if
        ! scalar_nrun
        if ( FP%fisher_par%want_scalar_nrun           ) then
            P%InitPower%n_run(1) = params_array(ind)
            ind = ind + 1
        end if
        ! tensor_spectral_index
        if ( FP%fisher_par%want_tensor_spectral_index ) then
            P%InitPower%ant(1) = params_array(ind)
            ind = ind + 1
        end if
        ! initial_ratio
        if ( FP%fisher_par%want_initial_ratio         ) then
            P%InitPower%rat(1) = params_array(ind)
            ind = ind + 1
        end if
        ! reionization: reionization optical depth
        if ( FP%fisher_par%want_re_optical_depth      ) then
            P%Reion%optical_depth = params_array(ind)
            ind = ind + 1
        end if

        ! Nuisance parameters for LSS:
        ! bias in every counts redshift windows:
        if ( FP%fisher_par%want_bias ) then
            ind2 = 0
            do i = 1, num_redshiftwindows
                if ( Redshift_w(i)%kind == window_counts ) then
                    ind2 = ind2 +1
                    FP%fisher_cls%bias(ind2) = params_array(ind)
                    Redshift_w(i)%bias = FP%fisher_cls%bias(ind2)
                    ind = ind +1
                end if
            end do
        end if

        ! nuissance parameters for supernovae:
        ! alpha_SN
        if ( FP%fisher_par%want_alpha_SN ) then
            FP%fisher_SN%alpha_SN = params_array(ind)
            ind = ind +1
        end if
        ! beta_SN
        if ( FP%fisher_par%want_beta_SN  ) then
            FP%fisher_SN%beta_SN = params_array(ind)
            ind = ind +1
        end if
        ! M0_SN
        if ( FP%fisher_par%want_M0_SN    ) then
            FP%fisher_SN%M0_SN = params_array(ind)
            ind = ind +1
        end if

#ifdef COSMICFISH_CAMB
        ! ppf parameters:
        ! w0_ppf
        if ( FP%fisher_par%want_w0_ppf   ) then
            P%w_lam = params_array(ind)
            ind = ind +1
        end if
        ! wa_ppf
        if ( FP%fisher_par%want_wa_ppf   ) then
            P%wa_ppf = params_array(ind)
            ind = ind +1
        end if
        ! cs_ppf
        if ( FP%fisher_par%want_cs_ppf   ) then
            P%cs2_lam = params_array(ind)
            call setcgammappf( P )
            ind = ind +1
        end if
        ! cT
        if ( FP%fisher_par%want_cT       ) then
            P%csT2 = params_array(ind)
            ind = ind +1
        end if
#endif

#ifdef COSMICFISH_EFTCAMB
        ! EFTCAMB parameters
        ! w_DE parameters:
        if ( P%EFTflag /= 0 .and. P%EFTwDE /= 0 ) then
            if ( P%EFTwDE == 1 ) then ! DE with constant Eos
                P%EFTw0 = params_array(ind)
                ind = ind +1
            else if ( P%EFTwDE == 2 ) then ! CPL parametrization
                P%EFTw0 = params_array(ind)
                P%EFTwa = params_array(ind+1)
                ind = ind +2
            else if ( P%EFTwDE == 3 ) then ! JBP parametrization
                P%EFTw0 = params_array(ind)
                P%EFTwa = params_array(ind+1)
                P%EFTwn = params_array(ind+2)
                ind = ind +3
            else if ( P%EFTwDE == 4 ) then ! turning point parametrization
                P%EFTw0 = params_array(ind)
                P%EFTwa = params_array(ind+1)
                P%EFTwat = params_array(ind+2)
                ind = ind +3
            else if ( P%EFTwDE == 5 ) then ! Taylor expansion
                P%EFTw0 = params_array(ind)
                P%EFTwa = params_array(ind+1)
                P%EFTw2 = params_array(ind+2)
                P%EFTw3 = params_array(ind+3)
                ind = ind +4
            else if ( P%EFTwDE == 6 ) then ! User defined
                stop 'You have to define the parametrization first'
            end if
        end if

        ! EFT functions parameters:
        if ( P%EFTflag == 1 ) then ! Pure EFT models
            ! EFTOmega:
            if ( P%PureEFTmodelOmega == 1 ) then
                P%EFTOmega0 = params_array(ind)
                ind = ind +1
            else if ( P%PureEFTmodelOmega == 2 ) then
                P%EFTOmega0 = params_array(ind)
                ind = ind +1
            else if ( P%PureEFTmodelOmega == 3 ) then
                P%EFTOmega0   = params_array(ind)
                P%EFTOmegaExp = params_array(ind+1)
                ind = ind +2
            else if ( P%PureEFTmodelOmega == 4 ) then
                P%EFTOmega0   = params_array(ind)
                P%EFTOmegaExp = params_array(ind+1)
                ind = ind +2
            else if ( P%PureEFTmodelOmega == 5 ) then
                stop 'You have to define the parametrization first'
            end if
            ! EFTGamma1
            if ( P%PureEFTmodelGamma1 == 1 ) then
                P%EFTGamma10 = params_array(ind)
                ind = ind +1
            else if ( P%PureEFTmodelGamma1 == 2 ) then
                P%EFTGamma10 = params_array(ind)
                ind = ind +1
            else if ( P%PureEFTmodelGamma1 == 3 ) then
                P%EFTGamma10   = params_array(ind)
                P%EFTGamma1Exp = params_array(ind+1)
                ind = ind +2
            else if ( P%PureEFTmodelGamma1 == 4 ) then
                P%EFTGamma10   = params_array(ind)
                P%EFTGamma1Exp = params_array(ind+1)
                ind = ind +2
            else if ( P%PureEFTmodelGamma1 == 5 ) then
                stop 'You have to define the parametrization first'
            end if
            ! EFTGamma2
            if ( P%PureEFTmodelGamma2 == 1 ) then
                P%EFTGamma20 = params_array(ind)
                ind = ind +1
            else if ( P%PureEFTmodelGamma2 == 2 ) then
                P%EFTGamma20 = params_array(ind)
                ind = ind +1
            else if ( P%PureEFTmodelGamma2 == 3 ) then
                P%EFTGamma20   = params_array(ind)
                P%EFTGamma2Exp = params_array(ind+1)
                ind = ind +2
            else if ( P%PureEFTmodelGamma2 == 4 ) then
                P%EFTGamma20   = params_array(ind)
                P%EFTGamma2Exp = params_array(ind+1)
                ind = ind +2
            else if ( P%PureEFTmodelGamma2 == 5 ) then
                stop 'You have to define the parametrization first'
            end if
            ! EFTGamma3
            if ( P%PureEFTmodelGamma3 == 1 ) then
                P%EFTGamma30 = params_array(ind)
                ind = ind +1
            else if ( P%PureEFTmodelGamma3 == 2 ) then
                P%EFTGamma30 = params_array(ind)
                ind = ind +1
            else if ( P%PureEFTmodelGamma3 == 3 ) then
                P%EFTGamma30   = params_array(ind)
                P%EFTGamma3Exp = params_array(ind+1)
                ind = ind +2
            else if ( P%PureEFTmodelGamma3 == 4 ) then
                P%EFTGamma30   = params_array(ind)
                P%EFTGamma3Exp = params_array(ind+1)
                ind = ind +2
            else if ( P%PureEFTmodelGamma3 == 5 ) then
                stop 'You have to define the parametrization first'
            end if
            ! EFTGamma4
            if ( P%PureEFTmodelGamma4 == 1 ) then
                P%EFTGamma40 = params_array(ind)
                ind = ind +1
            else if ( P%PureEFTmodelGamma4 == 2 ) then
                P%EFTGamma40 = params_array(ind)
                ind = ind +1
            else if ( P%PureEFTmodelGamma4 == 3 ) then
                P%EFTGamma40   = params_array(ind)
                P%EFTGamma4Exp = params_array(ind+1)
                ind = ind +2
            else if ( P%PureEFTmodelGamma4 == 4 ) then
                P%EFTGamma40   = params_array(ind)
                P%EFTGamma4Exp = params_array(ind+1)
                ind = ind +2
            else if ( P%PureEFTmodelGamma4 == 5 ) then
                stop 'You have to define the parametrization first'
            end if
            ! EFTGamma5
            if ( P%PureEFTmodelGamma5 == 1 ) then
                P%EFTGamma50 = params_array(ind)
                ind = ind +1
            else if ( P%PureEFTmodelGamma5 == 2 ) then
                P%EFTGamma50 = params_array(ind)
                ind = ind +1
            else if ( P%PureEFTmodelGamma5 == 3 ) then
                P%EFTGamma50   = params_array(ind)
                P%EFTGamma5Exp = params_array(ind+1)
                ind = ind +2
            else if ( P%PureEFTmodelGamma5 == 4 ) then
                P%EFTGamma50   = params_array(ind)
                P%EFTGamma5Exp = params_array(ind+1)
                ind = ind +2
            else if ( P%PureEFTmodelGamma5 == 5 ) then
                stop 'You have to define the parametrization first'
            end if
            ! EFTGamma6
            if ( P%PureEFTmodelGamma6 == 1 ) then
                P%EFTGamma60 = params_array(ind)
                ind = ind +1
            else if ( P%PureEFTmodelGamma6 == 2 ) then
                P%EFTGamma60 = params_array(ind)
                ind = ind +1
            else if ( P%PureEFTmodelGamma6 == 3 ) then
                P%EFTGamma60   = params_array(ind)
                P%EFTGamma6Exp = params_array(ind+1)
                ind = ind +2
            else if ( P%PureEFTmodelGamma6 == 4 ) then
                P%EFTGamma60   = params_array(ind)
                P%EFTGamma6Exp = params_array(ind+1)
                ind = ind +2
            else if ( P%PureEFTmodelGamma6 == 5 ) then
                stop 'You have to define the parametrization first'
            end if
        else if ( P%EFTflag == 2 ) then ! designer EFT models
            if ( P%DesignerEFTmodel == 1 ) then
                P%EFTB0 = params_array(ind)
                ind = ind +1
            else if ( P%DesignerEFTmodel == 2 ) then
            end if
        else if ( P%EFTflag == 3 ) then ! EFT alternative parametrizations
            if ( P%AltParEFTmodel == 1 ) then
                if ( P%RPHmassPmodel == 1 ) then
                    P%RPHmassP0 = params_array(ind)
                    ind = ind +1
                else if ( P%RPHmassPmodel == 2 ) then
                    P%RPHmassP0 = params_array(ind)
                    P%RPHmassPexp = params_array(ind+1)
                    ind = ind +2
                end if
                if ( P%RPHkineticitymodel == 1 ) then
                    P%RPHkineticity0 = params_array(ind)
                    ind = ind +1
                else if ( P%RPHkineticitymodel == 2 ) then
                    P%RPHkineticity0   = params_array(ind)
                    P%RPHkineticityexp = params_array(ind+1)
                    ind = ind +2
                end if
                if ( P%RPHbraidingmodel == 1 ) then
                    P%RPHbraiding0 = params_array(ind)
                    ind = ind +1
                else if ( P%RPHbraidingmodel == 2 ) then
                    P%RPHbraiding0   = params_array(ind)
                    P%RPHbraidingexp = params_array(ind+1)
                    ind = ind +2
                end if
                if ( P%RPHtensormodel == 1 ) then
                    P%RPHtensor0 = params_array(ind)
                    ind = ind +1
                else if ( P%RPHtensormodel == 2 ) then
                    P%RPHtensor0   = params_array(ind)
                    P%RPHtensorexp = params_array(ind+1)
                    ind = ind +2
                end if
            end if
        else if ( P%EFTflag == 4 ) then ! Full mapping EFT model
            if ( P%FullMappingEFTmodel == 1 ) then
                if ( P%HoravaSolarSystem .eqv. .True. ) then
                    P%Horava_lambda = params_array(ind)
                    P%Horava_eta    = params_array(ind+1)
                    ind = ind +2
                else
                    P%Horava_lambda = params_array(ind)
                    P%Horava_eta    = params_array(ind+1)
                    P%Horava_xi     = params_array(ind+2)
                    ind = ind +3
                end if
            end if
        end if
#endif

#ifdef COSMICFISH_MGCAMB
        if ( P%MGC_model == 1 ) then
            P%B1        = params_array(ind)
            P%B2        = params_array(ind+1)
            P%lambda1_2 = params_array(ind+2)
            P%lambda2_2 = params_array(ind+3)
            P%ss        = params_array(ind+4)
            ind = ind +5
        else if ( P%MGC_model == 2 ) then
            P%MGQfix      = params_array(ind)
            P%MGRfix      = params_array(ind+1)
            ind = ind +2
        else if ( P%MGC_model == 3 ) then
            P%Qnot = params_array(ind)
            P%Rnot = params_array(ind+1)
            P%sss  = params_array(ind+2)
            ind = ind +3
        else if ( P%MGC_model == 4 ) then
            P%B0   = params_array(ind)
            P%B1         = 4.d0/3.d0
            P%lambda1_2  = (P%B0*(299792458.d-3)**2)/(2.d0*p%H0**2)
            P%B2         = 0.5d0
            P%lambda2_2  = P%B1*P%lambda1_2
            P%ss         = 4.d0
            ind = ind +1
        else if ( P%MGC_model ==5 ) then
            P%B0         = params_array(ind)
            P%B1         = params_array(ind+1)
            P%lambda1_2  = (P%B0*(299792458.d-3)**2)/(2.d0*p%H0**2)
            P%B2         = 2.d0/P%B1 -1.d0
            P%lambda2_2  = P%B1*P%lambda1_2
            P%ss         = params_array(ind+2)
            ind = ind +3
        else if ( P%MGC_model ==6 ) then
            P%Linder_gamma = params_array(ind)
            ind = ind +1
        else if ( P%MGC_model == 7 ) then
            P%beta_star = params_array(ind)
            P%xi_star   = params_array(ind+1)
            P%a_star    = params_array(ind+2)
            P%GRtrans   = P%a_star
            ind = ind +3
        else if ( P%MGC_model == 8 ) then
            P%beta0 = params_array(ind)
            P%xi0   = params_array(ind+1)
            P%DilR  = params_array(ind+2)
            P%DilS  = params_array(ind+3)
            ind = ind +4
        else if ( P%MGC_model == 9 ) then
            P%F_R0  = params_array(ind)
            P%FRn   = params_array(ind+1)
            P%beta0 = 1.d0/sqrt(6.d0)
            ind = ind +2
        else if ( P%MGC_model ==10 ) then
            P%beta0 = params_array(ind)
            P%A_2   = params_array(ind+1)
            ind = ind +2
        else if ( P%MGC_model == 11 ) then
            P%E11_mg   = params_array(ind)
            P%E22_mg   = params_array(ind+1)
            ind = ind +2
            if ( FP%fisher_par%want_c1   ) then
               P%c1_mg = params_array(ind)
               ind = ind +1
            end if
            if ( FP%fisher_par%want_c2   ) then
               P%c2_mg = params_array(ind)
               ind = ind +1
            end if
            if ( FP%fisher_par%want_lambda   ) then
               P%lam_mg = params_array(ind)
               ind = ind +1
            end if
        else if ( P%MGC_model == 12 ) then
            P%E11_mg   = params_array(ind)
            P%E22_mg   = params_array(ind+1)
            P%E12_mg   = params_array(ind+2)
            P%E21_mg   = params_array(ind+3)
            ind = ind +4
            if ( FP%fisher_par%want_c1   ) then
               P%c1_mg = params_array(ind)
               ind = ind +1
            end if
            if ( FP%fisher_par%want_c2   ) then
               P%c2_mg = params_array(ind)
               ind = ind +1
            end if
            if ( FP%fisher_par%want_lambda   ) then
               P%lam_mg = params_array(ind)
               ind = ind +1
            end if
        end if
#endif

        ! Quality check:
        if ( ind-1 /= num_param ) stop 'Wrong assignment of parameters'
        if (.not. CAMB_ValidateParams(P)) stop 'Stopped due to parameter error'

    end subroutine params_array_to_CAMB_param

    ! ---------------------------------------------------------------------------------------------
    !> This subroutine saves the Fisher matrix to file taking care of putting into the file all the
    !> relevant informations.
    subroutine save_Fisher_to_file( P, FP, fish_dim, fisher_matrix, filename )

        implicit none

        Type(CAMBparams)                          :: P                        !< Input CAMBparams
        Type(cosmicfish_params)                    :: FP                       !< Input Cosmicfish params
        integer , intent(in)                      :: fish_dim                 !< Dimension of the input Fisher matrix
        real(dl), intent(in), dimension(fish_dim, fish_dim) :: fisher_matrix  !< Input Fisher matrix
        character(len=*), intent(in)              :: filename                 !< Name of the file where the Fisher matrix will be saved

        integer :: ind
        character(LEN=500) :: param_name
        character(LEN=500) :: param_name_latex
        real(dl), allocatable, dimension(:) :: parameters

        allocate( parameters(fish_dim) )
        call CAMB_params_to_params_array( P, FP, fish_dim, parameters )

        open(unit=666, FILE=filename, ACTION="write", STATUS="replace")

        ! write the header:
        write(666,'(a)') '#'
        write(666,'(a)') '# This file contains a Fisher matrix created with the CosmicFish code.'
        ! write the parameter names:
        write(666,'(a)') '#'
        write(666,'(a)') '# The parameters of this Fisher matrix are:'
        write(666,'(a)') '#'

        do ind = 1, fish_dim
            call Fisher_param_names(P, FP, ind, param_name, param_name_latex)
            write(666,'(a,I5,a,a,a,a,a,E25.16)') '#', ind, '    ', trim(param_name), '    ', trim(param_name_latex), '    ', parameters(ind)
        end do

        write(666,'(a)') '#'
        write(666,'(a)')

        ! write the Fisher matrix:
        do ind = 1, fish_dim
            write(666,'(2x,*(E25.16,2x))') fisher_matrix(ind, :)
        end do

        close(666)

    end subroutine save_Fisher_to_file

    ! ---------------------------------------------------------------------------------------------
    !> This subroutine saves a file containing the parameter names for the Fisher matrix run.
    subroutine save_paramnames_to_file( P, FP, fish_dim, filename )

        implicit none

        Type(CAMBparams)                          :: P                        !< Input CAMBparams
        Type(cosmicfish_params)                    :: FP                       !< Input Cosmicfish params
        integer , intent(in)                      :: fish_dim                 !< Dimension of the parameter space
        character(len=*), intent(in)              :: filename                 !< Name of the file where the Fisher matrix paramnames file will be saved

        integer :: ind
        character(LEN=500) :: param_name
        character(LEN=500) :: param_name_latex
        real(dl), allocatable, dimension(:) :: parameters

        allocate( parameters(fish_dim) )
        call CAMB_params_to_params_array( P, FP, fish_dim, parameters )

        open(unit=666, FILE=filename, ACTION="write", STATUS="replace")

        ! write the header:

        write(666,'(a)') '#'
        write(666,'(a)') '# This file contains the parameter names for a Fisher matrix.'
        write(666,'(a)') '#'
        ! write the parameter names:

        do ind = 1, fish_dim
            call Fisher_param_names(P, FP, ind, param_name, param_name_latex)
            write(666,'(a,a,a,a,E25.16)') trim(param_name), '    ', trim(param_name_latex), '    ', parameters(ind)
        end do

        close(666)

    end subroutine save_paramnames_to_file

end module cosmicfish_camb_interface
