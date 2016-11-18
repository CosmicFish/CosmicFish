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

!> @file 004_cosmicfish_types.f90
!! This file contains the definitions of derived data types used to store the informations
!! about the fisher calculation and experimental specifications

!----------------------------------------------------------------------------------------
!> This module contains the definitions of derived data types used to store the informations
!! about the fisher calculation and experimental specifications

!> @author Marco Raveri

module cosmicfish_types

    use precision

    implicit none

    !----------------------------------------------------------------------------------------
    !> This derived data type contains all the informations about the base parameters to run.
    type cosmicfish_param_fisher

        ! base parameters selection flags:
        logical   :: want_ombh2                  !< Decide wether to include ombh2 in the Fisher parameters.
        logical   :: want_omch2                  !< Decide wether to include omch2 in the Fisher parameters.
        logical   :: want_omnuh2                 !< Decide wether to include omnuh2 in the Fisher parameters.
        logical   :: want_hubble                 !< Decide wether to include hubble or OmegaV in the Fisher parameters.
        logical   :: want_helium_fraction        !< Decide wether to include helium_fraction in the Fisher parameters.
        logical   :: want_massless               !< Decide wether to include massless neutrinos in the Fisher parameters.
        ! primordial spectrum selection flags:
        logical   :: want_scalar_amp             !< Decide wether to include scalar amplitude in the Fisher parameters.
        logical   :: want_scalar_spectral_index  !< Decide wether to include scalar spectral index in the Fisher parameters.
        logical   :: want_scalar_nrun            !< Decide wether to include scalar running of the spectral index in the Fisher parameters.
        logical   :: want_tensor_spectral_index  !< Decide wether to include tensor spectral index in the Fisher parameters.
        logical   :: want_initial_ratio          !< Decide wether to include r in the Fisher parameters.
        ! reionization:
        logical   :: want_re_optical_depth       !< Decide wether to include reionization optical depth in the Fisher parameters.
        ! LSS nuissance parameters::
        logical   :: want_bias                   !< Decide wether to include bias in the Fisher parameters.
        ! SN nuissance parameters:
        logical   :: want_alpha_SN               !< Decide wether to include alpha_SN in the Fisher parameters.
        logical   :: want_beta_SN                !< Decide wether to include beta_SN in the Fisher parameters.
        logical   :: want_M0_SN                  !< Decide wether to include M0_SN in the Fisher parameters.
        ! ppf parameters (only with camb):
#ifdef COSMICFISH_CAMB
        logical   :: want_w0_ppf                 !< Decide wether to include w0_ppf in the Fisher parameters.
        logical   :: want_wa_ppf                 !< Decide wether to include wa_ppf in the Fisher parameters.
        logical   :: want_cs_ppf                 !< Decide wether to include cs_ppf in the Fisher parameters.
#endif
#ifdef COSMICFISH_MGCAMB
        logical   :: want_c1                     !< Decide wether to include c1 in the Fisher parameters.
        logical   :: want_c2                     !< Decide wether to include c2 in the Fisher parameters.
        logical   :: want_lambda                 !< Decide wether to include lambda in the Fisher parameters.
#endif

    end type cosmicfish_param_fisher

    !----------------------------------------------------------------------------------------
    !> This derived data type contains all the informations that the cosmicfish code needs to
    !! compute the Fisher matrix for Cls.
    type cosmicfish_fisher_cls

        ! selection flags:
        logical   :: Fisher_want_CMB_T       !< Decide wether to include CMB temperature in the Cls matrix.
        logical   :: Fisher_want_CMB_E       !< Decide wether to include CMB E mode polarization in the Cls matrix.
        logical   :: Fisher_want_CMB_B       !< Decide wether to include CMB B mode polarization in the Cls matrix.
        logical   :: Fisher_want_CMB_lensing !< Decide wether to include CMB lensing in the Cls matrix.
        logical   :: Fisher_want_LSS_lensing !< Decide wether to include LSS lensing windows in the Cls matrix.
        logical   :: Fisher_want_LSS_counts  !< Decide wether to include LSS number counts windows in the Cls matrix.
        logical   :: Fisher_want_XC          !< Decide wether to include cross correlations between the windows.

        ! CMB noise specifications:
        integer   :: CMB_n_channels          !< Number of frequency channels of the CMB experiment
        real(dl)  :: CMB_TT_fsky             !< Sky fraction for CMB temperature
        real(dl)  :: CMB_EE_fsky             !< Sky fraction for CMB E mode polarization
        real(dl)  :: CMB_BB_fsky             !< Sky fraction for CMB B mode polarization

        integer   :: l_max_TT                !< Maximum multipole considered for CMB temperature
        integer   :: l_max_EE                !< Maximum multipole considered for CMB E mode polarization
        integer   :: l_max_BB                !< Maximum multipole considered for CMB B mode polarization

        integer   :: l_min_TT                !< Minimum multipole considered for CMB temperature
        integer   :: l_min_EE                !< Minimum multipole considered for CMB E mode polarization
        integer   :: l_min_BB                !< Minimum multipole considered for CMB B mode polarization

        real(dl), dimension(:), allocatable :: CMB_temp_sens !< Temperature sensitivity in each frequency channel
        real(dl), dimension(:), allocatable :: CMB_pol_sens  !< Polarization sensitivity in each frequency channel
        real(dl), dimension(:), allocatable :: CMB_fwhm      !< Beam fwhm in each frequency channel

        ! LSS windows specifications:
        integer   :: LSS_number_windows          !< Number of LSS windows, either lensing or number counts

        real(dl), dimension(:), allocatable :: LSS_num_galaxies          !< Number of galaxies in each window
        real(dl), dimension(:), allocatable :: LSS_intrinsic_ellipticity !< Intrinsic ellipticity in  each window, relevant only for lensing
        real(dl), dimension(:), allocatable :: LSS_fsky                  !< Sky fraction in each window
        integer , dimension(:), allocatable :: LSS_lmax                  !< Maximum multipole in each window
        integer , dimension(:), allocatable :: LSS_lmin                  !< Minimum multipole in each window

        ! other window specifications:
        integer  :: window_type                  !< Specifies the window type:
                                                 !! 1) Gaussian window (default)
                                                 !! 2) Euclid-like window
        ! window parameters:
        real(dl) :: window_alpha                 !< Parameter of the window function:
                                                 !! - In gaussian window does nothing;
                                                 !! - In Euclid-like window it specifies the low-z shape of the window function;
        real(dl) :: window_beta                  !< Parameter of the window function:
                                                 !! - In gaussian window does nothing;
                                                 !! - In Euclid-like window it specifies the hi-z shape of the window function;
        real(dl) :: redshift_zero                !< Parameter of the window function:
                                                 !! - In gaussian window does nothing;
                                                 !! - In Euclid-like window it specifies the median redshift of the window;
        real(dl) :: photoz_error                 !< Parameter for the photo-z error:
                                                 !! - In gaussian window does nothing;
                                                 !! - In Euclid-like window it specifies the photo-z error;

        real(dl), dimension(:), allocatable :: bias     !< Fiducial bias parameters. Nuissance parameters for GC

    end type cosmicfish_fisher_cls

    !----------------------------------------------------------------------------------------
    !> This derived data type contains all the informations that the cosmicfish code needs to
    !! compute the Fisher matrix for supernovae.
    type cosmicfish_fisher_SN

        ! SN windows:
        integer  :: number_SN_windows                            !< Number of SN survey bins in redshift
        integer , dimension(:), allocatable :: SN_number         !< Number of SN in each bin
        real(dl), dimension(:), allocatable :: SN_redshift_start !< Bin starting redshift
        real(dl), dimension(:), allocatable :: SN_redshift_end   !< Bin ending redshift
        integer  :: total_SN_number                              !< Total number of supernovae

        ! SN Fisher Monte Carlo average:
        integer  :: SN_Fisher_MC_samples !< number of realizations over which the SN Fisher is computed

        ! fiducial SN parameters:
        real(dl) :: alpha_SN  !< fiducial value of the SN alpha
        real(dl) :: beta_SN   !< fiducial value of the SN beta
        real(dl) :: M0_SN     !< fiducial value of the SN M0

        ! parameters to generate the SN mock data:
        real(dl) :: color_dispersion   !< Dispersion of the color of the mock SN cathalog
        real(dl) :: stretch_dispersion !< Dispersion of the stretch of the mock SN cathalog

        ! parameters to generate the mock SN covariance:
        real(dl) :: magnitude_sigma  !< error on the absolute magnitude
        real(dl) :: c_sigmaz         !< error in the redshift determination
        real(dl) :: sigma_lens_0     !< lensing error

        real(dl) :: dcolor_offset  !< error in the color at redshift zero
        real(dl) :: dcolor_zcorr   !< coefficient for the redshift dependence of the color error
        real(dl) :: dshape_offset  !< error in the stretch at redshift zero
        real(dl) :: dshape_zcorr   !< coefficient for the redshift dependence of the stretch error

        real(dl) :: cov_ms_offset  !< covariance between the error in magnitude and stretch at redshift zero
        real(dl) :: cov_ms_zcorr   !< coefficient for the redshift dependence of the covariance between the error in magnitude and stretch
        real(dl) :: cov_mc_offset  !< covariance between the error in magnitude and color at redshift zero
        real(dl) :: cov_mc_zcorr   !< coefficient for the redshift dependence of the covariance between the error in magnitude and color
        real(dl) :: cov_sc_offset  !< covariance between the error in stretch and color at redshift zero
        real(dl) :: cov_sc_zcorr   !< coefficient for the redshift dependence of the covariance between the error in stretch and color

    end type cosmicfish_fisher_SN

    !----------------------------------------------------------------------------------------
    !> This derived data type contains all the informations that the cosmicfish code needs to
    !! compute the Fisher matrix for redshift drift.
    type cosmicfish_fisher_RD

        !RD specifications
        integer :: number_RD_redshifts
        integer, dimension(:), allocatable  :: RD_number         !< Number of observations in each bin
        real(dl), dimension(:), allocatable :: RD_redshift       !< center of redshift bins
        integer                             :: exptype           !< type of survey (will be 1=E-ELT, 2=SKA, 3=SKA second derivative? 4=CHIME?)
        real(dl)                            :: obs_time          !< time interval for RD observations
        real(dl)                            :: signoise          !< signal to noise of redshift observations


    end type cosmicfish_fisher_RD

    !----------------------------------------------------------------------------------------
    !> This derived data type contains all the informations that the cosmicfish code needs to
    !! compute the Fisher matrix for derived parameters.
    type cosmicfish_fisher_derived

        ! base derived parameters selection flags:
        logical   :: want_omegab                !< Decide wether to include omegab in the Fisher parameters.
        logical   :: want_omegac                !< Decide wether to include omegac in the Fisher parameters.
        logical   :: want_omegan                !< Decide wether to include omegan in the Fisher parameters.
        logical   :: want_omegav                !< Decide wether to include omegav in the Fisher parameters.
        logical   :: want_omegak                !< Decide wether to include omegak in the Fisher parameters.
        logical   :: want_omegam                !< Decide wether to include omegam in the Fisher parameters.
        logical   :: want_theta                 !< Decide wether to include theta in the Fisher parameters.
        logical   :: want_mnu                   !< Decide wether to include mnu in the Fisher parameters.
        logical   :: want_zre                   !< Decide wether to include zre in the Fisher parameters.
        logical   :: want_neff                  !< Decide wether to include Neff in the Fisher parameters.

        ! derived parameters in addition to the standard derived:
        logical   :: want_sigma8                !< Decide wether to include sigma8(z) (tomography) in the Fisher parameters.
        logical   :: want_loghubble             !< Decide wether to include H(z) (tomography) in the Fisher parameters.
        logical   :: want_logDA                 !< Decide wether to include DA(z) (tomography) in the Fisher parameters.

        ! derived parameters tomography:
        integer                             :: FD_num_redshift  !< Number of tomographic redshifts
        real(dl), dimension(:), allocatable :: FD_redshift      !< The values of the redshifts

#ifdef COSMICFISH_MGCAMB
        !derived parameters for planck parametrization
        logical   :: want_mu0                   !< Decide wether to include mu0 in the Fisher parameters.
        logical   :: want_gam0                  !< Decide wether to include gamma0 in the Fisher parameters.
#endif

    end type cosmicfish_fisher_derived

    !----------------------------------------------------------------------------------------
    !> This derived data type contains all the informations that the cosmicfish code needs to
    !! run. In particular it contains several derived types containing the informations about
    !! the experimental noise for different experiments.
    type cosmicfish_params

        ! General cosmicfish flags:
        logical   :: adaptivity             !< Decide wether all the derivatives are computed with a
                                            !! fixed stepsize or with an adaptive algorithm.

        integer   :: cosmicfish_feedback     !< Level of feedback from the cosmicfish code
                                            !! 0 = no feedback
                                            !! 1 = basic feedback (default)
                                            !! 2 = extended feedback
                                            !! 3 = debug feedback

        logical   :: cosmicfish_want_cls     !< Decide wether to compute the Fisher matrix for cls.
        logical   :: cosmicfish_want_SN      !< Decide wether to compute the Fisher matrix for supernovae.
        logical   :: cosmicfish_want_Mpk     !< Decide wether to compute the Fisher matrix for the matter power spectrum.
        logical   :: cosmicfish_want_RD      !< Decide wether to compute the Fisher matrix for redshift drift.
        logical   :: cosmicfish_want_derived !< Decide wether to compute derived parameters.

        ! other parameters:
        real(dl)  :: output_factor          !< Overall factor for the Cls.

        ! Specifications for the Cls Fisher matrix:
        type(cosmicfish_param_fisher)   :: fisher_par  !< contains what is needed to choose the parameters to run
        type(cosmicfish_fisher_cls)     :: fisher_cls  !< contains what is needed for Cls Fisher matrix
        type(cosmicfish_fisher_SN)      :: fisher_SN   !< contains what is needed for SN Fisher matrix
        type(cosmicfish_fisher_RD)      :: fisher_RD   !< contains what is needed for RD Fisher matrix
        type(cosmicfish_fisher_derived) :: fisher_der !< contains what is needed for derived Fisher matrix

        ! Output root:
        character(len=:), allocatable :: outroot    !< Root for all the filenames that will be produced by CosmicFish

    end type cosmicfish_params

end module cosmicfish_types
