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

!> @file 005_init_from_file.f90
!! This file contains the subroutines that are needed to initialize camb and cosmicfish
!! with a parameter file.

!----------------------------------------------------------------------------------------
!> This module contains the subroutine and functions to initialize camb and cosmicfish
!! from a parameter file.

!> @author Marco Raveri

module init_from_file

    use precision
    use cosmicfish_types
    use IniFile
    use ModelParams

    implicit none

    private

    public init_cosmicfish_from_file

contains

    ! ---------------------------------------------------------------------------------------------
    !> This subroutine initializes camb parameters from a parameter file.
    subroutine init_cosmicfish_from_file( P, FP, filename, param_out_name )

        use constants, only : f_21cm, COBE_CMBTemp
        use CAMB
        use CAMBmain, only : Alens
        use Bispectrum

        implicit none

        Type(CAMBparams)                       :: P                !< CAMBparams object to be filled with the parameters
        Type(cosmicfish_params)                :: FP               !< CosmicFish parameter object
        character(len=*), intent(in)           :: filename         !< name of the file frow which parameters are read
        character(len=*), intent(in), optional :: param_out_name   !< name of the file to which to write the parameters

        logical :: bad
        logical :: DoCounts = .false.
        integer :: i, ind
        Type (TRedWin), pointer :: RedWin
        character(LEN=Ini_max_string_len) numstr, S, outroot, version_check
        real(dl) :: nmassive
        character(LEN=:), allocatable :: ParamDir

        !MMmod
        character(LEN=Ini_max_string_len) red_ind

        ! open the parameter file:
        call Ini_Open(filename, 1, bad, .false.)
        if (bad) stop 'Error opening parameter file'
        Ini_fail_on_not_found = .false.
        ! output root:
        outroot = Ini_Read_String('output_root')
        if (outroot /= '') outroot = trim(outroot) // '_'
        ! set CAMB parameters to defaul at the beginning:
        call CAMB_SetDefParams(P)
        ! high l template:
        if (Ini_HasKey('highL_unlensed_cl_template')) then
            highL_unlensed_cl_template=  Ini_Read_String_Default('highL_unlensed_cl_template',highL_unlensed_cl_template)
        else
            ! assume that the file highL_unlensed_cl_template is in the parameter folder:
            i = scan( filename, '/' , .true. )
            ParamDir = ''
            if ( i .gt. 1 ) then
                ParamDir = filename(1:i)
            endif
            highL_unlensed_cl_template = concat(ParamDir, highL_unlensed_cl_template)
        end if
        ! read what to do:
        P%WantScalars = Ini_Read_Logical('get_scalar_cls',.true.)
        P%WantVectors = Ini_Read_Logical('get_vector_cls',.false.)
        P%WantTensors = Ini_Read_Logical('get_tensor_cls',.false.)

        P%Want_CMB =  Ini_Read_Logical('want_CMB',.true.)
        P%Want_CMB_lensing =  Ini_Read_Logical('want_CMB_lensing',.true.)
        P%Want_CMB_lensing =  P%Want_CMB .or. P%Want_CMB_lensing

        P%WantCls= P%WantScalars .or. P%WantTensors .or. P%WantVectors
        P%PK_WantTransfer=Ini_Read_Logical('get_transfer')
        ! read the window functions:
        if (P%WantScalars) then
            num_redshiftwindows = Ini_Read_Int('num_redshiftwindows',0)
        else
            num_redshiftwindows = 0
        end if
        limber_windows = Ini_Read_Logical('limber_windows',limber_windows)
        if (limber_windows) limber_phiphi = Ini_Read_Int('limber_phiphi',limber_phiphi)
        if (num_redshiftwindows>0) then
            DoRedshiftLensing = Ini_Read_Logical('DoRedshiftLensing',.false.)
            Kmax_Boost = Ini_Read_Double('Kmax_Boost',Kmax_Boost)
        end if
        Do21cm = Ini_Read_Logical('Do21cm', .false.)
        num_extra_redshiftwindows = 0

        ! initialize the window functions:
        do i=1, num_redshiftwindows
            RedWin => Redshift_w(i)
            call InitRedshiftWindow(RedWin)
            write (numstr,*) i
            numstr=adjustl(numstr)
            RedWin%Redshift = Ini_Read_Double('redshift('//trim(numstr)//')')
            S = Ini_Read_String('redshift_kind('//trim(numstr)//')')
            if (S=='21cm') then
                RedWin%kind = window_21cm
            elseif (S=='counts') then
                RedWin%kind = window_counts
            elseif (S=='lensing') then
                RedWin%kind = window_lensing
            else
                write (*,*) i, 'Error: unknown type of window '//trim(S)
                stop
            end if
            RedWin%a = 1/(1+RedWin%Redshift)
            if (RedWin%kind /= window_21cm) then
                RedWin%sigma = Ini_Read_Double('redshift_sigma('//trim(numstr)//')')
                RedWin%sigma_z = RedWin%sigma
            else
                Do21cm = .true.
                RedWin%sigma = Ini_Read_Double('redshift_sigma_Mhz('//trim(numstr)//')')
                if (RedWin%sigma < 0.003) then
                    write(*,*) 'WARNING:Window very narrow.'
                    write(*,*) ' --> use transfer functions and transfer_21cm_cl =T ?'
                end if
                !with 21cm widths are in Mhz, make dimensionless scale factor
                RedWin%sigma = RedWin%sigma/(f_21cm/1e6)
                RedWin%sigma_z = RedWin%sigma*(1+RedWin%RedShift)**2
                write(*,*) i,'delta_z = ', RedWin%sigma_z
            end if
            if (RedWin%kind == window_counts) then
                DoCounts = .true.
                RedWin%bias = Ini_Read_Double('redshift_bias('//trim(numstr)//')')
                RedWin%dlog10Ndm = Ini_Read_Double('redshift_dlog10Ndm('//trim(numstr)//')',0.d0)
                if (DoRedshiftLensing) then
                    num_extra_redshiftwindows=num_extra_redshiftwindows+1
                    RedWin%mag_index = num_extra_redshiftwindows
                end if
            end if
        end do

        if (Do21cm) then
            print*, 'Not yet implemented'
            stop
            line_basic = Ini_Read_Logical('line_basic')
            line_distortions = Ini_read_Logical('line_distortions')
            line_extra = Ini_Read_Logical('line_extra')

            line_phot_dipole = Ini_read_Logical('line_phot_dipole')
            line_phot_quadrupole = Ini_Read_Logical('line_phot_quadrupole')
            line_reionization = Ini_Read_Logical('line_reionization')

            use_mK = Ini_read_Logical('use_mK')
            if (DebugMsgs) then
                write (*,*) 'Doing 21cm'
                write (*,*) 'dipole = ',line_phot_dipole, ' quadrupole =', line_phot_quadrupole
            end if
        else
            line_extra = .false.
        end if

        if (DoCounts) then
            counts_density = Ini_read_Logical('counts_density')
            counts_redshift = Ini_read_Logical('counts_redshift')
            counts_radial = Ini_read_Logical('counts_radial')
            counts_evolve = Ini_read_Logical('counts_evolve')
            counts_timedelay = Ini_read_Logical('counts_timedelay')
            counts_ISW = Ini_read_Logical('counts_ISW')
            counts_potential = Ini_read_Logical('counts_potential')
            counts_velocity = Ini_read_Logical('counts_velocity')
        end if

        P%OutputNormalization=outNone

        AccuracyBoost  = Ini_Read_Double('accuracy_boost',AccuracyBoost)
        lAccuracyBoost = Ini_Read_Real('l_accuracy_boost',lAccuracyBoost)
        HighAccuracyDefault = Ini_Read_Logical('high_accuracy_default',HighAccuracyDefault)

        P%NonLinear = Ini_Read_Int('do_nonlinear',NonLinear_none)

        evolve_delta_xe = Ini_read_Logical('evolve_delta_xe', .false.)

        P%DoLensing = .false.
        if (P%WantCls) then
            if (P%WantScalars  .or. P%WantVectors) then
                P%Max_l = Ini_Read_Int('l_max_scalar')
                P%Max_eta_k = Ini_Read_Double('k_eta_max_scalar',P%Max_l*2._dl)
                if (P%WantScalars) then
                    P%DoLensing = Ini_Read_Logical('do_lensing',.false.)
                    if (P%DoLensing) lensing_method = Ini_Read_Int('lensing_method',1)
                end if
                if (P%WantVectors) then
                    if (P%WantScalars .or. P%WantTensors) stop 'Must generate vector modes on their own'
                    i = Ini_Read_Int('vector_mode')
                    if (i==0) then
                        vec_sig0 = 1
                        Magnetic = 0
                    else if (i==1) then
                        Magnetic = -1
                        vec_sig0 = 0
                    else
                        stop 'vector_mode must be 0 (regular) or 1 (magnetic)'
                    end if
                end if
            end if

            if (P%WantTensors) then
                P%Max_l_tensor = Ini_Read_Int('l_max_tensor')
                P%Max_eta_k_tensor =  Ini_Read_Double('k_eta_max_tensor',Max(500._dl,P%Max_l_tensor*2._dl))
            end if
        endif

        P%window_kmax_boost = Ini_Read_Double('window_kmax_boost', 1._dl)

#ifdef COSMICFISH_CAMB
        !  Read cosmological parameters.
        call DarkEnergy_ReadParams( P, DefIni )
#endif

        P%h0     = Ini_Read_Double('hubble')

        if (Ini_Read_Logical('use_physical',.false.)) then
            P%omegab = Ini_Read_Double('ombh2')/(P%H0/100)**2
            P%omegac = Ini_Read_Double('omch2')/(P%H0/100)**2
            P%omegan = Ini_Read_Double('omnuh2')/(P%H0/100)**2
            P%omegav = 1- Ini_Read_Double('omk') - P%omegab-P%omegac - P%omegan
        else
            P%omegab = Ini_Read_Double('omega_baryon')
            P%omegac = Ini_Read_Double('omega_cdm')
            P%omegav = Ini_Read_Double('omega_lambda')
            P%omegan = Ini_Read_Double('omega_neutrino')
        end if

        P%tcmb   = Ini_Read_Double('temp_cmb',COBE_CMBTemp)

#ifdef COSMICFISH_EFTCAMB
        !  Read EFTCAMB parameters.

        ! 1) Initialization of EFTCAMB flags.

        P%EFTflag = Ini_Read_Int('EFTflag',0)

        P%PureEFTmodelOmega  = Ini_Read_Int('PureEFTmodelOmega',0)
        P%PureEFTmodelGamma1 = Ini_Read_Int('PureEFTmodelGamma1',0)
        P%PureEFTmodelGamma2 = Ini_Read_Int('PureEFTmodelGamma2',0)
        P%PureEFTmodelGamma3 = Ini_Read_Int('PureEFTmodelGamma3',0)
        P%PureEFTmodelGamma4 = Ini_Read_Int('PureEFTmodelGamma4',0)
        P%PureEFTmodelGamma5 = Ini_Read_Int('PureEFTmodelGamma5',0)
        P%PureEFTmodelGamma6 = Ini_Read_Int('PureEFTmodelGamma6',0)

        P%DesignerEFTmodel = Ini_Read_Int('DesignerEFTmodel',1)
        P%AltParEFTmodel   = Ini_Read_Int('AltParEFTmodel',1)
        P%FullMappingEFTmodel = Ini_Read_Int('FullMappingEFTmodel',1)

        ! 2) Initialization of EFTCAMB model properties flags.

        ! read the DE eos model selection flag:
        P%EFTwDE = Ini_Read_Int('EFTwDE',0)
        ! read pure EFT Horndeski model selection flag:
        P%PureEFTHorndeski = Ini_Read_Logical('PureEFTHorndeski',.false.)
        ! read RPH model selection flags:
        P%RPHmassPmodel      = Ini_Read_Int('RPHmassPmodel',0)
        P%RPHkineticitymodel = Ini_Read_Int('RPHkineticitymodel',0)
        P%RPHbraidingmodel   = Ini_Read_Int('RPHbraidingmodel',0)
        P%RPHtensormodel     = Ini_Read_Int('RPHtensormodel',0)
        ! read the Horava Solar System Free flag:
        P%HoravaSolarSystem  = Ini_Read_Logical('HoravaSolarSystem',.false.)

        ! 3) Initialization of EFTCAMB stability flags:

        P%EFT_mathematical_stability = Ini_Read_Logical('EFT_mathematical_stability',.true.)
        P%EFT_physical_stability     = Ini_Read_Logical('EFT_physical_stability',.true.)
        P%EFTAdditionalPriors        = Ini_Read_Logical('EFTAdditionalPriors',.true.)
        P%MinkowskyPriors            = Ini_Read_Logical('MinkowskyPriors',.true.)

        ! 4) Initialization of EFTCAMB model parameters.

        ! read the DE eos parameters:
        P%EFTw0  = Ini_Read_Double('EFTw0',-1._dl)
        P%EFTwa  = Ini_Read_Double('EFTwa',0._dl)
        P%EFTwn  = Ini_Read_Double('EFTwn',2._dl)
        P%EFTwat = Ini_Read_Double('EFTwat',1._dl)
        P%EFtw2  = Ini_Read_Double('EFtw2',0._dl)
        P%EFTw3  = Ini_Read_Double('EFTw3',0._dl)
        ! read pure EFT parameters:
        P%EFTOmega0    = Ini_Read_Double('EFTOmega0', 0.0_dl)
        P%EFTOmegaExp  = Ini_Read_Double('EFTOmegaExp', 0.0_dl)
        P%EFTGamma10   = Ini_Read_Double('EFTGamma10', 0.0_dl)
        P%EFTGamma1Exp = Ini_Read_Double('EFTGamma1Exp', 0.0_dl)
        P%EFTGamma20   = Ini_Read_Double('EFTGamma20', 0.0_dl)
        P%EFTGamma2Exp = Ini_Read_Double('EFTGamma2Exp', 0.0_dl)
        P%EFTGamma30   = Ini_Read_Double('EFTGamma30', 0.0_dl)
        P%EFTGamma3Exp = Ini_Read_Double('EFTGamma3Exp', 0.0_dl)
        P%EFTGamma40   = Ini_Read_Double('EFTGamma40', 0.0_dl)
        P%EFTGamma4Exp = Ini_Read_Double('EFTGamma4Exp', 0.0_dl)
        P%EFTGamma50   = Ini_Read_Double('EFTGamma50', 0.0_dl)
        P%EFTGamma5Exp = Ini_Read_Double('EFTGamma5Exp', 0.0_dl)
        P%EFTGamma60   = Ini_Read_Double('EFTGamma60', 0.0_dl)
        P%EFTGamma6Exp = Ini_Read_Double('EFTGamma6Exp', 0.0_dl)
        ! read f(R) parameters:
        P%EFTB0 = Ini_Read_Double('EFTB0', 0.0_dl)
        ! read RPH parameters:
        P%RPHmassP0        = Ini_Read_Double('RPHmassP0', 0.0_dl)
        P%RPHmassPexp      = Ini_Read_Double('RPHmassPexp', 0.0_dl)
        P%RPHkineticity0   = Ini_Read_Double('RPHkineticity0', 0.0_dl)
        P%RPHkineticityexp = Ini_Read_Double('RPHkineticityexp', 0.0_dl)
        P%RPHbraiding0     = Ini_Read_Double('RPHbraiding0', 0.0_dl)
        P%RPHbraidingexp   = Ini_Read_Double('RPHbraidingexp', 0.0_dl)
        P%RPHtensor0       = Ini_Read_Double('RPHtensor0', 0.0_dl)
        P%RPHtensorexp     = Ini_Read_Double('RPHtensorexp', 0.0_dl)
        ! read Horava parameters:
        P%Horava_xi      = Ini_Read_Double('Horava_xi', 0.0_dl)
        P%Horava_lambda  = Ini_Read_Double('Horava_lambda', 0.0_dl)
        P%Horava_eta     = Ini_Read_Double('Horava_eta', 0.0_dl)
#endif

#ifdef COSMICFISH_MGCAMB
        P%MGC_model = Ini_Read_Int('model',0)
        P%GRtrans   = Ini_Read_Double('GRtrans',0.d0)

        if ( P%MGC_model == 1 ) then
            P%B1        = Ini_Read_Double('B1',0.d0)
            P%B2        = Ini_Read_Double('B2',0.d0)
            P%lambda1_2 = Ini_Read_Double('lambda1_2',0.d0)
            P%lambda2_2 = Ini_Read_Double('lambda2_2',0.d0)
            P%ss        = Ini_Read_Double('ss',0.d0)
        else if ( P%MGC_model == 2 ) then
            P%MGQfix      = Ini_Read_Double('MGQfix',1.d0)
            P%MGRfix      = Ini_Read_Double('MGRfix',1.d0)
        else if ( P%MGC_model == 3 ) then
            P%Qnot = Ini_Read_Double('Qnot',1.d0)
            P%Rnot = Ini_Read_Double('Rnot',1.d0)
            P%sss  = Ini_Read_Double('sss',0.d0)
        else if ( P%MGC_model == 4 ) then
            P%B0         = Ini_Read_Double('B0',0.d0)
            P%B1         = 4.d0/3.d0
            P%lambda1_2  = (P%B0*(299792458.d-3)**2)/(2.d0*p%H0**2)
            P%B2         = 0.5d0
            P%lambda2_2  = P%B1*P%lambda1_2
            P%ss         = 4.d0
        else if ( P%MGC_model ==5 ) then
            P%B0         = Ini_Read_Double('B0',0.d0)
            P%B1         = Ini_Read_Double('beta1',0.d0)
            P%lambda1_2  = (P%B0*(299792458.d-3)**2)/(2.d0*p%H0**2)
            P%B2         = 2.d0/P%B1 -1.d0
            P%lambda2_2  = P%B1*P%lambda1_2
            P%ss         = Ini_Read_Double('s',0.d0)
        else if ( P%MGC_model ==6 ) then
            P%Linder_gamma = Ini_Read_Double('Linder_gamma',0.d0)
        else if ( P%MGC_model == 7 ) then
            P%beta_star = Ini_Read_Double('beta_star', 0.d0)
            P%xi_star   = Ini_Read_Double('xi_star', 0.d0)
            P%a_star    = Ini_Read_Double('a_star', 0.d0)
            P%GRtrans   = P%a_star
        else if ( P%MGC_model == 8 ) then
            P%beta0 = Ini_Read_Double('beta0', 0.d0)
            P%xi0   = Ini_Read_Double('xi0', 0.d0)
            P%DilR  = Ini_Read_Double('DilR', 0.d0)
            P%DilS  = Ini_Read_Double('DilS', 0.d0)
        else if ( P%MGC_model == 9 ) then
            P%F_R0  = Ini_Read_Double('F_R0', 0.d0)
            P%FRn   = Ini_Read_Double('FRn', 0.d0)
            P%beta0 = 1.d0/sqrt(6.d0)
        else if ( P%MGC_model ==10 ) then
            P%beta0 = Ini_Read_Double('beta0', 0.d0)
            P%A_2   = Ini_Read_Double('A2',0.d0)
        else if ( P%MGC_model ==11 ) then
            P%E11_mg= Ini_Read_Double('E11', 0.d0)
            P%E22_mg= Ini_Read_Double('E22', 0.d0)
            P%c1_mg = Ini_Read_Double('c1', 0.d0)
            P%c2_mg = Ini_Read_Double('c2', 0.d0)
            P%lam_mg= Ini_Read_Double('lambda', 0.d0)
        else if ( P%MGC_model ==12 ) then
            P%E11_mg= Ini_Read_Double('E11', 0.d0)
            P%E22_mg= Ini_Read_Double('E22', 0.d0)
            P%E12_mg= Ini_Read_Double('E12', 0.d0)
            P%E21_mg= Ini_Read_Double('E21', 0.d0)
            P%c1_mg = Ini_Read_Double('c1', 0.d0)
            P%c2_mg = Ini_Read_Double('c2', 0.d0)
            P%lam_mg= Ini_Read_Double('lambda', 0.d0)
        else if ( P%MGC_model /= 0 ) then
            print*, '***MGCAMB: please choose a model***'
            stop
        end if
#endif

        P%yhe    = Ini_Read_Double('helium_fraction',0.24_dl)
        P%Num_Nu_massless  = Ini_Read_Double('massless_neutrinos')

        ! massive neutrinos:
        P%Nu_mass_eigenstates = Ini_Read_Int('nu_mass_eigenstates',1)
        if (P%Nu_mass_eigenstates > max_nu) stop 'too many mass eigenstates'

        numstr = Ini_Read_String('massive_neutrinos')
        read(numstr, *) nmassive
        if (abs(nmassive-nint(nmassive))>1e-6) stop 'massive_neutrinos should now be integer (or integer array)'
        read(numstr,*) P%Nu_Mass_numbers(1:P%Nu_mass_eigenstates)

        P%Num_Nu_massive = sum(P%Nu_Mass_numbers(1:P%Nu_mass_eigenstates))

        if (P%Num_Nu_massive>0) then
            P%share_delta_neff = Ini_Read_Logical('share_delta_neff', .true.)
            numstr = Ini_Read_String('nu_mass_degeneracies')
            if (P%share_delta_neff) then
                if (numstr/='') write (*,*) 'WARNING: nu_mass_degeneracies ignored when share_delta_neff'
            else
                if (numstr=='') stop 'must give degeneracies for each eigenstate if share_delta_neff=F'
                read(numstr,*) P%Nu_mass_degeneracies(1:P%Nu_mass_eigenstates)
            end if
            numstr = Ini_Read_String('nu_mass_fractions')
            if (numstr=='') then
                if (P%Nu_mass_eigenstates >1) stop 'must give nu_mass_fractions for the eigenstates'
                P%Nu_mass_fractions(1)=1
            else
                read(numstr,*) P%Nu_mass_fractions(1:P%Nu_mass_eigenstates)
            end if
        end if

        !JD 08/13 begin changes for nonlinear lensing of CMB + LSS compatibility
        !P%Transfer%redshifts -> P%Transfer%PK_redshifts and P%Transfer%num_redshifts -> P%Transfer%PK_num_redshifts
        !in the P%WantTransfer loop.
        if (((P%NonLinear==NonLinear_lens .or. P%NonLinear==NonLinear_both) .and. (P%DoLensing .or. num_redshiftwindows>0)) &
            .or. P%PK_WantTransfer) then
            P%Transfer%high_precision=  Ini_Read_Logical('transfer_high_precision',.false.)
        else
            P%transfer%high_precision = .false.
        endif

        if (P%PK_WantTransfer)  then
            P%WantTransfer  = .true.
            P%transfer%kmax          =  Ini_Read_Double('transfer_kmax')
            P%transfer%k_per_logint  =  Ini_Read_Int('transfer_k_per_logint')
            P%transfer%PK_num_redshifts =  Ini_Read_Int('transfer_num_redshifts')

            if (Do21cm) transfer_21cm_cl = Ini_Read_Logical('transfer_21cm_cl',.false.)
            if (transfer_21cm_cl .and. P%transfer%kmax > 800) then
                !Actually line widths are important at significantly larger scales too
                write (*,*) 'WARNING: kmax very large. '
                write(*,*) ' -- Neglected line width effects will dominate'
            end if

            transfer_interp_matterpower = Ini_Read_Logical('transfer_interp_matterpower ', transfer_interp_matterpower)
            transfer_power_var = Ini_read_int('transfer_power_var',transfer_power_var)
            if (P%transfer%PK_num_redshifts > max_transfer_redshifts) stop 'Too many redshifts'
            do i=1, P%transfer%PK_num_redshifts
                P%transfer%PK_redshifts(i)  = Ini_Read_Double_Array('transfer_redshift',i,0._dl)
            end do
        else
            P%Transfer%PK_num_redshifts = 1
            P%Transfer%PK_redshifts = 0
        end if

        if ((P%NonLinear==NonLinear_lens .or. P%NonLinear==NonLinear_both) .and. &
            (P%DoLensing .or. num_redshiftwindows > 0)) then
            P%WantTransfer  = .true.
            call Transfer_SetForNonlinearLensing(P%Transfer)
        end if

        call Transfer_SortAndIndexRedshifts(P%Transfer)
        !JD 08/13 end changes

        P%transfer%kmax=P%transfer%kmax*(P%h0/100._dl)

        Ini_fail_on_not_found = .false.

        DebugParam = Ini_Read_Double('DebugParam',DebugParam)
        ALens = Ini_Read_Double('Alens',Alens)

        call Reionization_ReadParams(P%Reion, DefIni)
        call InitialPower_ReadParams(P%InitPower, DefIni, P%WantTensors)
        call Recombination_ReadParams(P%Recomb, DefIni)
        if (Ini_HasKey('recombination')) then
            i = Ini_Read_Int('recombination',1)
            if (i/=1) stop 'recombination option deprecated'
        end if

        call Bispectrum_ReadParams(BispectrumParams, DefIni, outroot)

        if (P%WantScalars .or. P%WantTransfer) then
            P%Scalar_initial_condition = Ini_Read_Int('initial_condition',initial_adiabatic)
            if (P%Scalar_initial_condition == initial_vector) then
                P%InitialConditionVector=0
                numstr = Ini_Read_String('initial_vector',.true.)
                read (numstr,*) P%InitialConditionVector(1:initial_iso_neutrino_vel)
            end if
            if (P%Scalar_initial_condition/= initial_adiabatic) use_spline_template = .false.
        end if

        if (P%WantScalars) then
            has_cl_2D_array = .true.
        end if

        Ini_fail_on_not_found = .false.

        !optional parameters controlling the computation
        P%AccuratePolarization = Ini_Read_Logical('accurate_polarization',.true.)
        P%AccurateReionization = Ini_Read_Logical('accurate_reionization',.false.)
        P%AccurateBB = Ini_Read_Logical('accurate_BB',.false.)
        P%DerivedParameters = Ini_Read_Logical('derived_parameters',.true.)

        version_check = Ini_Read_String('version_check')
        if (version_check == '') then
            !tag the output used parameters .ini file with the version of CAMB being used now
            call TNameValueList_Add(DefIni%ReadValues, 'version_check', version)
        else if (version_check /= version) then
            write(*,*) 'WARNING: version_check does not match this CAMB version'
        end if
        !Mess here to fix typo with backwards compatibility
        if (Ini_HasKey('do_late_rad_trunction')) then
            DoLateRadTruncation = Ini_Read_Logical('do_late_rad_trunction',.true.)
            if (Ini_HasKey('do_late_rad_truncation')) stop 'check do_late_rad_xxxx'
        else
            DoLateRadTruncation = Ini_Read_Logical('do_late_rad_truncation',.true.)
        end if
        DoTensorNeutrinos = Ini_Read_Logical('do_tensor_neutrinos',DoTensorNeutrinos )
        FeedbackLevel = Ini_Read_Int('feedback_level',FeedbackLevel)

        P%MassiveNuMethod  = Ini_Read_Int('massive_nu_approx',Nu_best)

        ThreadNum      = Ini_Read_Int('number_of_threads',ThreadNum)
        use_spline_template = Ini_Read_Logical('use_spline_template',use_spline_template)

        DoTensorNeutrinos = DoTensorNeutrinos .or. HighAccuracyDefault
        if (do_bispectrum) then
            lSampleBoost   = 50
        else
            lSampleBoost   = Ini_Read_Double('l_sample_boost',lSampleBoost)
        end if

        ! From now on read cosmicfish parameters:

        FP%outroot = TRIM(outroot)

        ! initialize general stuff:
        FP%adaptivity         = Ini_Read_Logical('adaptivity',.false.)
        FP%cosmicfish_feedback = Ini_Read_Int('cosmicfish_feedback',0)

        FP%cosmicfish_want_cls  = Ini_Read_Logical('cosmicfish_want_cls',.false.)
        FP%cosmicfish_want_SN   = Ini_Read_Logical('cosmicfish_want_SN' ,.false.)
        FP%cosmicfish_want_Mpk  = Ini_Read_Logical('cosmicfish_want_Mpk',.false.)
        FP%cosmicfish_want_RD   = Ini_Read_Logical('cosmicfish_want_RD' ,.false.)
        FP%cosmicfish_want_derived = Ini_Read_Logical('cosmicfish_want_derived',.false.)

        FP%output_factor = Ini_Read_Double('CMB_outputscale',1.d0)

        ! read parameter selection:
        FP%fisher_par%want_ombh2                   = Ini_Read_Logical( 'param[ombh2]'   ,.false. )
        FP%fisher_par%want_omch2                   = Ini_Read_Logical( 'param[omch2]'   ,.false. )
        FP%fisher_par%want_omnuh2                  = Ini_Read_Logical( 'param[omnuh2]'  ,.false. )
        FP%fisher_par%want_hubble                  = Ini_Read_Logical( 'param[hubble]'  ,.false. )
        FP%fisher_par%want_helium_fraction         = Ini_Read_Logical( 'param[Ye]'      ,.false. )
        FP%fisher_par%want_massless                = Ini_Read_Logical( 'param[nom_nu]'  ,.false. )
        FP%fisher_par%want_scalar_amp              = Ini_Read_Logical( 'param[As]'      ,.false. )
        FP%fisher_par%want_scalar_spectral_index   = Ini_Read_Logical( 'param[ns]'      ,.false. )
        FP%fisher_par%want_scalar_nrun             = Ini_Read_Logical( 'param[nsrun]'   ,.false. )
        FP%fisher_par%want_tensor_spectral_index   = Ini_Read_Logical( 'param[nt]'      ,.false. )
        FP%fisher_par%want_initial_ratio           = Ini_Read_Logical( 'param[r]'       ,.false. )
        FP%fisher_par%want_re_optical_depth        = Ini_Read_Logical( 'param[tau]'     ,.false. )
        FP%fisher_par%want_bias                    = Ini_Read_Logical( 'param[bias]'    ,.false. )
        FP%fisher_par%want_alpha_SN                = Ini_Read_Logical( 'param[alpha_SN]',.false. )
        FP%fisher_par%want_beta_SN                 = Ini_Read_Logical( 'param[beta_SN]' ,.false. )
        FP%fisher_par%want_M0_SN                   = Ini_Read_Logical( 'param[M0_SN]'   ,.false. )

#ifdef COSMICFISH_CAMB
        FP%fisher_par%want_w0_ppf                  = Ini_Read_Logical( 'param[w0_ppf]'   ,.false. )
        FP%fisher_par%want_wa_ppf                  = Ini_Read_Logical( 'param[wa_ppf]'   ,.false. )
        FP%fisher_par%want_cs_ppf                  = Ini_Read_Logical( 'param[cs_ppf]'   ,.false. )
        FP%fisher_par%want_cT                      = Ini_Read_Logical( 'param[cT]'       ,.false. )
#endif

#ifdef COSMICFISH_MGCAMB
        FP%fisher_par%want_c1                      = Ini_Read_Logical( 'param[c1]'       ,.false. )
        FP%fisher_par%want_c2                      = Ini_Read_Logical( 'param[c2]'       ,.false. )
        FP%fisher_par%want_lambda                  = Ini_Read_Logical( 'param[lambda]'   ,.false. )
#endif

        ! read Cls parameters:
        if ( FP%cosmicfish_want_cls ) then

            FP%fisher_cls%Fisher_want_CMB_T       = Ini_Read_Logical( 'Fisher_want_CMB_T'       , .False. )
            FP%fisher_cls%Fisher_want_CMB_E       = Ini_Read_Logical( 'Fisher_want_CMB_E'       , .False. )
            FP%fisher_cls%Fisher_want_CMB_B       = Ini_Read_Logical( 'Fisher_want_CMB_B'       , .False. )
            FP%fisher_cls%Fisher_want_CMB_lensing = Ini_Read_Logical( 'Fisher_want_CMB_lensing' , .False. )
            FP%fisher_cls%Fisher_want_LSS_lensing = Ini_Read_Logical( 'Fisher_want_LSS_lensing' , .False. )
            FP%fisher_cls%Fisher_want_LSS_counts  = Ini_Read_Logical( 'Fisher_want_LSS_counts'  , .False. )
            FP%fisher_cls%Fisher_want_XC          = Ini_Read_Logical( 'Fisher_want_XC'          , .False. )

            FP%fisher_cls%CMB_n_channels = Ini_Read_Int( 'CMB_n_channels' , 0 )
            FP%fisher_cls%CMB_TT_fsky    = Ini_Read_Double( 'CMB_TT_fsky' , 1._dl )
            FP%fisher_cls%CMB_EE_fsky    = Ini_Read_Double( 'CMB_EE_fsky' , 1._dl )
            FP%fisher_cls%CMB_BB_fsky    = Ini_Read_Double( 'CMB_BB_fsky' , 1._dl )

            FP%fisher_cls%l_max_TT = Ini_Read_Int( 'l_max_TT' , 0 )
            FP%fisher_cls%l_max_EE = Ini_Read_Int( 'l_max_EE' , 0 )
            FP%fisher_cls%l_max_BB = Ini_Read_Int( 'l_max_BB' , 0 )

            FP%fisher_cls%l_min_TT = Ini_Read_Int( 'l_min_TT' , 2 )
            FP%fisher_cls%l_min_EE = Ini_Read_Int( 'l_min_EE' , 2 )
            FP%fisher_cls%l_min_BB = Ini_Read_Int( 'l_min_BB' , 2 )

            allocate( FP%fisher_cls%CMB_temp_sens( FP%fisher_cls%CMB_n_channels ), &
                & FP%fisher_cls%CMB_pol_sens( FP%fisher_cls%CMB_n_channels ), &
                & FP%fisher_cls%CMB_fwhm( FP%fisher_cls%CMB_n_channels ) )

            do i=1, FP%fisher_cls%CMB_n_channels
                write (numstr,*) i
                numstr=adjustl(numstr)
                FP%fisher_cls%CMB_temp_sens(i) = Ini_Read_Double( 'CMB_temp_sens('//trim(numstr)//')', 0._dl )
                FP%fisher_cls%CMB_pol_sens(i)  = Ini_Read_Double( 'CMB_pol_sens('//trim(numstr)//')', 0._dl )
                FP%fisher_cls%CMB_fwhm(i)      = Ini_Read_Double( 'CMB_fwhm('//trim(numstr)//')', 0._dl )
            end do

            FP%fisher_cls%LSS_number_windows = Ini_Read_Int( 'LSS_number_windows' , 0 )

            allocate( FP%fisher_cls%LSS_num_galaxies( FP%fisher_cls%LSS_number_windows ), &
                & FP%fisher_cls%LSS_intrinsic_ellipticity( FP%fisher_cls%LSS_number_windows ), &
                & FP%fisher_cls%LSS_fsky( FP%fisher_cls%LSS_number_windows ), &
                & FP%fisher_cls%LSS_lmax( FP%fisher_cls%LSS_number_windows ), &
                & FP%fisher_cls%LSS_lmin( FP%fisher_cls%LSS_number_windows ) )

            do i=1, FP%fisher_cls%LSS_number_windows
                write (numstr,*) i
                numstr=adjustl(numstr)
                FP%fisher_cls%LSS_num_galaxies(i) = Ini_Read_Double('LSS_num_galaxies('//trim(numstr)//')', 0._dl )
                FP%fisher_cls%LSS_intrinsic_ellipticity(i) = Ini_Read_Double('LSS_intrinsic_ellipticity('//trim(numstr)//')', 0._dl )
                FP%fisher_cls%LSS_fsky(i)         = Ini_Read_Double('LSS_fsky('//trim(numstr)//')', 0._dl )
                FP%fisher_cls%LSS_lmax(i)         = Ini_Read_Int('LSS_lmax('//trim(numstr)//')', 0 )
                FP%fisher_cls%LSS_lmin(i)         = Ini_Read_Int('LSS_lmin('//trim(numstr)//')', 2 )
            end do

            FP%fisher_cls%window_type   = Ini_Read_Int( 'window_type' , 1 )
            FP%fisher_cls%window_alpha  = Ini_Read_Double('window_alpha' , 0._dl )
            FP%fisher_cls%window_beta   = Ini_Read_Double('window_beta'  , 0._dl )
            FP%fisher_cls%redshift_zero = Ini_Read_Double('redshift_zero', 0._dl )
            FP%fisher_cls%photoz_error  = Ini_Read_Double('photoz_error' , 0._dl )

            if ( FP%fisher_par%want_bias ) then
                ! number of GC windows:
                ind = 0
                do i = 1, FP%fisher_cls%LSS_number_windows
                    if ( Redshift_w(i)%kind == window_counts ) then
                        ind = ind +1
                    end if
                end do
                ! allocate the bias array:
                allocate( FP%fisher_cls%bias( ind ) )
                ! write the fiducial:
                ind = 0
                do i = 1, FP%fisher_cls%LSS_number_windows
                    if ( Redshift_w(i)%kind == window_counts ) then
                        ind = ind +1
                        FP%fisher_cls%bias( ind ) = Redshift_w(i)%bias
                    end if
                end do
            end if

        end if

        ! read SN parameters:

        ! fiducial SN parameters:
        FP%fisher_SN%alpha_SN  = Ini_Read_Double('alpha_SN' , 0._dl )
        FP%fisher_SN%beta_SN   = Ini_Read_Double('beta_SN' , 0._dl )
        FP%fisher_SN%M0_SN     = Ini_Read_Double('M0_SN' , 0._dl )

        if ( FP%cosmicfish_want_SN ) then

            FP%fisher_SN%number_SN_windows = Ini_Read_Int( 'number_SN_windows', 0 )

            allocate( FP%fisher_SN%SN_number(FP%fisher_SN%number_SN_windows), &
                & FP%fisher_SN%SN_redshift_start(FP%fisher_SN%number_SN_windows), &
                & FP%fisher_SN%SN_redshift_end(FP%fisher_SN%number_SN_windows) )

            do i=1, FP%fisher_SN%number_SN_windows
                write (numstr,*) i
                numstr=adjustl(numstr)
                FP%fisher_SN%SN_number(i)          = Ini_Read_Int( 'SN_number('//trim(numstr)//')', 0 )
                FP%fisher_SN%SN_redshift_start(i)  = Ini_Read_Double('SN_redshift_start('//trim(numstr)//')', 0._dl )
                FP%fisher_SN%SN_redshift_end(i)    = Ini_Read_Double('SN_redshift_end('//trim(numstr)//')', 0._dl )
            end do

            FP%fisher_SN%total_SN_number = sum( FP%fisher_SN%SN_number )

            ! SN Fisher Monte Carlo average:
            FP%fisher_SN%SN_Fisher_MC_samples = Ini_Read_Int( 'SN_Fisher_MC_samples', 1 )

            ! parameters to generate the SN mock data:
            FP%fisher_SN%color_dispersion   = Ini_Read_Double('color_dispersion' , 0._dl )
            FP%fisher_SN%stretch_dispersion = Ini_Read_Double('stretch_dispersion' , 0._dl )

            ! parameters to generate the mock SN covariance:
            FP%fisher_SN%magnitude_sigma  = Ini_Read_Double('magnitude_sigma' , 0._dl )
            FP%fisher_SN%c_sigmaz         = Ini_Read_Double('c_sigmaz' , 0._dl )
            FP%fisher_SN%sigma_lens_0     = Ini_Read_Double('sigma_lens_0' , 0._dl )

            FP%fisher_SN%dcolor_offset  = Ini_Read_Double('dcolor_offset' , 0._dl )
            FP%fisher_SN%dcolor_zcorr   = Ini_Read_Double('dcolor_zcorr' , 0._dl )
            FP%fisher_SN%dshape_offset  = Ini_Read_Double('dshape_offset' , 0._dl )
            FP%fisher_SN%dshape_zcorr   = Ini_Read_Double('dshape_zcorr' , 0._dl )

            FP%fisher_SN%cov_ms_offset  = Ini_Read_Double('cov_ms_offset' , 0._dl )
            FP%fisher_SN%cov_ms_zcorr   = Ini_Read_Double('cov_ms_zcorr' , 0._dl )
            FP%fisher_SN%cov_mc_offset  = Ini_Read_Double('cov_mc_offset' , 0._dl )
            FP%fisher_SN%cov_mc_zcorr   = Ini_Read_Double('cov_mc_zcorr' , 0._dl )
            FP%fisher_SN%cov_sc_offset  = Ini_Read_Double('cov_sc_offset' , 0._dl )
            FP%fisher_SN%cov_sc_zcorr   = Ini_Read_Double('cov_sc_zcorr' , 0._dl )

        end if

        ! read RD parameters:
        if ( FP%cosmicfish_want_RD ) then
            FP%fisher_RD%exptype = Ini_Read_Int( 'RD_exptype', 0 )
            if (FP%fisher_RD%exptype.eq.1) then
                FP%fisher_RD%number_RD_redshifts = Ini_Read_Int( 'number_RD_redshifts', 0 )
                allocate(FP%fisher_RD%RD_redshift(FP%fisher_RD%number_RD_redshifts),FP%fisher_RD%RD_number(FP%fisher_RD%number_RD_redshifts))
                do i=1,FP%fisher_RD%number_RD_redshifts
                    write(red_ind,*) i
                    FP%fisher_RD%RD_redshift(i) = Ini_Read_Double( 'RD_redshift('//trim(adjustl(red_ind))//')',0._dl)
                    FP%fisher_RD%RD_number(i)   = Ini_Read_Int( 'RD_source_number('//trim(adjustl(red_ind))//')', 0 ) !MM can be improved
                end do
                FP%fisher_RD%obs_time            = Ini_Read_Int( 'delta_time', 0 )
                FP%fisher_RD%signoise            = Ini_Read_Int( 'RD_sig_to_noise', 0 )
            else if (FP%fisher_RD%exptype.eq.2) then
                write(0,*)'not implemented yet'
                stop
            else if (FP%fisher_RD%exptype.eq.3) then
                write(0,*)'not implemented yet'
                stop
            end if
        end if

        ! read derived parameters:
        if ( FP%cosmicfish_want_derived )  then
            ! fixed parameters:
            FP%fisher_der%want_omegab = Ini_Read_Logical( 'param[omegab]'    ,.false. )
            FP%fisher_der%want_omegac = Ini_Read_Logical( 'param[omegac]'    ,.false. )
            FP%fisher_der%want_omegan = Ini_Read_Logical( 'param[omegan]'    ,.false. )
            FP%fisher_der%want_omegav = Ini_Read_Logical( 'param[omegav]'    ,.false. )
            FP%fisher_der%want_omegak = Ini_Read_Logical( 'param[omegak]'    ,.false. )
            FP%fisher_der%want_omegam = Ini_Read_Logical( 'param[omegam]'    ,.false. )
            FP%fisher_der%want_theta  = Ini_Read_Logical( 'param[theta]'     ,.false. )
            FP%fisher_der%want_mnu    = Ini_Read_Logical( 'param[mnu]'       ,.false. )
            FP%fisher_der%want_zre    = Ini_Read_Logical( 'param[zre]'       ,.false. )
            FP%fisher_der%want_neff   = Ini_Read_Logical( 'param[neff]'      ,.false. )
            ! tomographic parameters:
            FP%fisher_der%want_sigma8     = Ini_Read_Logical( 'param[sigma8]'    ,.false. )
            FP%fisher_der%want_loghubble  = Ini_Read_Logical( 'param[loghubble]' ,.false. )
            FP%fisher_der%want_logDA      = Ini_Read_Logical( 'param[logDA]'     ,.false. )
            FP%fisher_der%want_S8         = Ini_Read_Logical( 'param[S8]'        ,.false. )

            FP%fisher_der%FD_num_redshift = Ini_Read_Int( 'FD_num_redshift', 0 )

            allocate( FP%fisher_der%FD_redshift( FP%fisher_der%FD_num_redshift ) )

            do i=1, FP%fisher_der%FD_num_redshift
                write (numstr,*) i
                numstr=adjustl(numstr)

                FP%fisher_der%FD_redshift(i) = Ini_Read_Double( 'FD_redshift('//trim(numstr)//')', 0._dl )
            end do

#ifdef COSMICFISH_MGCAMB
            FP%fisher_der%want_mu0 = Ini_Read_Logical( 'param[mu0]'    ,.false. )
            FP%fisher_der%want_gam0 = Ini_Read_Logical( 'param[gamma0]'    ,.false. )
#endif

        end if

        ! read Mpk parameters:
        if ( FP%cosmicfish_want_Mpk ) then

        end if

        ! dump read parameters to file:
        if ( present(param_out_name) ) then
            call Ini_SaveReadValues(trim(param_out_name),1)
        else if ( FP%outroot /= '' ) then
            if (filename /= trim(FP%outroot) //'params.ini') then
                call Ini_SaveReadValues(trim(FP%outroot) //'params.ini',1)
            else
                write(*,*) 'WARNING: Output _params.ini not created as would overwrite input parameters.'
            end if
        end if

        ! finalize:
        call Ini_Close
        if (.not. CAMB_ValidateParams(P)) stop 'Stopped due to parameter error'
        if (global_error_flag/=0) stop 'Error in initialization.'

        ! check the consistency of CosmicFish parameters
        call check_parameters_consistency(P, FP)

    end subroutine init_cosmicfish_from_file

    ! ---------------------------------------------------------------------------------------------
    !> This subroutine checks the consistency of the cosmicfish parameters with the camb parameters
    subroutine check_parameters_consistency( P, FP )

        implicit none

        Type(CAMBparams)                       :: P                !< Input CAMBparams object
        Type(cosmicfish_params)                 :: FP               !< Input CosmicFish parameter object

        ! general check:
        if  ( .not. FP%cosmicfish_want_cls .and. &
            & .not. FP%cosmicfish_want_SN  .and. &
            & .not. FP%cosmicfish_want_RD  .and. &
            & .not. FP%cosmicfish_want_Mpk .and. &
            & .not. FP%cosmicfish_want_derived ) then
            write(*,*) 'WARNING: you want no Fisher!.'
        end if

        ! check the cls Fisher matrix:
        if ( FP%cosmicfish_want_cls ) then
            ! l_max check:
            if  ( P%Max_l+20 < FP%fisher_cls%l_max_TT .or. &
                & P%Max_l+20 < FP%fisher_cls%l_max_EE .or. &
                & P%Max_l+20 < FP%fisher_cls%l_max_BB ) then
                write(*,*) 'ERROR: l_max+20 is smaller that l_max_TT or l_max_EE or l_max_BB.'
                stop
            end if
        end if

        write(*,*) 'WARNING: check_parameters_consistency: not fully implemented.'

    end subroutine check_parameters_consistency

end module init_from_file
