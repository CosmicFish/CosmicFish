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

program test_params_io

    use precision
    use cosmicfish_types
    use init_from_file
    use ModelParams
    use AMLutils
    use IniFile

    implicit none

    Type(CAMBparams)       :: P
    Type(cosmicfish_params) :: FP
    character(LEN=Ini_max_string_len) paramfile

    paramfile = ''
    if (GetParamCount() /= 0)  paramfile = GetParam(1)
    if (paramfile == '') stop 'No parameter input file'

    call init_cosmicfish_from_file( P, FP, paramfile )

    open (unit = 100, file = "output")

    write(100,*) P%h0
    write(100,*) P%omegab
    write(100,*) P%omegac
    write(100,*) P%omegan
    write(100,*) P%omegav
    write(100,*) P%tcmb

#ifdef COSMICFISH_EFTCAMB

    write(100,*) P%EFTflag
    write(100,*) P%EFTwDE
    write(100,*) P%PureEFTmodelOmega
    write(100,*) P%PureEFTmodelGamma1

#endif

    write(100,*) FP%outroot
    write(100,*) FP%adaptivity
    write(100,*) FP%cosmicfish_feedback
    write(100,*) FP%cosmicfish_want_cls
    write(100,*) FP%cosmicfish_want_SN
    write(100,*) FP%cosmicfish_want_Mpk
    write(100,*) FP%output_factor

    write(100,*) FP%fisher_par%want_ombh2
    write(100,*) FP%fisher_par%want_omch2
    write(100,*) FP%fisher_par%want_omnuh2
    write(100,*) FP%fisher_par%want_hubble
    write(100,*) FP%fisher_par%want_helium_fraction
    write(100,*) FP%fisher_par%want_scalar_amp
    write(100,*) FP%fisher_par%want_scalar_spectral_index
    write(100,*) FP%fisher_par%want_scalar_nrun
    write(100,*) FP%fisher_par%want_tensor_spectral_index
    write(100,*) FP%fisher_par%want_initial_ratio
    write(100,*) FP%fisher_par%want_re_optical_depth
    write(100,*) FP%fisher_par%want_bias
    write(100,*) FP%fisher_par%want_alpha_SN
    write(100,*) FP%fisher_par%want_beta_SN
    write(100,*) FP%fisher_par%want_M0_SN

    close(100)

    ! all the successful test should return the following.
    stop 100

end program test_params_io
