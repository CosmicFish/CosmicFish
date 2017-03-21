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

    close(100)


    ! all the successful test should return the following.
    stop 100

end program test_params_io
