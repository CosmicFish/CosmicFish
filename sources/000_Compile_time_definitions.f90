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

!> @file 000_Compile_time_definitions.f90
!! This file contains the CosmicFish compile time variables that are used by the code.

!----------------------------------------------------------------------------------------
!> This module contains the CosmicFish compile time variables that are used by the code.

!> @author Marco Raveri and Matteo Martinelli

module cosmicfish_def

    use precision

    implicit none

    ! CosmicFish version flag:
    character(LEN=*), parameter :: CosmicFish_version = 'V1.1 Mar17'

end module cosmicfish_def

!----------------------------------------------------------------------------------------
