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

!> @file 008_Fisher_calculator_Cls_windows.f90
!! This file contains the definition of the window functions that are used by CAMB sources.
!! If the considered experiement has a window that is not already implemennted here is the
!! place to code it.

!----------------------------------------------------------------------------------------
!> This module contains the definition of the window functions that are used by CAMB sources.

!> @author Marco Raveri

module Fisher_calculator_Cls_windows

    use precision
    use constants, only: const_pi
    use RedshiftSpaceData
    use cosmicfish_types

    implicit none

    type(cosmicfish_params), save :: Window_FP

    private

    public init_camb_sources_windows, gaussian_window, bin_window, flat_window, flat_smooth_window

contains

    ! ---------------------------------------------------------------------------------------------
    !> This function initialises the camb sources number counts and lensing window function.
    !! The default is the gaussian window that is included in camb sources.
    !! If another window is chosen it will be initialised by the CosmicFish value.
    subroutine init_camb_sources_windows( FP )

        implicit none

        Type(cosmicfish_params), intent(in)  :: FP !< Input Cosmicfish params

        Window_FP = FP

        if ( FP%fisher_cls%window_type == 1 ) then
            window_function_calc => default_window
        else if ( FP%fisher_cls%window_type == 2 ) then
            window_function_calc => bin_window
        else if ( FP%fisher_cls%window_type == 3 ) then
            window_function_calc => flat_window
        else if ( FP%fisher_cls%window_type == 4 ) then
            window_function_calc => flat_smooth_window
        else
            window_function_calc => Null()
        end if

    end subroutine init_camb_sources_windows

    ! ---------------------------------------------------------------------------------------------
    !> This function computes a gaussian window function.
    !! Note this is the total count distribution observed, not a fractional selection function
    !! on an underlying distribution.
    function gaussian_window(Win, z, winamp)

        implicit none

        Type(TRedWin)        :: Win      !< Window function object
        real(dl), intent(in) :: z        !< Input redshift
        real(dl)             :: winamp   !< Output amplitude of the window

        real(dl) :: gaussian_window

        real(dl) :: dz

        gaussian_window = 0._dl

        dz     = z - Win%Redshift
        winamp =  exp(-0.5_dl*(dz/Win%sigma)**2)
        gaussian_window = winamp/Win%sigma/sqrt(2._dl*const_pi)

    end function gaussian_window

    ! ---------------------------------------------------------------------------------------------
    !> This function computes a Gaussian smoothed binned window.
    !! For more details refer to the documentation.
    function bin_window(Win, z, winamp)

        implicit none

        Type(TRedWin)        :: Win      !< Window function object
        real(dl), intent(in) :: z        !< Input redshift
        real(dl)             :: winamp   !< Output amplitude of the window

        real(dl) :: bin_window

        real(dl) :: alpha, beta, z0, sigma, total_number_z

        alpha = Window_FP%fisher_cls%window_alpha
        beta  = Window_FP%fisher_cls%window_beta
        z0    = Window_FP%fisher_cls%redshift_zero
        sigma = Window_FP%fisher_cls%photoz_error*(1._dl+z)

        total_number_z = beta*z**alpha*Exp(-(z/z0)**beta)/(z0**(alpha+1._dl))/Gamma((alpha+1._dl)/beta)

        winamp = 0.5_dl*total_number_z*( &
            & +Erf( (Win%Redshift + Win%sigma - z)/(sqrt(2._dl))/sigma ) &
            & -Erf( (Win%Redshift - z)/(sqrt(2._dl))/sigma ) )

        if ( winamp<1.d-10 ) winamp = 0._dl

        bin_window = winamp

    end function bin_window

    ! ---------------------------------------------------------------------------------------------
    !> This function computes a flat binned window.
    !! For more details refer to the documentation.
    function flat_window(Win, z, winamp)

        implicit none

        Type(TRedWin)        :: Win      !< Window function object
        real(dl), intent(in) :: z        !< Input redshift
        real(dl)             :: winamp   !< Output amplitude of the window

        real(dl) :: flat_window

        real(dl) :: alpha, beta, z0, sigma, total_number_z

        winamp = 0._dl

        if ( Win%Redshift<z .and. z<Win%Redshift + Win%sigma ) then
            winamp = 1._dl
        end if

        flat_window = winamp

    end function flat_window

    ! ---------------------------------------------------------------------------------------------
    !> This function computes a flat binned window function smoothed by photo z errors.
    !! For more details refer to the documentation.
    function flat_smooth_window(Win, z, winamp)

        implicit none

        Type(TRedWin)        :: Win      !< Window function object
        real(dl), intent(in) :: z        !< Input redshift
        real(dl)             :: winamp   !< Output amplitude of the window

        real(dl) :: flat_smooth_window

        real(dl) :: alpha, beta, z0, sigma, total_number_z

        winamp = 0._dl

        if ( Win%Redshift<z .and. z<Win%Redshift + Win%sigma ) then
            winamp = 1._dl
        end if

        sigma = Window_FP%fisher_cls%photoz_error*(1._dl+z)

        winamp = 0.5_dl*winamp*( &
            & +Erf( (Win%Redshift + Win%sigma - z)/(sqrt(2._dl))/sigma ) &
            & -Erf( (Win%Redshift - z)/(sqrt(2._dl))/sigma ) )

        if ( winamp<1.d-10 ) winamp = 0._dl

        flat_smooth_window = winamp

    end function flat_smooth_window

end module Fisher_calculator_Cls_windows
