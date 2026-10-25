!@descr: CTF-related math routines
module simple_math_ctf
use simple_memoize_ft_maps, only: ft_map_astigang, ft_map_spaFreqSq
use simple_defs, only: PI
implicit none

contains

    ! local fast CTF kernel (numerically equivalent refactor)
    pure elemental real function ft_map_ctf_kernel( h, k, sum_df, diff_df, angast, phshift, &
        &amp_contr_const, wl, half_wl2_cs ) result( tval )
        integer, intent(in) :: h, k
        real,    intent(in) :: sum_df, diff_df, angast, phshift, amp_contr_const, wl, half_wl2_cs
        ! phshift is canonical radians in [0,2pi); canonicalize once before entering this hot kernel.
        real :: df, PhSh, cterm, pi_wl_s2
        ! Defocus term: exactly the same expression as eval_df,
        ! but factored the sum & difference terms are memoized
        cterm   = cos( 2.0 * (ft_map_astigang(h,k) - angast) )
        df      = 0.5 * ( sum_df + cterm * diff_df )
        ! Phase shift term: preserve the same evaluation structure
        pi_wl_s2 = PI * wl * ft_map_spaFreqSq(h,k)
        PhSh     = pi_wl_s2 * (df - half_wl2_cs * ft_map_spaFreqSq(h,k))
        tval     = sin( PhSh + phshift + amp_contr_const )
    end function ft_map_ctf_kernel

    !> Relativistic electron wavelength (A) at the accelerating voltage kv (kV), the one
    !! home of the formula SIMPLE's CTF has always used (ctf%kV2wl calls it)
    elemental real function electron_wavelength( kv ) result( wavelength )
        real, intent(in) :: kv
        wavelength = 12.26 / sqrt((1000.0 * kv) + (0.9784*(1000.0 * kv)**2)/(10.0**6.0))
    end function electron_wavelength

    !> Interaction constant sigma (rad / (V A)) at kv: the phase a potential of one volt imparts
    !! over one angstrom, 2 pi m e lambda / h^2 with the relativistic mass (Kirkland, eq. 5.6)
    elemental real function electron_interaction_constant( kv ) result( sigma )
        real, intent(in) :: kv
        real(kind=8), parameter :: PLANCK            = 6.62607015d-34
        real(kind=8), parameter :: ELECTRON_MASS     = 9.1093837139d-31
        real(kind=8), parameter :: ELEMENTARY_CHARGE = 1.602176634d-19
        real(kind=8), parameter :: LIGHT_SPEED       = 2.99792458d8
        real(kind=8) :: voltage, wavelength_m, relativistic_mass
        voltage           = 1000.0d0 * real(kv, kind=8)
        wavelength_m      = 1.0d-10 * real(electron_wavelength(kv), kind=8)
        relativistic_mass = ELECTRON_MASS * (1.0d0 + ELEMENTARY_CHARGE * voltage / (ELECTRON_MASS * LIGHT_SPEED**2))
        sigma = real(2.0d0 * real(PI,kind=8) * relativistic_mass * ELEMENTARY_CHARGE * wavelength_m / PLANCK**2 * 1.0d-10)
    end function electron_interaction_constant

end module simple_math_ctf
