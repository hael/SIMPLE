!@descr: CTF-related math routines
module simple_math_ctf
use simple_memoize_ft_maps, only: ft_map_astigang, ft_map_spaFreqSq
use simple_defs, only: PI
implicit none

contains

    ! local fast CTF kernel (numerically equivalent refactor)
    pure elemental real function ft_map_ctf_kernel( h, k, sum_df, diff_df, angast, amp_contr_const, wl, half_wl2_cs ) result( tval )
        integer, intent(in) :: h, k
        real,    intent(in) :: sum_df, diff_df, angast, amp_contr_const, wl, half_wl2_cs
        real :: df, PhSh, cterm, pi_wl_s2
        ! Defocus term: exactly the same expression as eval_df,
        ! but factored the sum & difference terms are memoized
        cterm   = cos( 2.0 * (ft_map_astigang(h,k) - angast) )
        df      = 0.5 * ( sum_df + cterm * diff_df )
        ! Phase shift term: preserve the same evaluation structure
        pi_wl_s2 = PI * wl * ft_map_spaFreqSq(h,k)
        PhSh     = pi_wl_s2 * (df - half_wl2_cs * ft_map_spaFreqSq(h,k))
        tval     = sin( PhSh + amp_contr_const )
    end function ft_map_ctf_kernel

end module simple_math_ctf