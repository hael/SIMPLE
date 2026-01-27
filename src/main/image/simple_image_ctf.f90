!@descr: for applying CTF to images
submodule (simple_image) simple_image_ctf
use simple_memoize_ft_maps
use simple_math_ctf, only: ft_map_ctf_kernel
implicit none
#include "simple_local_flags.inc"

contains
 
    !>  \brief  is for generating an image of CTF
    module subroutine ctf2img( self, tfun, dfx_in, dfy_in, angast_in )
        class(image), intent(inout) :: self
        class(ctf),   intent(inout) :: tfun
        real,         intent(in)    :: dfx_in, dfy_in, angast_in
        type(ctfparams) :: ctfvars
        integer :: h,k,physh,physk
        real    :: sum_df, diff_df, angast, amp_contr_const, wl, half_wl2_cs, tval
        ! flag image as FT
        call self%set_ft(.true.)
        ! init
        call tfun%init(dfx_in, dfy_in, angast_in) ! conversions
        ctfvars = tfun%get_ctfvars()
        wl              = ctfvars%wl
        half_wl2_cs     = 0.5 * wl * wl * ctfvars%cs
        sum_df          = ctfvars%dfx + ctfvars%dfy
        diff_df         = ctfvars%dfx - ctfvars%dfy
        angast          = ctfvars%angast
        amp_contr_const = ctfvars%amp_contr_const
        do h = ft_map_lims(1,1), ft_map_lims(1,2)
            do k = ft_map_lims(2,1), ft_map_lims(2,2)
                ! calculate CTF
                tval = ft_map_ctf_kernel(h, k, sum_df, diff_df, angast, amp_contr_const, wl, half_wl2_cs)
                ! set cmat
                physh = ft_map_phys_addrh(h,k)
                physk = ft_map_phys_addrk(h,k)
                call self%set_cmat_at(physh,physk,1, cmplx(tval,0.))
            end do
        end do
    end subroutine ctf2img

    !>  \brief  is for applying CTF to an image
    module subroutine apply_ctf_wpad( self, tfun, dfx, mode, dfy, angast, bfac )
        class(image),     intent(inout) :: self   !< instance
        class(ctf),       intent(inout) :: tfun   !< CTF object
        real,             intent(in)    :: dfx    !< defocus x-axis
        character(len=*), intent(in)    :: mode   !< abs, ctf, flip, flipneg, neg, square
        real, optional,   intent(in)    :: dfy    !< defocus y-axis
        real, optional,   intent(in)    :: angast !< angle of astigmatism
        real, optional,   intent(in)    :: bfac   !< bfactor
        integer         :: ldim(3), ldim_pd(3)
        type(ctfparams) :: ctfvars
        type(image)     :: img_pd
        ldim = self%get_ldim()
        if( self%is_3d() )then
            write(logfhandle,*) 'ldim: ', ldim
            THROW_HARD('Only 4 2D images; apply_ctf_wpad')
        endif
        ctfvars%smpd   = self%smpd
        ctfvars%dfx    = dfx
        ! defaults
        ctfvars%dfy    = ctfvars%dfx
        ctfvars%angast = 0.
        ! optionals
        if(present(dfy))    ctfvars%dfy    = dfy
        if(present(angast)) ctfvars%angast = angast
        if( self%is_ft() )then
            call self%apply_ctf(tfun, mode, ctfvars)
        else
            ldim_pd(1:2) = 2*ldim(1:2)
            ldim_pd(3)   = 1
            call img_pd%new(ldim_pd, self%smpd)
            call self%pad_mirr(img_pd)
            call img_pd%fft()
            call img_pd%apply_ctf(tfun, mode, ctfvars)
            call img_pd%ifft()
            call img_pd%clip(self)
            call img_pd%kill()
        endif
        if( present(bfac) ) call self%apply_bfac(bfac)
    end subroutine apply_ctf_wpad

    !>  \brief  is for optimised serial application of CTF
    !!          modes: abs, ctf, flip, flipneg, neg, square
    module subroutine apply_ctf( self, tfun, mode, ctfparms )
        class(image),     intent(inout) :: self     !< instance
        class(ctf),       intent(inout) :: tfun     !< CTF object
        character(len=*), intent(in)    :: mode     !< abs, ctf, flip, flipneg, neg, square
        type(ctfparams),  intent(in)    :: ctfparms !< CTF parameters
        type(ctfparams) :: ctfvars
        real, parameter :: ZERO = 0., ONE = 1.0, MIN_SQUARE = 0.001
        integer :: h,k,physh,physk
        real    :: sum_df, diff_df, dfy, angast, amp_contr_const, wl, half_wl2_cs, tval, t
        logical :: is_abs,is_ctf,is_flip,is_flipneg,is_neg,is_square
        ! Convert mode string to logical flags (SIMD-compatible)
        is_abs     = (mode == 'abs')
        is_ctf     = (mode == 'ctf')
        is_flip    = (mode == 'flip')
        is_flipneg = (mode == 'flipneg')
        is_neg     = (mode == 'neg')
        is_square  = (mode == 'square')
        ! initialize
        call tfun%init(ctfparms%dfx, ctfparms%dfy, ctfparms%angast) ! conversions
        ctfvars = tfun%get_ctfvars()
        wl              = ctfvars%wl
        half_wl2_cs     = 0.5 * wl * wl * ctfvars%cs
        sum_df          = ctfvars%dfx + ctfvars%dfy
        diff_df         = ctfvars%dfx - ctfvars%dfy
        angast          = ctfvars%angast
        amp_contr_const = ctfvars%amp_contr_const
        do h = ft_map_lims(1,1), ft_map_lims(1,2)
            do k = ft_map_lims(2,1), ft_map_lims(2,2)
                ! calculate CTF
                tval = ft_map_ctf_kernel(h, k, sum_df, diff_df, angast, amp_contr_const, wl, half_wl2_cs)
                ! SIMD-compatible mode handling using masked operations
                t = merge(       abs(tval),                     ZERO, is_abs)     + &
                    merge(           tval,                      ZERO, is_ctf)     + &
                    merge( sign(ONE, tval),                     ZERO, is_flip)    + &
                    merge(-sign(ONE, tval),                     ZERO, is_flipneg) + &
                    merge(          -tval,                      ZERO, is_neg)     + &
                    merge(min(ONE, max(tval*tval, MIN_SQUARE)), ZERO, is_square)
                ! Multiply (this is the key operation)
                physh = ft_map_phys_addrh(h,k)
                physk = ft_map_phys_addrk(h,k)
                self%cmat(physh,physk,1) = t * self%cmat(physh,physk,1)
            end do
        end do
    end subroutine apply_ctf

    ! apply CTF to image, CTF values are also returned
    module subroutine eval_and_apply_ctf( self, tfun, imode, tvalsdims, tvals, dfx_in, dfy_in, angast_in )
        class(image), intent(inout) :: self         !< instance
        class(ctf),   intent(inout) :: tfun         !< CTF object
        integer,      intent(in)    :: imode        !< CTFFLAG_FLIP=abs CTFFLAG_YES=ctf CTFFLAG_NO=no
        integer,      intent(in)    :: tvalsdims(2) !< tvals dimensions
        real,         intent(out)   :: tvals(1:tvalsdims(1),1:tvalsdims(2))
        real,         intent(in)    :: dfx_in       !< defocus x-axis
        real,         intent(in)    :: dfy_in       !< defocus y-axis
        real,         intent(in)    :: angast_in    !< angle of astigmatism
        type(ctfparams) :: ctfvars !< CTF parameters
        real    :: angast, amp_contr_const, wl, half_wl2_cs
        real    :: sum_df, diff_df,tval,t
        integer :: h,k, physh,physk
        logical :: l_flip
        if( imode == CTFFLAG_NO )then
            tvals = 1.0
            return
        endif
        ! initialize
        call tfun%init(dfx_in, dfy_in, angast_in) ! conversions
        ctfvars         = tfun%get_ctfvars()
        wl              = ctfvars%wl
        half_wl2_cs     = 0.5 * wl * wl * ctfvars%cs
        angast          = ctfvars%angast
        amp_contr_const = ctfvars%amp_contr_const    
        l_flip          = imode == CTFFLAG_FLIP
        sum_df          = ctfvars%dfx + ctfvars%dfy
        diff_df         = ctfvars%dfx - ctfvars%dfy
        do h = ft_map_lims(1,1), ft_map_lims(1,2)
            do k = ft_map_lims(2,1), ft_map_lims(2,2)
                ! calculate CTF
                tval = ft_map_ctf_kernel(h, k, sum_df, diff_df, angast, amp_contr_const, wl, half_wl2_cs)
                t    = merge(abs(tval), tval, l_flip)
                ! store tval and multiply image with tval
                physh = ft_map_phys_addrh(h,k)
                physk = ft_map_phys_addrk(h,k)
                tvals(physh, physk)      = t
                self%cmat(physh,physk,1) = t * self%cmat(physh,physk,1)
            end do
        end do
    end subroutine eval_and_apply_ctf

     ! Calculate Joe's Ice Fraction Score for images, not intended for micrographs
    module subroutine calc_ice_frac( self, tfun, ctfparms, score )
        class(image),     intent(in)    :: self
        class(ctf),       intent(inout) :: tfun
        class(ctfparams), intent(in)    :: ctfparms
        real,             intent(out)   :: score
        real, allocatable :: res(:)
        real    :: phshift, start_freq, end_freq, ice_avg, band_avg, mag_max, mag, g
        integer :: lims(3,2), box, sh, h,k, start_find, end_find, ice_maxind, hmax, kmax, cnt
        score = 0.
        if( abs(self%smpd-ctfparms%smpd) > 1.d-4) THROW_HARD('Inconsistent SMPD; calc_ice_frac')
        if( self%smpd > (ICE_BAND1/2.) ) return
        if( .not.self%is_ft() ) THROW_HARD('Image input must be in the Fourier domain!; calc_ice_frac')
        box = self%get_box()
        res = get_resarr(box, self%smpd)
        lims = self%loop_lims(2)
        call tfun%init(ctfparms%dfx, ctfparms%dfy, ctfparms%angast)
        phshift    = 0.
        start_freq = sqrt(tfun%SpaFreqSqAtNthZero(1, phshift, deg2rad(ctfparms%angast)))
        end_freq   = sqrt(tfun%SpaFreqSqAtNthZero(2, phshift, deg2rad(ctfparms%angast)))
        ice_maxind = get_find_at_res(res, ICE_BAND1)
        start_find = max(1,     ice_maxind - 3)
        end_find   = min(box/2, ice_maxind + 3)
        hmax     = -1
        kmax     = -1
        mag_max  = -1.
        band_avg = 0.
        cnt      = 0
        do k = lims(2,1),lims(2,2)
            do h = lims(1,1),lims(1,2)
                sh = nint(hyp(h,k))
                g  = real(sh) / real(box)
                if( g > start_freq .and. g < end_freq )then
                    band_avg = band_avg + csq_fast(self%get_fcomp2D(h,k))
                    cnt      = cnt + 1
                else if( sh >= start_find .and. sh < end_find )then
                    mag = csq_fast(self%get_fcomp2D(h,k))
                    if( mag > mag_max )then
                        hmax = h
                        kmax = k
                        mag_max = mag
                    endif
                endif
            end do
        end do
        if( cnt < 1 ) return
        band_avg = band_avg / real(cnt)
        ice_avg  = 0.
        do k = kmax-1,kmax+1
            do h = hmax-1,hmax+1
                ice_avg = ice_avg + csq_fast(self%get_fcomp2D(h,k))
            enddo
        enddo
        ice_avg = ice_avg / 9.
        score   = ice_avg / band_avg
        if( score > 0.5 ) score = (ice_avg - 0.5*band_avg) / band_avg
    end subroutine calc_ice_frac

    ! KEEP SERIAL
    module subroutine ctf_dens_correct( self_sum, self_rho )
        class(image), intent(inout) :: self_sum
        class(image), intent(inout) :: self_rho
        integer :: h, k, l, lims(3,2), phys(3), nyq, sh
        real    :: denom
        ! set constants
        lims = self_sum%loop_lims(2)
        nyq  = self_sum%get_lfny(1)
        do h=lims(1,1),lims(1,2)
            do k=lims(2,1),lims(2,2)
                do l=lims(3,1),lims(3,2)
                    sh    = nint(hyp(h,k,l))
                    phys  = self_sum%comp_addr_phys(h,k,l)
                    denom = real(self_rho%cmat(phys(1),phys(2),phys(3)))
                    if(sh <= nyq .and. abs(denom) > 1.e-10 )then
                        self_sum%cmat(phys(1),phys(2),phys(3)) = self_sum%cmat(phys(1),phys(2),phys(3)) / denom
                    else
                        self_sum%cmat(phys(1),phys(2),phys(3)) = cmplx(0.,0.)
                    endif
                end do
            end do
        end do
    end subroutine ctf_dens_correct

    module subroutine ctf_dens_correct_wiener( self_sum, self_rho, ssnr )
        class(image), intent(inout) :: self_sum
        class(image), intent(in)    :: self_rho
        real,         intent(in)    :: ssnr(:)
        integer :: h, k, l, lims(3,2), phys(3), nyq, sh
        real    :: denom
        ! set constants
        lims = self_sum%loop_lims(2)
        nyq  = self_sum%get_lfny(1)
        do h=lims(1,1),lims(1,2)
            do k=lims(2,1),lims(2,2)
                do l=lims(3,1),lims(3,2)
                    sh   = nint(hyp(h,k,l))
                    phys = self_sum%comp_addr_phys(h,k,l)
                    if(sh > nyq )then
                        self_sum%cmat(phys(1),phys(2),phys(3)) = cmplx(0.,0.)
                    else
                        if( sh==0 )then
                            denom = real(self_rho%cmat(phys(1),phys(2),phys(3))) + 1.
                            self_sum%cmat(phys(1),phys(2),phys(3)) = self_sum%cmat(phys(1),phys(2),phys(3))/denom
                        else
                            denom = ssnr(sh)*real(self_rho%cmat(phys(1),phys(2),phys(3))) + 1.
                            self_sum%cmat(phys(1),phys(2),phys(3)) = ssnr(sh)*self_sum%cmat(phys(1),phys(2),phys(3))/denom
                        endif
                    endif
                end do
            end do
        end do
    end subroutine ctf_dens_correct_wiener

end submodule simple_image_ctf
