! Fourier filtering routines
module simple_filterer
use simple_defs
use simple_image,     only: image
use simple_projector, only: projector
use simple_syslib,    only: alloc_errchk
implicit none

real, private, parameter :: SHTHRESH=0.0001

contains

    !> DOSE FILTERING (Grant, Grigorieff eLife 2015)
    !! input is template image, accumulative dose (in e/A2) and acceleration voltage
    !!         output is filter coefficients
    !! \f$  \mathrm{dose}_\mathrm{acc} = \int^{N}_{1} \mathrm{dose\_weight}(a,F,V),\ n_\mathrm{e}/\si{\angstrom\squared}  \f$
    !! \param acc_dose accumulative dose (in \f$n_\mathrm{e}^- per \si{\angstrom\squared}\f$)
    function acc_dose2filter( img, acc_dose, kV ) result( filter )
        type(image), intent(in) :: img           !< input image
        real,        intent(in) :: acc_dose, kV  !< acceleration voltage
        real, allocatable       :: filter(:)
        integer :: find, sz
        sz = img%get_filtsz()
        allocate(filter(sz),stat=alloc_stat)
        call alloc_errchk("simple_filterer::acc_dose2filter ",alloc_stat)
        do find=1,sz
            filter(find) = dose_weight(acc_dose, img%get_spat_freq(find), kV)
        end do
    end function acc_dose2filter

    !>  \brief Calculate dose weight. Input is accumulative dose (in e/A2) and
    !>  spatial frequency (in 1/A)
    !!         output is resolution dependent weight applied to individual frames
    !!         before correlation search and averaging
    !!
    !! \f$  \mathrm{dose\_weight}(a,F,V) = \exp\left(- \frac{A_\mathrm{dose}}{k\times (2.0\times A\times f^B + C)} \right), \f$
    !! where \f$k\f$ is 0.75 for \f$V<200\f$ kV, 1.0 for \f$200 \leqslant  V \leqslant  300\f$
    real function dose_weight( acc_dose, spat_freq, kV )
        real, intent(in) :: acc_dose                !< accumulative dose (in e/A2)
        real, intent(in) :: spat_freq               !< spatial frequency (in 1/A)
        real, intent(in) :: kV                      !< accelleration voltage
        real, parameter  :: A=0.245, B=-1.665, C=2.81, kV_factor=0.75
        real             :: critical_exp !< critical exposure (only depends on spatial frequency)
        critical_exp = A*(spat_freq**B)+C
        if( abs(kV-300.) < 0.001 )then
            ! critical exposure does not need modification
        else if( abs(kV-200.) < 0.001 )then
            ! critical exposure at 200 kV expected to be ~25% lower
            critical_exp = critical_exp*kV_factor
        else
            stop 'unsupported kV (acceleration voltage); simple_filterer :: dose_weight'
        endif
        dose_weight = exp(-acc_dose/(2.0*critical_exp))
    end function dose_weight

    !> \brief  re-samples a filter array
    function resample_filter( filt_orig, res_orig, res_new ) result( filt_resamp )
        use simple_math, only: find
        real, intent(in)  :: filt_orig(:), res_orig(:), res_new(:)
        real, allocatable :: filt_resamp(:)                 !< output filter array
        integer :: filtsz_orig, filtsz_resamp, k, ind
        real    :: dist
        filtsz_orig   = size(filt_orig)
        filtsz_resamp = size(res_new)
        allocate(filt_resamp(filtsz_resamp),stat=alloc_stat)
        call alloc_errchk("simple_filterer::resample_filter ",alloc_stat)
        do k=1,filtsz_resamp
            call find(res_orig, filtsz_orig, res_new(k), ind, dist)
            filt_resamp(k) = filt_orig(ind)
        end do
    end function resample_filter

    !> WIENER RESTORATION ROUTINES

    !>  \brief does the Wiener restoration of aligned images in 2D
    !!   only for testing
    subroutine wiener_restore2D( img_set, o_set, tfplan, img_rec, msk )
        use simple_oris,  only: oris
        use simple_ori,   only: ori
        class(image),     intent(inout) :: img_set(:) !< input images
        class(oris),      intent(inout) :: o_set      !< set of oris objects
        type(ctfplan),    intent(in)    :: tfplan     !< CTF plan
        class(image),     intent(inout) :: img_rec     !< reconstructed image
        real,             intent(in)    :: msk         !< mask
        integer           :: ldim(3), ldim_pad(3), nimgs, iptcl
        type(ori)         :: o
        type(image)       :: ctfsqsum
        real              :: smpd         !< sampling distance
        if( o_set%get_noris() /= size(img_set) )&
        stop 'nr of imgs and oris not consistent; simple_filterer :: wiener_restore2D_1'
        ! set constants
        ldim     = img_set(1)%get_ldim()
        smpd     = img_set(1)%get_smpd()
        ldim_pad = img_rec%get_ldim()
        nimgs    = size(img_set)
        ! create & init objs
        call img_rec%new(ldim_pad, smpd)
        call ctfsqsum%new(ldim_pad, smpd)
        ctfsqsum = cmplx(0.,0.)
        ! average in the assumption of infinite signal
        do iptcl=1,nimgs
            o = o_set%get_ori(iptcl)
            call wiener_restore2D_online(img_set(iptcl), o,&
            &tfplan, img_rec, ctfsqsum, msk)
        end do
        ! do the density correction
        call img_rec%fwd_ft
        call img_rec%ctf_dens_correct(ctfsqsum)
        call img_rec%bwd_ft
        ! destroy objects
        call ctfsqsum%kill
    end subroutine wiener_restore2D

    !>  \brief does the online Wiener restoration of 2D images, including shift+rotations
    !!         the image is left shifted and Fourier transformed on output
    subroutine wiener_restore2D_online( img, o, tfplan, img_rec, ctfsqsum, msk, add )
        use simple_ori,            only: ori
        use simple_ctf,            only: ctf
        use simple_projector_hlev, only: rotimg
        class(image),      intent(inout) :: img       !< input image
        class(ori),        intent(inout) :: o         !< ori object
        type(ctfplan),     intent(in)    :: tfplan    !< CTF plan object
        class(image),      intent(inout) :: img_rec   !< reconstructed image
        class(image),      intent(inout) :: ctfsqsum  !< instrument CTF filter
        real,              intent(in)    :: msk       !< mask
        logical, optional, intent(in)    :: add       !< add or subtr  rotation and ctfsumsq
        type(image) :: roimg, ctfsq
        type(ctf)   :: tfun
        integer     :: ldim(3)
        real        :: angast, dfx, dfy, x, y, smpd, w
        logical     :: aadd
        aadd = .true.
        if( present(add) ) aadd = add
        ! set constants
        ldim = img%get_ldim()
        smpd = img%get_smpd()
        ! create & init objs
        call ctfsq%new(ldim, smpd)
        call ctfsq%set_ft(.true.)
        if( tfplan%flag .ne. 'no' )&
        tfun = ctf(img%get_smpd(), o%get('kv'), o%get('cs'), o%get('fraca'))
        ! set CTF and shift parameters
        select case(tfplan%mode)
            case('astig') ! astigmatic CTF
                dfx    = o%get('dfx')
                dfy    = o%get('dfy')
                angast = o%get('angast')
            case('noastig') ! non-astigmatic CTF
                dfx    = o%get('dfx')
                dfy    = dfx
                angast = 0.
        end select
        x = -o%get('x')
        y = -o%get('y')
        w = o%get('w')
        ! apply
        call img%fwd_ft
        ! take care of the nominator
        select case(tfplan%flag)
            case('yes')  ! multiply with CTF
                call tfun%apply_and_shift(img, ctfsq, x, y, dfx, 'ctf', dfy, angast)
            case('flip') ! multiply with abs(CTF)
                call tfun%apply_and_shift(img, ctfsq, x, y, dfx, 'abs', dfy, angast)
            case('mul','no')
                call tfun%apply_and_shift(img, ctfsq, x, y, dfx, '', dfy, angast)
        end select
        ! griding-based image rotation and filtering
        call rotimg(img, -o%e3get(), msk, roimg)
        ! assemble img_rec sum
        if( aadd )then
            call img_rec%add(roimg, w)
            call ctfsqsum%add(ctfsq, w)
        else
            call img_rec%subtr(roimg, w)
            call ctfsqsum%subtr(ctfsq, w)
        endif
        call roimg%kill
        call ctfsq%kill
    end subroutine wiener_restore2D_online

    subroutine fit_bfac( img_ref, img_ptcl, o, tfplan, bfac_range, lp, msk, bfac_best )
        use simple_ctf, only: ctf
        use simple_ori, only: ori
        class(image),  intent(in)    :: img_ref, img_ptcl
        class(ori),    intent(inout) :: o
        type(ctfplan), intent(in)    :: tfplan
        real,          intent(in)    :: bfac_range(2), lp, msk
        real,          intent(out)   :: bfac_best
        character(len=STDLEN) :: mode
        type(ctf)             :: tfun
        type(image)           :: ref, ptcl, maskimg, ref_tmp
        real                  :: dfx, dfy, angast, bfac_range_refine(2), smpd, corr_best, bfac, cc
        integer               :: npix, ldim(3), imode
        real, parameter       :: BFAC_STEPSZ = 5.0, BFAC_STEPSZ_REFINE = 1.0
        if( .not. (img_ref.eqdims.img_ptcl)) stop 'ref & ptcl imgs not of same dims; filterer :: fit_bfac'
        ! extract image info
        ldim = img_ref%get_ldim()
        smpd = img_ref%get_smpd()
        ! extract CTF info
        select case(tfplan%flag)
            case('yes')  ! multiply with CTF
                mode  ='ctf'
                imode = 1
            case('flip') ! multiply with abs(CTF)
                mode  = 'abs'
                imode = 2
            case('mul','no')
                mode  = ''
                imode = 3
        end select
        dfx    = o%get('dfx')
        dfy    = o%get('dfy')
        angast = o%get('angast')
        ! make hard mask
        call maskimg%disc(ldim, smpd, msk, npix)
        ! make CTF object
        tfun = ctf(smpd, o%get('kv'), o%get('cs'), o%get('fraca'))
        ! prepare particle image
        ptcl = img_ptcl
        call ptcl%bp(0.,lp)
        call ptcl%bwd_ft
        ! prepare reference image
        ref = img_ref
        call ref%fwd_ft
        ! init
        corr_best = -1.
        bfac_best =  0.
        ! coarse search
        bfac = bfac_range(1)
        do while( bfac < bfac_range(2) )
            call iteration
            bfac = bfac + BFAC_STEPSZ
        end do
        ! refinement
        bfac_range_refine(1) = bfac_best - BFAC_STEPSZ
        bfac_range_refine(2) = bfac_best + BFAC_STEPSZ
        bfac = bfac_range_refine(1)
        do while( bfac < bfac_range_refine(2) )
            call iteration
             bfac = bfac + BFAC_STEPSZ_REFINE
        end do

        contains

            subroutine iteration
                ref_tmp = ref
                if( imode < 3 )then
                    call tfun%apply(ref_tmp, dfx, trim(mode), dfy, angast, bfac)
                else
                    call ref_tmp%apply_bfac(bfac)
                endif
                call ref_tmp%bp(0.,lp)
                call ref_tmp%bwd_ft
                cc = ref_tmp%real_corr(ptcl, maskimg)
                if( cc > corr_best )then
                    corr_best = cc
                    bfac_best = bfac
                endif
            end subroutine iteration

    end subroutine fit_bfac

end module simple_filterer
