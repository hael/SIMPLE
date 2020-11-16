! common PRIME2D/PRIME3D routines used primarily by the Hadamard matchers
module simple_strategy2D3D_common
include 'simple_lib.f08'
use simple_image,      only: image
use simple_cmdline,    only: cmdline
use simple_builder,    only: build_glob
use simple_parameters, only: params_glob
implicit none

public :: read_img, read_imgbatch, set_bp_range, set_bp_range2D, grid_ptcl, prepimg4align,&
&norm_struct_facts, calcrefvolshift_and_mapshifts2ptcls, preprefvol,&
&prep2Dref, preprecvols, killrecvols, gen_projection_frcs, prepimgbatch,&
&build_pftcc_particles
private
#include "simple_local_flags.inc"

interface read_img
    module procedure read_img_1
    module procedure read_img_2
end interface read_img

interface read_imgbatch
    module procedure read_imgbatch_1
    module procedure read_imgbatch_2
end interface read_imgbatch

interface grid_ptcl
    module procedure grid_ptcl_1
    module procedure grid_ptcl_2
end interface grid_ptcl

real, parameter :: SHTHRESH  = 0.001
real, parameter :: CENTHRESH = 0.5    ! threshold for performing volume/cavg centering in pixels

contains

    subroutine read_img_1( iptcl )
        integer,       intent(in)     :: iptcl
        character(len=:), allocatable :: stkname
        integer :: ind_in_stk
        call build_glob%spproj%get_stkname_and_ind(params_glob%oritype, iptcl, stkname, ind_in_stk)
        call build_glob%img%read(stkname, ind_in_stk)
    end subroutine read_img_1

    subroutine read_img_2( iptcl, img )
        integer,        intent(in)    :: iptcl
        class(image),   intent(inout) :: img
        character(len=:), allocatable :: stkname
        integer :: ind_in_stk
        call build_glob%spproj%get_stkname_and_ind(params_glob%oritype, iptcl, stkname, ind_in_stk)
        call img%read(stkname, ind_in_stk)
    end subroutine read_img_2

    subroutine read_imgbatch_1( fromptop, ptcl_mask )
        integer,           intent(in) :: fromptop(2)
        logical, optional, intent(in) :: ptcl_mask(params_glob%fromp:params_glob%top)
        character(len=:), allocatable :: stkname
        integer :: iptcl, ind_in_batch, ind_in_stk
        if( present(ptcl_mask) )then
            do iptcl=fromptop(1),fromptop(2)
                if( ptcl_mask(iptcl) )then
                    ind_in_batch = iptcl - fromptop(1) + 1
                    call build_glob%spproj%get_stkname_and_ind(params_glob%oritype, iptcl, stkname, ind_in_stk)
                    call build_glob%imgbatch(ind_in_batch)%read(stkname, ind_in_stk)
                endif
            end do
        else
            do iptcl=fromptop(1),fromptop(2)
                ind_in_batch = iptcl - fromptop(1) + 1
                call build_glob%spproj%get_stkname_and_ind(params_glob%oritype, &
                    iptcl, stkname, ind_in_stk)
                call build_glob%imgbatch(ind_in_batch)%read(stkname, ind_in_stk)
            end do
        endif
    end subroutine read_imgbatch_1

    subroutine read_imgbatch_2( n, pinds, batchlims )
        integer,          intent(in)  :: n, pinds(n), batchlims(2)
        character(len=:), allocatable :: stkname
        integer :: ind_in_stk, i, ii
        do i=batchlims(1),batchlims(2)
            ii = i - batchlims(1) + 1
            call build_glob%spproj%get_stkname_and_ind(params_glob%oritype, &
                pinds(i), stkname, ind_in_stk)
            call build_glob%imgbatch(ii)%read(stkname, ind_in_stk)
        end do
    end subroutine read_imgbatch_2

    subroutine set_bp_range( cline )
        class(cmdline), intent(inout) :: cline
        real, allocatable     :: resarr(:), fsc_arr(:)
        real                  :: fsc0143, fsc05
        real                  :: mapres(params_glob%nstates)
        integer               :: s, loc(1), lp_ind, k_nyq
        character(len=STDLEN) :: fsc_fname
        logical               :: fsc_bin_exists(params_glob%nstates), all_fsc_bin_exist
        ! Nyqvist index
        k_nyq = calc_fourier_index(2.*params_glob%smpd, params_glob%boxmatch, params_glob%smpd)
        if( params_glob%l_lpset )then
            ! set Fourier index range
            params_glob%kfromto(1) = max(2, calc_fourier_index( params_glob%hp, &
                params_glob%boxmatch, params_glob%smpd))
            params_glob%kfromto(2) = calc_fourier_index(params_glob%lp, &
                params_glob%boxmatch, params_glob%smpd)
            if( cline%defined('lpstop') )then
                params_glob%kfromto(2) = min(params_glob%kfromto(2), &
                    calc_fourier_index(params_glob%lpstop, &
                    params_glob%boxmatch, params_glob%smpd))
            endif
            params_glob%kstop = params_glob%kfromto(2)
            if( params_glob%l_needs_sigma ) params_glob%kfromto(2) = k_nyq
        else
            ! check all fsc_state*.bin exist
            all_fsc_bin_exist = .true.
            fsc_bin_exists    = .false.
            do s=1,params_glob%nstates
                if( params_glob%nstates > 1 )then
                    fsc_fname = trim(CLUSTER3D_FSC)
                else
                    fsc_fname = trim(FSC_FBODY)//int2str_pad(s,2)//BIN_EXT
                endif
                fsc_bin_exists( s ) = file_exists(trim(adjustl(fsc_fname)))
                if( build_glob%spproj_field%get_pop(s, 'state') > 0 .and. .not.fsc_bin_exists(s))&
                    & all_fsc_bin_exist = .false.
            enddo
            if(build_glob%spproj%is_virgin_field(params_glob%oritype)) &
                all_fsc_bin_exist = (count(fsc_bin_exists)==params_glob%nstates)
            ! set low-pass Fourier index limit
            if( all_fsc_bin_exist )then
                resarr = build_glob%img%get_res()
                do s=1,params_glob%nstates
                    if( fsc_bin_exists(s) )then
                        ! these are the 'classical' resolution measures
                        if( params_glob%nstates > 1 )then
                            fsc_fname = trim(CLUSTER3D_FSC) ! mixed model FSC
                        else
                            fsc_fname = trim(FSC_FBODY)//int2str_pad(s,2)//BIN_EXT
                        endif
                        fsc_arr = file2rarr(trim(adjustl(fsc_fname)))
                        build_glob%fsc(s,:) = fsc_arr(:)
                        deallocate(fsc_arr)
                        call get_resolution(build_glob%fsc(s,:), resarr, fsc05, fsc0143)
                        mapres(s) = fsc0143
                    else
                        ! empty state
                        mapres(s)           = 0.
                        build_glob%fsc(s,:) = 0.
                    endif
                end do
                loc = maxloc(mapres) ! worst resolved
                if( params_glob%nstates == 1 )then
                    ! get median updatecnt
                    if( build_glob%spproj_field%median('updatecnt') > 1.0 )then ! more than half have been updated
                        lp_ind = get_lplim_at_corr(build_glob%fsc(1,:), &
                            params_glob%lplim_crit)
                    else
                        lp_ind = get_lplim_at_corr(build_glob%fsc(1,:), 0.5) ! more conservative limit @ start
                    endif
                else
                    lp_ind = get_lplim_at_corr(build_glob%fsc(loc(1),:), params_glob%lplim_crit)
                endif
                ! interpolation limit is NOT Nyqvist in correlation search
                params_glob%kfromto(2) = calc_fourier_index(resarr(lp_ind), &
                    params_glob%boxmatch, params_glob%smpd)
            else if( params_glob%l_lpset )then
                params_glob%kfromto(2) = calc_fourier_index(params_glob%lp, &
                    params_glob%boxmatch, params_glob%smpd)
            else if( build_glob%spproj_field%isthere(params_glob%fromp,'lp') )then
                params_glob%kfromto(2) = calc_fourier_index(&
                    build_glob%spproj_field%get(params_glob%fromp,'lp'), &
                    params_glob%boxmatch, params_glob%smpd)
            else
                THROW_HARD('no method available for setting the low-pass limit. Need fsc file or lp find; set_bp_range')
            endif
            ! lpstop overrides any other method for setting the low-pass limit
            if( cline%defined('lpstop') )then
                params_glob%kfromto(2) = min(params_glob%kfromto(2), &
                    calc_fourier_index(params_glob%lpstop, params_glob%boxmatch, params_glob%smpd))
            endif
            ! low-pass limit equals interpolation limit for correlation search
            params_glob%kstop = params_glob%kfromto(2)
            if( params_glob%l_needs_sigma ) params_glob%kfromto(2) = k_nyq
            ! set high-pass Fourier index limit
            params_glob%kfromto(1) = max(2,calc_fourier_index( params_glob%hp, &
                params_glob%boxmatch, params_glob%smpd))
            ! re-set the low-pass limit
            params_glob%lp = calc_lowpass_lim(params_glob%kstop, &
                params_glob%boxmatch, params_glob%smpd)
        endif
        call build_glob%spproj_field%set_all2single('lp',params_glob%lp)
    end subroutine set_bp_range

    subroutine set_bp_range2D( cline, which_iter, frac_srch_space )
        class(cmdline), intent(inout) :: cline
        integer,        intent(in)    :: which_iter
        real,           intent(in)    :: frac_srch_space
        real    :: lplim
        integer :: lpstart_find, k_nyq
        ! Nyqvist index
        k_nyq = calc_fourier_index(2.*params_glob%smpd, params_glob%boxmatch, params_glob%smpd)
        ! High-pass index
        params_glob%kfromto(1) = max(2, calc_fourier_index(params_glob%hp, params_glob%boxmatch, params_glob%smpd))
        if( params_glob%l_lpset )then
            lplim = params_glob%lp
            params_glob%kfromto(2) = calc_fourier_index(lplim, params_glob%boxmatch, params_glob%smpd)
        else
            if( file_exists(params_glob%frcs) .and. which_iter > LPLIM1ITERBOUND )then
                lplim = build_glob%projfrcs%estimate_lp_for_align()
            else
                if( which_iter <= LPLIM1ITERBOUND )then
                    lplim = params_glob%lplims2D(1)
                else if( frac_srch_space >= FRAC_SH_LIM .and. which_iter > LPLIM3ITERBOUND )then
                    lplim = params_glob%lplims2D(3)
                else
                    lplim = params_glob%lplims2D(2)
                endif
            endif
            params_glob%kfromto(2) = calc_fourier_index(lplim, params_glob%boxmatch, params_glob%smpd)
            ! to avoid pathological cases, fall-back on lpstart
            lpstart_find = calc_fourier_index(params_glob%lpstart, params_glob%boxmatch, params_glob%smpd)
            if( lpstart_find > params_glob%kfromto(2) ) params_glob%kfromto(2) = lpstart_find
        endif
        params_glob%kstop = params_glob%kfromto(2)
        if( params_glob%l_needs_sigma ) params_glob%kfromto(2) = k_nyq
        call build_glob%spproj_field%set_all2single('lp',lplim)
    end subroutine set_bp_range2D

    !>  \brief  grids one particle image to the volume
    subroutine grid_ptcl_1( fpl, se, o )
        use simple_fplane, only: fplane
        use simple_sym,    only: sym
        use simple_ori,    only: ori
        class(fplane),   intent(in) :: fpl
        class(sym),      intent(inout) :: se
        class(ori),      intent(inout) :: o
        real      :: pw
        integer   :: s, eo
        ! state flag
        s = o%get_state()
        if( s == 0 ) return
        ! eo flag
        eo = nint(o%get('eo'))
        ! particle-weight
        pw = 1.0
        if( o%isthere('w') ) pw = o%get('w')
        if( pw > TINY ) call build_glob%eorecvols(s)%grid_planes(se, o, fpl, eo, pwght=pw)
    end subroutine grid_ptcl_1

    !>  \brief  grids one particle image to the volume (distribution of weigted oris)
    subroutine grid_ptcl_2( fpl, se, o, os )
        use simple_fplane, only: fplane
        use simple_sym,    only: sym
        use simple_ori,    only: ori
        use simple_oris,   only: oris
        class(fplane), intent(in)    :: fpl
        class(sym),    intent(inout) :: se
        class(ori),    intent(inout) :: o
        class(oris),   intent(inout) :: os
        real, allocatable :: states(:)
        real    :: pw
        integer :: s, eo
        ! eo flag
        eo = nint(o%get('eo'))
        ! particle-weight
        pw = 1.0
        if( o%isthere('w') ) pw = o%get('w')
        if( pw > TINY )then
            ! gridding
            if( params_glob%nstates == 1 )then
                call build_glob%eorecvols(1)%grid_planes(se, os, fpl, eo, pwght=pw)
            else
                states = os%get_all('state')
                do s=1,params_glob%nstates
                    if( count(nint(states) == s) > 0 )then
                        call build_glob%eorecvols(s)%grid_planes(se, os, fpl, eo, pwght=pw, state=s)
                    endif
                end do
            endif
        endif
    end subroutine grid_ptcl_2

    !>  \brief  prepares all particle images for alignment
    subroutine build_pftcc_particles( pftcc, batchsz_max, match_imgs, ptcl_mask )
        use simple_polarft_corrcalc, only: polarft_corrcalc
        use simple_polarizer,        only: polarizer
        class(polarft_corrcalc), intent(inout) :: pftcc
        integer,                 intent(in)    :: batchsz_max
        class(polarizer),        intent(inout) :: match_imgs(batchsz_max)
        logical, optional,       intent(in)    :: ptcl_mask(params_glob%fromp:params_glob%top)
        logical :: mask_here(params_glob%fromp:params_glob%top)
        integer :: iptcl_batch, batchlims(2), imatch, iptcl
        if( present(ptcl_mask) )then
            mask_here = ptcl_mask
        else
            mask_here = .true.
        endif
        if( .not. params_glob%l_distr_exec ) write(logfhandle,'(A)') '>>> BUILDING PARTICLES'
        call prepimgbatch( batchsz_max )
        do iptcl_batch=params_glob%fromp,params_glob%top,batchsz_max
            batchlims = [iptcl_batch,min(params_glob%top,iptcl_batch + batchsz_max - 1)]
            call read_imgbatch( batchlims, mask_here )
            !$omp parallel do default(shared) private(iptcl,imatch)&
            !$omp schedule(static) proc_bind(close)
            do iptcl=batchlims(1),batchlims(2)
                if( .not. mask_here(iptcl) ) cycle
                imatch = iptcl - batchlims(1) + 1
                call prepimg4align( iptcl, build_glob%imgbatch(imatch), match_imgs(imatch))
                ! transfer to polar coordinates
                call match_imgs(imatch)%polarize(pftcc, iptcl, .true., .true., mask=build_glob%l_resmsk)
            end do
            !$omp end parallel do
        end do
    end subroutine build_pftcc_particles

    !>  \brief  prepares one particle image for alignment
    !!          serial routine
    subroutine prepimg4align( iptcl, img_in, img_out )
        use simple_polarizer,     only: polarizer
        use simple_estimate_ssnr, only: fsc2optlp_sub
        use simple_ctf,           only: ctf
        integer,          intent(in)    :: iptcl
        class(image),     intent(inout) :: img_in
        class(polarizer), intent(inout) :: img_out
        type(ctf)       :: tfun
        type(ctfparams) :: ctfparms
        real            :: x, y, sdev_noise
        x = build_glob%spproj_field%get(iptcl, 'x')
        y = build_glob%spproj_field%get(iptcl, 'y')
        ! CTF parameters
        ctfparms = build_glob%spproj%get_ctfparams(params_glob%oritype, iptcl)
        ! normalise
        call img_in%noise_norm(build_glob%lmsk, sdev_noise)
        ! move to Fourier space
        call img_in%fft()
        ! Shift image to rotational origin & phase-flipping
        if(abs(x) > SHTHRESH .or. abs(y) > SHTHRESH) call img_in%shift2Dserial([-x,-y])
        select case(ctfparms%ctfflag)
            case(CTFFLAG_NO, CTFFLAG_FLIP)
                ! all good
            case(CTFFLAG_YES) ! phase flip
                tfun = ctf(ctfparms%smpd, ctfparms%kv, ctfparms%cs, ctfparms%fraca)
                call tfun%apply_serial(img_in, 'flip', ctfparms)
            case DEFAULT
                THROW_HARD('unsupported CTF flag: '//int2str(ctfparms%ctfflag)//' prepimg4align')
        end select
        ! back to real-space
        call img_in%ifft()
        ! clip image if needed
        call img_in%clip(img_out)
        ! soft-edged mask
        if( params_glob%l_innermsk )then
            call img_out%mask(params_glob%msk, 'soft', inner=params_glob%inner, width=params_glob%width)
        else
            if( params_glob%l_focusmsk )then
                call img_out%mask(params_glob%focusmsk, 'soft')
            else
                call img_out%mask(params_glob%msk, 'soft')
            endif
        endif
        ! gridding prep
        if( params_glob%griddev.eq.'yes' )then
            call img_out%div_by_instrfun
        endif
        ! return in Fourier space
        call img_out%fft()
        call img_out%mul(real(params_glob%boxmatch**2) / real(params_glob%box**2))
    end subroutine prepimg4align

    !>  \brief  prepares one cluster centre image for alignment
    subroutine prep2Dref( pftcc, img_in, img_out, icls, center, xyz_in, xyz_out )
        use simple_polarft_corrcalc, only: polarft_corrcalc
        use simple_estimate_ssnr,    only: fsc2optlp_sub, subsample_optlp, subsample_filter
        use simple_polarizer,        only: polarizer
        class(polarft_corrcalc), intent(inout) :: pftcc
        class(image),            intent(inout) :: img_in
        class(polarizer),        intent(inout) :: img_out
        integer,           intent(in)    :: icls
        logical, optional, intent(in)    :: center
        real,    optional, intent(in)    :: xyz_in(3)
        real,    optional, intent(out)   :: xyz_out(3)
        real    :: frc(build_glob%projfrcs%get_filtsz()), filter(build_glob%projfrcs%get_filtsz())
        real    :: subfilter(build_glob%img_match%get_filtsz())
        real    :: xyz(3), sharg
        integer :: filtsz
        logical :: do_center
        filtsz    = build_glob%projfrcs%get_filtsz()
        do_center = (params_glob%center .eq. 'yes')
        ! centering only performed if params_glob%center.eq.'yes'
        if( present(center) ) do_center = do_center .and. center
        if( do_center )then
            if( present(xyz_in) )then
                sharg = arg(xyz_in)
                if( sharg > CENTHRESH )then
                    ! apply shift and do NOT update the corresponding class parameters
                    call img_in%fft()
                    call img_in%shift2Dserial(xyz_in(1:2))
                endif
            else
                xyz = img_in%calc_shiftcen_serial(params_glob%cenlp, params_glob%msk)
                sharg = arg(xyz)
                if( sharg > CENTHRESH )then
                    ! apply shift and update the corresponding class parameters
                    call img_in%fft()
                    call img_in%shift2Dserial(xyz(1:2))
                    call build_glob%spproj_field%add_shift2class(icls, -xyz(1:2))
                endif
                if( present(xyz_out) ) xyz_out = xyz
            endif
        endif
        ! filter
        if( params_glob%l_pssnr )then
            call build_glob%projpssnrs%frc_getter(icls, params_glob%hpind_fsc, params_glob%l_phaseplate, filter)
            call subsample_filter(filtsz, build_glob%img_match%get_filtsz(), filter, subfilter)
            call pftcc%set_ref_optlp(icls, subfilter(params_glob%kfromto(1):params_glob%kstop))
        else
            call build_glob%projfrcs%frc_getter(icls, params_glob%hpind_fsc, params_glob%l_phaseplate, frc)
            if( any(frc > 0.143) )then
                call fsc2optlp_sub(build_glob%projfrcs%get_filtsz(), frc, filter)
                if( params_glob%l_match_filt )then
                    call subsample_optlp(filtsz, build_glob%img_match%get_filtsz(), filter, subfilter)
                    call pftcc%set_ref_optlp(icls, subfilter(params_glob%kfromto(1):params_glob%kstop))
                else
                    call img_in%fft() ! needs to be here in case the shift was never applied (above)
                    call img_in%apply_filter_serial(filter)
                endif
            endif
        endif
        ! ensure we are in real-space before clipping
        call img_in%ifft()
        ! clip image if needed
        call img_in%clip(img_out)
        ! apply mask
        if( params_glob%l_innermsk )then
            call img_out%mask(params_glob%msk, 'soft', &
                inner=params_glob%inner, width=params_glob%width)
        else
            call img_out%mask(params_glob%msk, 'soft')
        endif
        ! gridding prep
        if( params_glob%griddev.eq.'yes' )then
            call img_out%div_by_instrfun
        endif
        ! move to Fourier space
        call img_out%fft()
        call img_out%mul(real(params_glob%boxmatch**2) / real(params_glob%box**2))
    end subroutine prep2Dref

    !>  \brief  initializes all volumes for reconstruction
    subroutine preprecvols( wcluster )
        real, optional, intent(in)    :: wcluster
        character(len=:), allocatable :: part_str
        real,    allocatable :: resarr(:)
        integer, allocatable :: pops(:)
        real    :: lplim_rec, fsc05, fsc0143
        integer :: istate
        allocate(part_str, source=int2str_pad(params_glob%part,params_glob%numlen))
        call build_glob%spproj_field%get_pops(pops, 'state')
        lplim_rec = huge(lplim_rec)
        resarr    = build_glob%img%get_res()
        do istate = 1, params_glob%nstates
            if( pops(istate) > 0)then
                call build_glob%eorecvols(istate)%new( build_glob%spproj)
                call build_glob%eorecvols(istate)%reset_all
                if( params_glob%l_frac_update )then
                    call build_glob%eorecvols(istate)%read_eos(trim(VOL_FBODY)//&
                        int2str_pad(istate,2)//'_part'//part_str)
                    call build_glob%eorecvols(istate)%expand_exp
                    call build_glob%eorecvols(istate)%apply_weight(1.0 - &
                        params_glob%update_frac)
                endif
                ! determining resolution for low-pass limited reconstruction
                if( any(build_glob%fsc(istate,:) > 0.143) )then
                    call get_resolution(build_glob%fsc(istate,:), resarr, fsc05, fsc0143)
                    lplim_rec = min(lplim_rec, fsc0143)
                endif
            endif
        end do
        deallocate(pops,resarr)
    end subroutine preprecvols

    !>  \brief  destructs all volumes for reconstruction
    subroutine killrecvols
        integer :: istate
        do istate = 1, params_glob%nstates
            call build_glob%eorecvols(istate)%kill
        end do
    end subroutine killrecvols

    !>  \brief  prepares a batch of image
    subroutine prepimgbatch( batchsz )
        integer,        intent(in)    :: batchsz
        integer :: currsz, ibatch
        logical :: doprep
        if( .not. allocated(build_glob%imgbatch) )then
            doprep = .true.
        else
            currsz = size(build_glob%imgbatch)
            if( batchsz > currsz )then
                do ibatch=1,currsz
                    call build_glob%imgbatch(ibatch)%kill
                end do
                deallocate(build_glob%imgbatch)
                doprep = .true.
            else
                doprep = .false.
            endif
        endif
        if( doprep )then
            allocate(build_glob%imgbatch(batchsz))
            do ibatch=1,batchsz
                call build_glob%imgbatch(ibatch)%new([params_glob%box,params_glob%box,1], &
                    params_glob%smpd, wthreads=.false.)
            end do
        endif
    end subroutine prepimgbatch

    !>  \brief  determines the reference volume shift and map shifts back to particles
    !>          reference volume shifting is performed in preprefvol
    subroutine calcrefvolshift_and_mapshifts2ptcls(cline, s, volfname, do_center, xyz )
        class(cmdline),   intent(inout) :: cline
        integer,          intent(in)    :: s
        character(len=*), intent(in)    :: volfname
        logical,          intent(out)   :: do_center
        real,             intent(out)   :: xyz(3)
        logical :: has_been_searched
        do_center = .true.
        ! ensure correct build_glob%vol dim
        call build_glob%vol%new([params_glob%box,params_glob%box,params_glob%box],params_glob%smpd)
        ! centering
        if( params_glob%center .eq. 'no' .or. params_glob%nstates > 1 .or. &
            .not. params_glob%l_doshift .or. params_glob%pgrp(:1) .ne. 'c' .or. &
            cline%defined('mskfile') .or. params_glob%l_frac_update )then
            do_center = .false.
            xyz       = 0.
            return
        endif
        call build_glob%vol%read(volfname)
        xyz = build_glob%vol%calc_shiftcen(params_glob%cenlp,params_glob%msk)
        if( params_glob%pgrp .ne. 'c1' ) xyz(1:2) = 0.     ! shifts only along z-axis for C2 and above
        if( arg(xyz) <= CENTHRESH )then
            do_center = .false.
            xyz       = 0.
            return
        endif
        ! map back to particle oritentations
        has_been_searched = .not.build_glob%spproj%is_virgin_field(params_glob%oritype)
        if( has_been_searched ) call build_glob%spproj_field%map3dshift22d(-xyz(:), state=s)
    end subroutine calcrefvolshift_and_mapshifts2ptcls

    !>  \brief  prepares one volume for references extraction
    subroutine preprefvol( pftcc, cline, s, volfname, do_center, xyz, iseven )
        use simple_polarft_corrcalc, only: polarft_corrcalc
        use simple_estimate_ssnr,    only: fsc2optlp_sub, subsample_optlp, subsample_filter
        use simple_ori,              only: ori
        class(polarft_corrcalc), intent(inout) :: pftcc
        class(cmdline),          intent(inout) :: cline
        integer,                 intent(in)    :: s
        character(len=*),        intent(in)    :: volfname
        logical,                 intent(in)    :: do_center
        real,                    intent(in)    :: xyz(3)
        logical,                 intent(in)    :: iseven
        type(image)                   :: mskvol
        type(ori)                     :: o
        character(len=:), allocatable :: fname_vol_filter
        real,             allocatable :: pssnr(:)
        real    :: subfilter(build_glob%img_match%get_filtsz())
        real    :: filter(build_glob%projfrcs%get_filtsz()), frc(build_glob%projfrcs%get_filtsz())
        integer :: iref, iproj, iprojred, filtsz, subfiltsz
        ! ensure correct build_glob%vol dim
        call build_glob%vol%new([params_glob%box,params_glob%box,params_glob%box],params_glob%smpd)
        call build_glob%vol%read(volfname)
        if( do_center )then
            call build_glob%vol%fft()
            call build_glob%vol%shift([xyz(1),xyz(2),xyz(3)])
        endif
        ! Volume filtering
        if( .not.params_glob%l_lpset )then
            filtsz    = build_glob%projfrcs%get_filtsz()
            subfiltsz = build_glob%img_match%get_filtsz()
            if( params_glob%l_match_filt )then
                ! stores filters in pftcc
                if( params_glob%l_pssnr )then
                    allocate(fname_vol_filter, source=PSSNR_FBODY//int2str_pad(s,2)//BIN_EXT)
                    subfilter = 1.
                    if( any(build_glob%fsc(s,:) > 0.143) .and. file_exists(fname_vol_filter))then
                        pssnr = file2rarr(fname_vol_filter)
                        call subsample_filter(filtsz, subfiltsz, pssnr, subfilter)
                    endif
                    do iref = (s-1)*params_glob%nspace+1, s*params_glob%nspace
                        call pftcc%set_ref_optlp(iref, subfilter(params_glob%kfromto(1):params_glob%kstop))
                    enddo
                else
                    if( file_exists(params_glob%frcs) )then
                        if( params_glob%clsfrcs.eq.'yes')then
                            iproj = 0
                            do iref = 1,2*build_glob%projfrcs%get_nprojs()
                                iproj = iproj+1
                                if( iproj > build_glob%projfrcs%get_nprojs() ) iproj = 1
                                call build_glob%projfrcs%frc_getter(iproj, params_glob%hpind_fsc, params_glob%l_phaseplate, frc)
                                call fsc2optlp_sub(filtsz, frc, filter)
                                call subsample_optlp(filtsz, subfiltsz, filter, subfilter)
                                call pftcc%set_ref_optlp(iref, subfilter(params_glob%kfromto(1):params_glob%kstop))
                            enddo
                        else
                            iproj = 0
                            do iref = (s-1)*params_glob%nspace+1, s*params_glob%nspace
                                iproj    = iproj+1
                                call build_glob%eulspace%get_ori(iproj, o)
                                iprojred = build_glob%eulspace_red%find_closest_proj(o)
                                call build_glob%projfrcs%frc_getter(iprojred, params_glob%hpind_fsc, params_glob%l_phaseplate, frc)
                                call fsc2optlp_sub(filtsz, frc, filter)
                                call subsample_optlp(filtsz, subfiltsz, filter, subfilter)
                                call pftcc%set_ref_optlp(iref, subfilter(params_glob%kfromto(1):params_glob%kstop))
                            enddo
                        endif
                    else
                        if( any(build_glob%fsc(s,:) > 0.143) )then
                            call fsc2optlp_sub(filtsz, build_glob%fsc(s,:), filter)
                            call subsample_optlp(filtsz, subfiltsz, filter, subfilter)
                        else
                            subfilter = 1.
                        endif
                        do iref = (s-1)*params_glob%nspace+1, s*params_glob%nspace
                            call pftcc%set_ref_optlp(iref, subfilter(params_glob%kfromto(1):params_glob%kstop))
                        enddo
                    endif
                endif
            else
                if( params_glob%cc_objfun == OBJFUN_EUCLID )then
                    call build_glob%vol%fft() ! needs to be here in case the shift was never applied (above)
                    if( params_glob%nstates == 1 )then
                        allocate(fname_vol_filter, source=trim(ANISOLP_FBODY)//int2str_pad(s,2)//trim(params_glob%ext))
                    else
                        allocate(fname_vol_filter, source=trim(CLUSTER3D_ANISOLP)//trim(params_glob%ext))
                    endif
                    if( file_exists(fname_vol_filter) )then
                        ! anisotropic filter
                        call build_glob%vol2%read(fname_vol_filter)
                        call build_glob%vol%apply_filter(build_glob%vol2)
                    else
                        ! match filter based on Rosenthal & Henderson, 2003
                        if( any(build_glob%fsc(s,:) > 0.143) )then
                            call fsc2optlp_sub(filtsz,build_glob%fsc(s,:),filter)
                            call build_glob%vol%apply_filter(filter)
                        endif
                    endif
                else
                    call build_glob%vol%fft() ! needs to be here in case the shift was never applied (above)
                    if( params_glob%nstates == 1 )then
                        allocate(fname_vol_filter, source=trim(ANISOLP_FBODY)//int2str_pad(s,2)//trim(params_glob%ext))
                    else
                        allocate(fname_vol_filter, source=trim(CLUSTER3D_ANISOLP)//trim(params_glob%ext))
                    endif
                    if( file_exists(fname_vol_filter) )then
                        ! anisotropic filter
                        call build_glob%vol2%read(fname_vol_filter)
                        call build_glob%vol%apply_filter(build_glob%vol2)
                    else
                        ! match filter based on Rosenthal & Henderson, 2003
                        if( any(build_glob%fsc(s,:) > 0.143) )then
                            call fsc2optlp_sub(filtsz,build_glob%fsc(s,:),filter)
                            call build_glob%vol%apply_filter(filter)
                        endif
                    endif
                endif
            endif
        endif
        ! back to real space
        call build_glob%vol%ifft()
        ! clip
        if( params_glob%boxmatch < params_glob%box ) &
            call build_glob%vol%clip_inplace(&
            [params_glob%boxmatch,params_glob%boxmatch,params_glob%boxmatch])
        ! masking
        if( cline%defined('mskfile') )then
            ! mask provided
            call mskvol%new([params_glob%box, params_glob%box, params_glob%box], &
                params_glob%smpd)
            call mskvol%read(params_glob%mskfile)
            if( params_glob%boxmatch < params_glob%box )&
                call mskvol%clip_inplace(&
                [params_glob%boxmatch,params_glob%boxmatch,params_glob%boxmatch])
            if( cline%defined('lp_backgr') )then
                call build_glob%vol%lp_background(mskvol, params_glob%lp_backgr)
            else
                call build_glob%vol%zero_env_background(mskvol)
                call build_glob%vol%mul(mskvol)
            endif
            call mskvol%kill
        else
            ! circular masking
            if( params_glob%l_innermsk )then
                call build_glob%vol%mask(params_glob%msk, 'soft', &
                    inner=params_glob%inner, width=params_glob%width)
            else
                call build_glob%vol%mask(params_glob%msk, 'soft')
            endif
        endif
        ! gridding prep
        if( params_glob%griddev.eq.'yes' )then
            call build_glob%vol%div_w_instrfun(params_glob%alpha)
        endif
        ! FT volume
        call build_glob%vol%fft()
        ! expand for fast interpolation
        call build_glob%vol%expand_cmat(params_glob%alpha,norm4proj=.true.)
        call o%kill
    end subroutine preprefvol

    subroutine norm_struct_facts( cline, which_iter )
        use simple_filterer, only: gen_anisotropic_optlp
        class(cmdline),    intent(inout) :: cline
        integer, optional, intent(in)    :: which_iter
        integer               :: s, find4eoavg
        real                  :: res05s(params_glob%nstates), res0143s(params_glob%nstates)
        character(len=STDLEN) :: pprocvol
        call build_glob%vol%new([params_glob%box,params_glob%box,params_glob%box],params_glob%smpd)
        call build_glob%vol2%new([params_glob%box,params_glob%box,params_glob%box],params_glob%smpd)
        ! init
        res0143s = 0.
        res05s   = 0.
        ! cycle through states
        do s=1,params_glob%nstates
            if( build_glob%spproj_field%get_pop(s, 'state') == 0 )then
                ! empty state
                if( present(which_iter) ) build_glob%fsc(s,:) = 0.
                cycle
            endif
            call build_glob%eorecvols(s)%compress_exp
            if( params_glob%l_distr_exec )then
                call build_glob%eorecvols(s)%write_eos('recvol_state'//int2str_pad(s,2)//'_part'//&
                    int2str_pad(params_glob%part,params_glob%numlen))
            else
                if( present(which_iter) )then
                    params_glob%vols(s) = 'recvol_state'//int2str_pad(s,2)//'_iter'//&
                        int2str_pad(which_iter,3)//params_glob%ext
                else
                    params_glob%vols(s) = 'startvol_state'//int2str_pad(s,2)//params_glob%ext
                endif
                params_glob%vols_even(s) = add2fbody(params_glob%vols(s), params_glob%ext, '_even')
                params_glob%vols_odd(s)  = add2fbody(params_glob%vols(s), params_glob%ext, '_odd')
                call build_glob%eorecvols(s)%sum_eos
                call build_glob%eorecvols(s)%sampl_dens_correct_eos(s, params_glob%vols_even(s), &
                    &params_glob%vols_odd(s), find4eoavg)
                call gen_projection_frcs(cline,  params_glob%vols_even(s), params_glob%vols_odd(s), &
                    params_glob%mskfile, s, build_glob%projfrcs)
                call build_glob%projfrcs%write(FRCS_FILE)
                call gen_anisotropic_optlp(build_glob%vol2, build_glob%projfrcs, &
                    build_glob%eulspace_red, s, params_glob%pgrp, params_glob%hpind_fsc, &
                    params_glob%l_phaseplate)
                call build_glob%vol2%write('aniso_optlp_state'//int2str_pad(s,2)//params_glob%ext)
                call build_glob%eorecvols(s)%get_res(res05s(s), res0143s(s))
                call build_glob%eorecvols(s)%sampl_dens_correct_sum(build_glob%vol)
                call build_glob%vol%write(params_glob%vols(s), del_if_exists=.true.)
                 ! need to put the sum back at lowres for the eo pairs
                call build_glob%vol%fft()
                call build_glob%vol2%zero_and_unflag_ft
                call build_glob%vol2%read(params_glob%vols_even(s))
                call build_glob%vol2%fft()
                call build_glob%vol2%insert_lowres(build_glob%vol, find4eoavg)
                call build_glob%vol2%ifft()
                call build_glob%vol2%write(params_glob%vols_even(s), del_if_exists=.true.)
                call build_glob%vol2%zero_and_unflag_ft
                call build_glob%vol2%read(params_glob%vols_odd(s))
                call build_glob%vol2%fft()
                call build_glob%vol2%insert_lowres(build_glob%vol, find4eoavg)
                call build_glob%vol2%ifft()
                call build_glob%vol2%write(params_glob%vols_odd(s), del_if_exists=.true.)
                if( present(which_iter) )then
                    ! post-process volume
                    pprocvol   = add2fbody(trim(params_glob%vols(s)), params_glob%ext, PPROC_SUFFIX)
                    build_glob%fsc(s,:) = file2rarr('fsc_state'//int2str_pad(s,2)//'.bin')
                    ! low-pass filter
                    call build_glob%vol%bp(0., params_glob%lp)
                    call build_glob%vol%ifft()
                    ! mask
                    call build_glob%vol%mask(params_glob%msk, 'soft')
                    call build_glob%vol%write(pprocvol)
                else
                    call build_glob%vol%zero_and_unflag_ft
                endif
            endif
        end do
        if( .not. params_glob%l_distr_exec )then
            ! set the resolution limit according to the worst resolved model
            params_glob%lp = min(params_glob%lp,max(params_glob%lpstop,maxval(res0143s)))
        endif
        call build_glob%vol%kill
        call build_glob%vol2%kill
    end subroutine norm_struct_facts

    !>  \brief generate projection FRCs from even/odd pairs
    subroutine gen_projection_frcs( cline, ename, oname, resmskname, state, projfrcs )
        use simple_oris,            only: oris
        use simple_projector_hlev,  only: reproject
        use simple_projection_frcs, only: projection_frcs
        class(cmdline),         intent(inout) :: cline
        character(len=*),       intent(in)    :: ename, oname, resmskname
        integer,                intent(in)    :: state
        class(projection_frcs), intent(inout) :: projfrcs
        type(oris)               :: e_space
        type(image)              :: mskvol
        type(image), allocatable :: even_imgs(:), odd_imgs(:)
        real,        allocatable :: frc(:)
        integer :: iproj, find_plate
        ! ensure correct build_glob%vol dim
        call build_glob%vol%new([params_glob%box,params_glob%box,params_glob%box],params_glob%smpd)
        ! read & prep even/odd pair
        call build_glob%vol%read(ename)
        call build_glob%vol2%read(oname)
        call mskvol%new([params_glob%box, params_glob%box, params_glob%box], params_glob%smpd)
        call prepeovol(build_glob%vol)
        call prepeovol(build_glob%vol2)
        ! create e_space
        call e_space%new(NSPACE_REDUCED)
        call build_glob%pgrpsyms%build_refspiral(e_space)
        ! generate even/odd projections
        even_imgs = reproject(build_glob%vol,  e_space)
        odd_imgs  = reproject(build_glob%vol2, e_space)
        ! calculate FRCs and fill-in projfrcs object
        allocate(frc(even_imgs(1)%get_filtsz()))
        !$omp parallel do default(shared) private(iproj,frc) schedule(static) proc_bind(close)
        do iproj=1,NSPACE_REDUCED
            call even_imgs(iproj)%fft()
            call odd_imgs(iproj)%fft()
            call even_imgs(iproj)%fsc(odd_imgs(iproj), frc)
            if( params_glob%l_phaseplate ) call phaseplate_correct_fsc(frc, find_plate)
            call projfrcs%set_frc(iproj, frc, state)
            call even_imgs(iproj)%kill
            call odd_imgs(iproj)%kill
        end do
        !$omp end parallel do
        deallocate(even_imgs, odd_imgs)
        call e_space%kill
        call mskvol%kill

        contains

            !>  \brief  prepares even/odd volume for FSC/FRC calcualtion
            subroutine prepeovol( vol )
                class(image), intent(inout) :: vol
                if( params_glob%l_envfsc .and. cline%defined('mskfile') )then
                    ! mask provided
                    call mskvol%read(resmskname)
                    call vol%zero_env_background(mskvol)
                    call vol%mul(mskvol)
                else
                    ! circular masking
                    if( params_glob%l_innermsk )then
                        call vol%mask(params_glob%msk, 'soft', &
                            inner=params_glob%inner, width=params_glob%width)
                    else
                        call vol%mask(params_glob%msk, 'soft')
                    endif
                endif
        end subroutine prepeovol

    end subroutine gen_projection_frcs

end module simple_strategy2D3D_common
