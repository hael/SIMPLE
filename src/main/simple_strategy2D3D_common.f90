! common PRIME2D/PRIME3D routines used primarily by the Hadamard matchers
module simple_strategy2D3D_common
include 'simple_lib.f08'
use simple_image,      only: image
use simple_cmdline,    only: cmdline
use simple_builder,    only: build_glob
use simple_parameters, only: params_glob
implicit none

public :: read_img, read_imgbatch, set_bp_range, set_bp_range2D, grid_ptcl, prepimg4align,&
&eonorm_struct_facts, norm_struct_facts, calcrefvolshift_and_mapshifts2ptcls, preprefvol,&
&prep2Dref, gen2Dclassdoc, preprecvols, killrecvols, gen_projection_frcs, prepimgbatch,&
&build_pftcc_particles
private
#include "simple_local_flags.inc"

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

    subroutine read_img( iptcl )
        integer,       intent(in)     :: iptcl
        character(len=:), allocatable :: stkname
        integer :: ind_in_stk
        call build_glob%spproj%get_stkname_and_ind(params_glob%oritype, iptcl, stkname, ind_in_stk)
        call build_glob%img%read(stkname, ind_in_stk)
    end subroutine read_img

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
                call build_glob%spproj%get_stkname_and_ind(params_glob%oritype, iptcl, stkname, ind_in_stk)
                call build_glob%imgbatch(ind_in_batch)%read(stkname, ind_in_stk)
            end do
        endif
    end subroutine read_imgbatch_1

    subroutine read_imgbatch_2( n, pinds, batchlims )
        integer,       intent(in)     :: n, pinds(n), batchlims(2)
        character(len=:), allocatable :: stkname
        integer :: ind_in_stk, i, ii
        do i=batchlims(1),batchlims(2)
            ii = i - batchlims(1) + 1
            call build_glob%spproj%get_stkname_and_ind(params_glob%oritype, pinds(i), stkname, ind_in_stk)
            call build_glob%imgbatch(ii)%read(stkname, ind_in_stk)
        end do
    end subroutine read_imgbatch_2

    subroutine set_bp_range( cline )
        class(cmdline), intent(inout) :: cline
        real, allocatable     :: resarr(:), fsc_arr(:)
        real                  :: fsc0143, fsc05, mapres(params_glob%nstates)
        integer               :: s, loc(1), lp_ind, k_nyq
        character(len=STDLEN) :: fsc_fname
        logical               :: fsc_bin_exists(params_glob%nstates), all_fsc_bin_exist
        ! Nyqvist index
        k_nyq = calc_fourier_index(2.*params_glob%smpd, params_glob%boxmatch, params_glob%smpd)
        select case(params_glob%eo)
            case('yes','aniso')
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
                if(build_glob%spproj%is_virgin_field(params_glob%oritype)) all_fsc_bin_exist = (count(fsc_bin_exists)==params_glob%nstates)
                ! Fourier index range for corr_valid
                params_glob%kfromto_valid(1) = calc_fourier_index(HP_CORR_VALID, params_glob%boxmatch, params_glob%smpd)
                params_glob%kfromto_valid(2) = min(calc_fourier_index(LP_CORR_VALID, params_glob%boxmatch, params_glob%smpd), k_nyq)
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
                            lp_ind = get_lplim_at_corr(build_glob%fsc(1,:), params_glob%lplim_crit)
                        else
                            lp_ind = get_lplim_at_corr(build_glob%fsc(1,:), 0.5) ! more conservative limit @ start
                        endif
                    else
                        lp_ind = get_lplim_at_corr(build_glob%fsc(loc(1),:), params_glob%lplim_crit)
                    endif
                    ! interpolation limit is NOT Nyqvist in correlation search
                    params_glob%kfromto(2) = calc_fourier_index(resarr(lp_ind), params_glob%boxmatch, params_glob%smpd)
                else if( cline%defined('lp') )then
                    params_glob%kfromto(2) = calc_fourier_index(params_glob%lp, params_glob%boxmatch, params_glob%smpd)
                else if( build_glob%spproj_field%isthere(params_glob%fromp,'lp') )then
                    params_glob%kfromto(2) = calc_fourier_index(build_glob%spproj_field%get(params_glob%fromp,'lp'), params_glob%boxmatch, params_glob%smpd)
                else
                    THROW_HARD('no method available for setting the low-pass limit. Need fsc file or lp find; set_bp_range')
                endif
                ! lpstop overrides any other method for setting the low-pass limit
                if( cline%defined('lpstop') )then
                    params_glob%kfromto(2) = min(params_glob%kfromto(2), calc_fourier_index(params_glob%lpstop, params_glob%boxmatch, params_glob%smpd))
                endif
                ! low-pass limit equals interpolation limit for correlation search
                params_glob%kstop = params_glob%kfromto(2)
                ! possible extension of interpolation limit to accomodate corr_valid
                params_glob%kfromto(2) = max(params_glob%kfromto(2), params_glob%kfromto_valid(2))
                ! for the euclidean distance case all frequencies need to be extracted
                if( params_glob%cc_objfun .eq. OBJFUN_EUCLID) params_glob%kfromto(2) = k_nyq
                ! set high-pass Fourier index limit
                params_glob%kfromto(1) = max(2,calc_fourier_index( params_glob%hp, params_glob%boxmatch, params_glob%smpd))
                ! re-set the low-pass limit
                params_glob%lp = calc_lowpass_lim(params_glob%kstop, params_glob%boxmatch, params_glob%smpd)
                call build_glob%spproj_field%set_all2single('lp',params_glob%lp)
            case('no')
                ! set Fourier index range
                params_glob%kfromto(1) = max(2, calc_fourier_index( params_glob%hp, params_glob%boxmatch, params_glob%smpd))
                params_glob%kfromto(2) = calc_fourier_index(params_glob%lp, params_glob%boxmatch, params_glob%smpd)
                if( cline%defined('lpstop') )then
                    params_glob%kfromto(2) = min(params_glob%kfromto(2), calc_fourier_index(params_glob%lpstop, params_glob%boxmatch, params_glob%smpd))
                endif
                params_glob%kstop = params_glob%kfromto(2)
                call build_glob%spproj_field%set_all2single('lp',params_glob%lp)
            case DEFAULT
                THROW_HARD('Unsupported eo flag')
        end select
        DebugPrint '*** simple_strategy2D3D_common ***: did set Fourier index range'
    end subroutine set_bp_range

    subroutine set_bp_range2D( cline, which_iter, frac_srch_space )
        class(cmdline), intent(inout) :: cline
        integer,        intent(in)    :: which_iter
        real,           intent(in)    :: frac_srch_space
        real    :: lplim
        integer :: lpstart_find, k_nyq
        ! Nyqvist index
        k_nyq = calc_fourier_index(2.*params_glob%smpd, params_glob%boxmatch, params_glob%smpd)
        ! Fourier index range for corr_valid
        params_glob%kfromto_valid(1) = calc_fourier_index(HP_CORR_VALID, params_glob%boxmatch, params_glob%smpd)
        params_glob%kfromto_valid(2) = min(calc_fourier_index(LP_CORR_VALID, params_glob%boxmatch, params_glob%smpd), k_nyq)
        ! High-pass index
        params_glob%kfromto(1) = max(2, calc_fourier_index(params_glob%hp, params_glob%boxmatch, params_glob%smpd))
        if( cline%defined('lp') )then
            params_glob%kfromto(2) = calc_fourier_index(params_glob%lp, params_glob%boxmatch, params_glob%smpd)
            call build_glob%spproj_field%set_all2single('lp',params_glob%lp)
            params_glob%kstop      = params_glob%kfromto(2)
            ! possible extension of interpolation limit to accomodate corr_valid
            params_glob%kfromto(2) = max(params_glob%kfromto(2), params_glob%kfromto_valid(2))
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
            params_glob%kstop = params_glob%kfromto(2)
            ! possible extension of interpolation limit to accomodate corr_valid
            params_glob%kfromto(2) = max(params_glob%kfromto(2), params_glob%kfromto_valid(2))
            call build_glob%spproj_field%set_all2single('lp',lplim)
        endif
        DebugPrint  '*** simple_strategy2D3D_common ***: did set Fourier index range'
    end subroutine set_bp_range2D

    !>  \brief  grids one particle image to the volume
    subroutine grid_ptcl_1( img, se, o, ctfvars )
        use simple_sym, only: sym
        use simple_ori,      only: ori
        class(image),    intent(inout) :: img
        class(sym),      intent(inout) :: se
        class(ori),      intent(inout) :: o
        type(ctfparams), intent(in)    :: ctfvars
        real      :: pw, bfac_rec
        integer   :: s, eo
        eo = 0
        if( params_glob%eo .ne. 'no' ) eo = nint(o%get('eo'))
        if( params_glob%shellw .eq. 'yes' )then
            ! shell-weighted reconstruction
            bfac_rec = 0.
            if( o%isthere('bfac_rec') ) bfac_rec = max(0., o%get('bfac_rec'))
            ! state flag
            s = o%get_state()
            ! fwd ft
            call img%fft()
            ! gridding
            if( params_glob%eo .ne. 'no' )then
                call build_glob%eorecvols(s)%grid_fplane(se, o, ctfvars, img, eo, pwght=1., bfac=bfac_rec)
            else
                call build_glob%recvols(s)%insert_fplane(se, o, ctfvars, img, pwght=1., bfac=bfac_rec)
            endif
        else
            ! particle-weighted reconstruction
            pw = 1.0
            if( o%isthere('w') ) pw = o%get('w')
            if( pw > TINY )then
                ! state flag
                s = o%get_state()
                ! fwd ft
                call img%fft()
                ! gridding
                if( params_glob%eo .ne. 'no' )then
                    call build_glob%eorecvols(s)%grid_fplane(se, o, ctfvars, img, eo, pwght=pw)
                else
                    call build_glob%recvols(s)%insert_fplane(se, o, ctfvars, img, pwght=pw)
                endif
            endif
        endif
    end subroutine grid_ptcl_1

    !>  \brief  grids one particle image to the volume (distribution of weigted oris)
    subroutine grid_ptcl_2( img, se, o, os, ctfvars )
        use simple_sym,  only: sym
        use simple_ori,  only: ori
        use simple_oris, only: oris
        class(image),    intent(inout) :: img
        class(sym),      intent(inout) :: se
        class(ori),      intent(inout) :: o
        class(oris),     intent(inout) :: os
        type(ctfparams), intent(in)    :: ctfvars
        real, allocatable :: states(:)
        real    :: pw, bfac_rec
        integer :: s, eo
        eo = 0
        if( params_glob%eo .ne. 'no' ) eo = nint(o%get('eo'))
        if( params_glob%shellw.eq.'yes' )then
            ! shell-weighted reconstruction
            bfac_rec = 0.
            if( o%isthere('bfac_rec') ) bfac_rec = o%get('bfac_rec')
            ! fwd ft
            call img%fft()
            ! gridding
            if( params_glob%nstates == 1 )then
                if( params_glob%eo .ne. 'no' )then
                    call build_glob%eorecvols(1)%grid_fplane(se, os, ctfvars, img, eo, pwght=1., bfac=bfac_rec)
                else
                    call build_glob%recvols(1)%insert_fplane(se, os, ctfvars, img, pwght=1., bfac=bfac_rec)
                endif
            else
                states = os%get_all('state')
                do s=1,params_glob%nstates
                    if( count(nint(states) == s) > 0 )then
                        if( params_glob%eo .ne. 'no' )then
                            call build_glob%eorecvols(s)%grid_fplane(se, os, ctfvars, img, eo, pwght=1., bfac=bfac_rec, state=s)
                        else
                            call build_glob%recvols(s)%insert_fplane(se, os, ctfvars, img, pwght=1., bfac=bfac_rec, state=s)
                        endif
                    endif
                end do
            endif
        else
            ! particle-weighted reconstruction
            pw = 1.0
            if( o%isthere('w') ) pw = o%get('w')
            if( pw > TINY )then
                ! fwd ft
                call img%fft()
                ! gridding
                if( params_glob%nstates == 1 )then
                    if( params_glob%eo .ne. 'no' )then
                        call build_glob%eorecvols(1)%grid_fplane(se, os, ctfvars, img, eo, pwght=pw)
                    else
                        call build_glob%recvols(1)%insert_fplane(se, os, ctfvars, img, pwght=pw)
                    endif
                else
                    states = os%get_all('state')
                    do s=1,params_glob%nstates
                        if( count(nint(states) == s) > 0 )then
                            if( params_glob%eo .ne. 'no' )then
                                call build_glob%eorecvols(s)%grid_fplane(se, os, ctfvars, img, eo, pwght=pw, state=s)
                            else
                                call build_glob%recvols(s)%insert_fplane(se, os, ctfvars, img, pwght=pw, state=s)
                            endif
                        endif
                    end do
                endif
            endif
        endif
    end subroutine grid_ptcl_2

    !>  \brief  prepares all particle images for alignment
    subroutine build_pftcc_particles( pftcc, batchsz_max, match_imgs, is3D,  ptcl_mask )
        use simple_polarft_corrcalc, only: polarft_corrcalc
        use simple_polarizer,        only: polarizer
        class(polarft_corrcalc), intent(inout) :: pftcc
        integer,                 intent(in)    :: batchsz_max
        class(polarizer),        intent(inout) :: match_imgs(batchsz_max)
        logical,                 intent(in)    :: is3D
        logical, optional,       intent(in)    :: ptcl_mask(params_glob%fromp:params_glob%top)
        logical :: mask_here(params_glob%fromp:params_glob%top)
        integer :: iptcl_batch, batchlims(2), imatch, iptcl
        if( present(ptcl_mask) )then
            mask_here = ptcl_mask
        else
            mask_here = .true.
        endif
        if( .not. params_glob%l_distr_exec ) write(*,'(A)') '>>> BUILDING PARTICLES'
        call prepimgbatch( batchsz_max )
        do iptcl_batch=params_glob%fromp,params_glob%top,batchsz_max
            batchlims = [iptcl_batch,min(params_glob%top,iptcl_batch + batchsz_max - 1)]
            call read_imgbatch( batchlims, mask_here )
            !$omp parallel do default(shared) private(iptcl,imatch)&
            !$omp schedule(static) proc_bind(close)
            do iptcl=batchlims(1),batchlims(2)
                if( .not. mask_here(iptcl) ) cycle
                imatch = iptcl - batchlims(1) + 1
                call prepimg4align( iptcl, build_glob%imgbatch(imatch), match_imgs(imatch), is3D=is3D)
                ! transfer to polar coordinates
                call match_imgs(imatch)%polarize(pftcc, iptcl, .true., .true.)
            end do
            !$omp end parallel do
        end do
    end subroutine build_pftcc_particles

    !>  \brief  prepares one particle image for alignment
    !!          serial routine
    subroutine prepimg4align( iptcl, img_in, img_out, is3D )
        use simple_polarizer,     only: polarizer
        use simple_estimate_ssnr, only: fsc2optlp_sub
        use simple_ctf,           only: ctf
        integer,          intent(in)    :: iptcl
        class(image),     intent(inout) :: img_in
        class(polarizer), intent(inout) :: img_out
        logical,          intent(in)    :: is3D
        type(ctf)       :: tfun
        type(ctfparams) :: ctfparms
        real            :: frc(build_glob%projfrcs%get_filtsz()), filter(build_glob%projfrcs%get_filtsz()), x, y
        integer         :: ifrc
        ! shift
        x = build_glob%spproj_field%get(iptcl, 'x')
        y = build_glob%spproj_field%get(iptcl, 'y')
        ! CTF parameters
        ctfparms = build_glob%spproj%get_ctfparams(params_glob%oritype, iptcl)
        ! normalise
        call img_in%norm()
        ! move to Fourier space
        call img_in%fft()
        ! shell normalization & filter
        if( is3D )then
            if( (trim(params_glob%eo).ne.'no') .and. (params_glob%nstates==1) )then
                if( trim(params_glob%clsfrcs).eq.'yes' )then
                    ifrc = nint(build_glob%spproj_field%get(iptcl,'class'))
                else
                    ifrc = build_glob%eulspace_red%find_closest_proj( build_glob%spproj_field%get_ori(iptcl) )
                endif
            else
                ifrc = 0 ! turns off the matched filter
            endif
        else
            if( params_glob%l_match_filt )then
                ifrc = nint(build_glob%spproj_field%get(iptcl,'class'))
            else
                ifrc = 0 ! turns off the matched filter
            endif
        endif
        if( ifrc > 0 .and. params_glob%l_match_filt )then
            call build_glob%projfrcs%frc_getter(ifrc, params_glob%hpind_fsc, params_glob%l_phaseplate, frc)
            if( any(frc > 0.143) )then
                call fsc2optlp_sub(build_glob%projfrcs%get_filtsz(), frc, filter)
                call img_in%shellnorm_and_apply_filter_serial(filter)
            endif
        endif
        ! Shift image to rotational origin & phase-flipping
        select case(ctfparms%ctfflag)
            case(CTFFLAG_NO, CTFFLAG_FLIP)  ! shift only
                if(abs(x) > SHTHRESH .or. abs(y) > SHTHRESH) call img_in%shift2Dserial([-x,-y])
            case(CTFFLAG_YES)               ! phase flip & shift
                tfun = ctf(ctfparms%smpd, ctfparms%kv, ctfparms%cs, ctfparms%fraca)
                if(abs(x) > SHTHRESH .or. abs(y) > SHTHRESH)then
                    call tfun%phaseflip_and_shift_serial(img_in, -x, -y, ctfparms)
                else
                    call tfun%apply_serial(img_in, 'flip', ctfparms)
                endif
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
        ! return in Fourier space
        call img_out%fft()
    end subroutine prepimg4align

    !>  \brief  prepares one cluster centre image for alignment
    subroutine prep2Dref( img_in, img_out, icls, center, xyz_in, xyz_out )
        use simple_estimate_ssnr, only: fsc2optlp_sub
        use simple_polarizer,     only: polarizer
        class(image),      intent(inout) :: img_in
        class(polarizer),  intent(inout) :: img_out
        integer,           intent(in)    :: icls
        logical, optional, intent(in)    :: center
        real,    optional, intent(in)    :: xyz_in(3)
        real,    optional, intent(out)   :: xyz_out(3)
        real    :: frc(build_glob%projfrcs%get_filtsz()), filter(build_glob%projfrcs%get_filtsz())
        real    :: xyz(3), sharg
        logical :: do_center
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
                xyz = img_in%center_serial(params_glob%cenlp, params_glob%msk)
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
        if( params_glob%l_match_filt )then
            ! anisotropic matched filter
            call build_glob%projfrcs%frc_getter(icls, params_glob%hpind_fsc, params_glob%l_phaseplate, frc)
            if( any(frc > 0.143) )then
                call img_in%fft() ! needs to be here in case the shift was never applied (above)
                call fsc2optlp_sub(build_glob%projfrcs%get_filtsz(), frc, filter)
                call img_in%shellnorm_and_apply_filter_serial(filter)
            endif
        endif
        ! ensure we are in real-space before clipping
        call img_in%ifft()
        ! clip image if needed
        call img_in%clip(img_out)
        ! apply mask
        if( params_glob%l_innermsk )then
            call img_out%mask(params_glob%msk, 'soft', inner=params_glob%inner, width=params_glob%width)
        else
            call img_out%mask(params_glob%msk, 'soft')
        endif
        ! move to Fourier space
        call img_out%fft()
    end subroutine prep2Dref

    !>  \brief prepares a 2D class document with class index, resolution,
    !!         poulation, average correlation and weight
    subroutine gen2Dclassdoc
        integer, allocatable :: pops(:)
        integer    :: icls, pop
        real       :: frc05, frc0143
        call build_glob%spproj%os_cls2D%new(params_glob%ncls)
        call build_glob%spproj_field%get_pops(pops, 'class', maxn=params_glob%ncls)
        do icls=1,params_glob%ncls
            call build_glob%projfrcs%estimate_res(icls, frc05, frc0143)
            pop = pops(icls)
            call build_glob%spproj%os_cls2D%set(icls, 'class', real(icls))
            call build_glob%spproj%os_cls2D%set(icls, 'pop',   real(pop))
            call build_glob%spproj%os_cls2D%set(icls, 'res',   frc0143)
            if( pop > 0 )then
                call build_glob%spproj%os_cls2D%set(icls, 'state', 1.0) ! needs to be default val if no selection has been done
            else
                call build_glob%spproj%os_cls2D%set(icls, 'state', 0.0) ! exclusion
            endif
            if( pop > 1 )then
                call build_glob%spproj%os_cls2D%set(icls, 'corr',  build_glob%spproj_field%get_avg('corr', class=icls))
                call build_glob%spproj%os_cls2D%set(icls, 'w',     build_glob%spproj_field%get_avg('w',    class=icls))
            else
                call build_glob%spproj%os_cls2D%set(icls, 'corr', -1.0)
                call build_glob%spproj%os_cls2D%set(icls, 'w',     0.0)
            endif
        end do
        deallocate(pops)
    end subroutine gen2Dclassdoc

    !>  \brief  initializes all volumes for reconstruction
    subroutine preprecvols( wcluster )
        real, optional, intent(in)    :: wcluster
        character(len=:), allocatable :: recname, rhoname, part_str
        real,    allocatable :: resarr(:)
        integer, allocatable :: pops(:)
        real    :: lplim_rec, fsc05, fsc0143
        integer :: istate
        allocate(part_str, source=int2str_pad(params_glob%part,params_glob%numlen))
        call build_glob%spproj_field%get_pops(pops, 'state')
        select case(params_glob%eo)
            case('yes','aniso')
                lplim_rec = huge(lplim_rec)
                resarr    = build_glob%img%get_res()
                do istate = 1, params_glob%nstates
                    if( pops(istate) > 0)then
                        call build_glob%eorecvols(istate)%new( build_glob%spproj)
                        call build_glob%eorecvols(istate)%reset_all
                        if( params_glob%l_frac_update )then
                            call build_glob%eorecvols(istate)%read_eos(trim(VOL_FBODY)//int2str_pad(istate,2)//'_part'//part_str)
                            call build_glob%eorecvols(istate)%expand_exp
                            call build_glob%eorecvols(istate)%apply_weight(1.0 - params_glob%update_frac)
                        endif
                        ! determining resolution for low-pass limited reconstruction
                        if( any(build_glob%fsc(istate,:) > 0.143) )then
                            call get_resolution(build_glob%fsc(istate,:), resarr, fsc05, fsc0143)
                            lplim_rec = min(lplim_rec, fsc0143)
                        endif
                    endif
                end do
                !!!!!!!! TAKE THIS OUT TO TURN OFF LOW-PASS LIMITED REC !!!!!!!!!!!!!!!!!
                ! do istate = 1, params_glob%nstates
                !     if( pops(istate) > 0)call build_glob%eorecvols(istate)%set_lplim(lplim_rec)
                ! end do
                !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                deallocate(resarr)
            case DEFAULT
                lplim_rec = params_glob%lp
                do istate = 1, params_glob%nstates
                    if( pops(istate) > 0)then
                        call build_glob%recvols(istate)%new([params_glob%boxpd, params_glob%boxpd, params_glob%boxpd], params_glob%smpd)
                        call build_glob%recvols(istate)%alloc_rho(build_glob%spproj)
                        call build_glob%recvols(istate)%reset
                        call build_glob%recvols(istate)%reset_exp
                        call build_glob%recvols(istate)%set_lplim(lplim_rec)
                        if( params_glob%l_frac_update )then
                            allocate(recname, source=trim(VOL_FBODY)//int2str_pad(istate,2)//'_part'//part_str//params_glob%ext)
                            allocate(rhoname, source='rho_'//trim(VOL_FBODY)//int2str_pad(istate,2)//'_part'//part_str//params_glob%ext)
                            if( file_exists(recname) .and. file_exists(rhoname) )then
                                call build_glob%recvols(istate)%read(recname)
                                call build_glob%recvols(istate)%read_rho(rhoname)
                                call build_glob%recvols(istate)%expand_exp
                                call build_glob%recvols(istate)%apply_weight(1.0 - params_glob%update_frac)
                            endif
                            deallocate(recname, rhoname)
                        endif
                    endif
                end do
        end select
        deallocate(pops)
    end subroutine preprecvols

    !>  \brief  destructs all volumes for reconstruction
    subroutine killrecvols
        integer :: istate
        if( params_glob%eo .ne. 'no' )then
            do istate = 1, params_glob%nstates
                call build_glob%eorecvols(istate)%kill
            end do
        else
            do istate = 1, params_glob%nstates
                call build_glob%recvols(istate)%dealloc_rho
                call build_glob%recvols(istate)%kill
            end do
        endif
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
                call build_glob%imgbatch(ibatch)%new([params_glob%box,params_glob%box,1], params_glob%smpd)
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
        if( params_glob%center .eq. 'no' .or. params_glob%nstates > 1 .or. .not. params_glob%l_doshift .or.&
        &params_glob%pgrp(:1) .ne. 'c' .or. cline%defined('mskfile') .or. params_glob%l_frac_update )then
            do_center = .false.
            xyz       = 0.
            return
        endif
        call build_glob%vol%read(volfname)
        xyz = build_glob%vol%center(params_glob%cenlp,params_glob%msk) ! find center of mass shift
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
    subroutine preprefvol( cline, s, volfname, do_center, xyz )
        use simple_estimate_ssnr, only: fsc2optlp
        class(cmdline),   intent(inout) :: cline
        integer,          intent(in)    :: s
        character(len=*), intent(in)    :: volfname
        logical,          intent(in)    :: do_center
        real,             intent(in)    :: xyz(3)
        type(image)                     :: mskvol
        real,             allocatable   :: filter(:)
        character(len=:), allocatable   :: fname_vol_filter
        ! ensure correct build_glob%vol dim
        call build_glob%vol%new([params_glob%box,params_glob%box,params_glob%box],params_glob%smpd)
        call build_glob%vol%read(volfname)
        if( do_center )then
            call build_glob%vol%fft()
            call build_glob%vol%shift([xyz(1),xyz(2),xyz(3)])
        endif
        ! Volume filtering
        if( params_glob%eo.ne.'no' )then
            ! anisotropic matched filter
            if( params_glob%nstates.eq.1 )then
                allocate(fname_vol_filter, source=trim(ANISOLP_FBODY)//int2str_pad(s,2)//trim(params_glob%ext))
                call build_glob%vol%fft() ! needs to be here in case the shift was never applied (above)
                if( file_exists(fname_vol_filter) )then
                    call build_glob%vol2%read(fname_vol_filter)
                    if( params_glob%l_match_filt )then
                        call build_glob%vol%shellnorm_and_apply_filter(build_glob%vol2)
                    else
                        call build_glob%vol%apply_filter(build_glob%vol2)
                    endif
                else
                    ! matched filter based on Rosenthal & Henderson, 2003
                    if( any(build_glob%fsc(s,:) > 0.143) )then
                        filter = fsc2optlp(build_glob%fsc(s,:))
                        if( params_glob%l_match_filt )then
                            call build_glob%vol%shellnorm_and_apply_filter(filter)
                        else
                            call build_glob%vol%apply_filter(filter)
                        endif
                    endif
                endif
            else
                ! no filtering involved for multiple states
            endif
        endif
        ! back to real space
        call build_glob%vol%ifft()
        ! clip
        if( params_glob%boxmatch < params_glob%box ) call build_glob%vol%clip_inplace([params_glob%boxmatch,params_glob%boxmatch,params_glob%boxmatch])
        ! masking
        if( cline%defined('mskfile') )then
            ! mask provided
            call mskvol%new([params_glob%box, params_glob%box, params_glob%box], params_glob%smpd)
            call mskvol%read(params_glob%mskfile)
            if( params_glob%boxmatch < params_glob%box ) call mskvol%clip_inplace([params_glob%boxmatch,params_glob%boxmatch,params_glob%boxmatch])
            call build_glob%vol%zero_background
            call build_glob%vol%mul(mskvol)
            call mskvol%kill
        else
            ! circular masking
            if( params_glob%l_innermsk )then
                call build_glob%vol%mask(params_glob%msk, 'soft', inner=params_glob%inner, width=params_glob%width)
            else
                call build_glob%vol%mask(params_glob%msk, 'soft')
            endif
        endif
        ! FT volume
        call build_glob%vol%fft()
        ! expand for fast interpolation
        call build_glob%vol%expand_cmat(params_glob%alpha)
    end subroutine preprefvol

    subroutine norm_struct_facts( which_iter )
        integer, optional, intent(in)    :: which_iter
        integer :: s
        character(len=:), allocatable :: fbody
        character(len=STDLEN) :: pprocvol
        call build_glob%vol%new([params_glob%box,params_glob%box,params_glob%box],params_glob%smpd)
        do s=1,params_glob%nstates
            if( build_glob%spproj_field%get_pop(s, 'state') == 0 )then
                ! empty space
                cycle
            endif
            if( params_glob%l_distr_exec )then
                allocate(fbody, source='recvol_state'//int2str_pad(s,2)//'_part'//int2str_pad(params_glob%part,params_glob%numlen))
                params_glob%vols(s) = trim(adjustl(fbody))//params_glob%ext
                call build_glob%recvols(s)%compress_exp
                call build_glob%recvols(s)%write(params_glob%vols(s), del_if_exists=.true.)
                call build_glob%recvols(s)%write_rho('rho_'//trim(adjustl(fbody))//params_glob%ext)
                deallocate(fbody)
            else
                if( params_glob%refine .eq. 'snhc' )then
                     params_glob%vols(s) = trim(SNHCVOL)//int2str_pad(s,2)//params_glob%ext
                else
                    if( present(which_iter) )then
                        params_glob%vols(s) = 'recvol_state'//int2str_pad(s,2)//'_iter'//int2str_pad(which_iter,3)//params_glob%ext
                    else
                        params_glob%vols(s) = 'startvol_state'//int2str_pad(s,2)//params_glob%ext
                    endif
                endif
                call build_glob%recvols(s)%compress_exp
                call build_glob%recvols(s)%sampl_dens_correct
                call build_glob%recvols(s)%ifft()
                call build_glob%recvols(s)%clip(build_glob%vol)
                call build_glob%vol%write(params_glob%vols(s), del_if_exists=.true.)
                if( present(which_iter) )then
                    ! post-process volume
                    pprocvol = add2fbody(trim(params_glob%vols(s)), params_glob%ext, PPROC_SUFFIX)
                    call build_glob%vol%fft()
                    ! low-pass filter
                    call build_glob%vol%bp(0., params_glob%lp)
                    call build_glob%vol%ifft()
                    ! mask
                    call build_glob%vol%mask(params_glob%msk, 'soft')
                    call build_glob%vol%write(pprocvol)
                endif
            endif
        end do
        call build_glob%vol%kill
    end subroutine norm_struct_facts

    subroutine eonorm_struct_facts( cline, which_iter )
        use simple_filterer, only: gen_anisotropic_optlp
        class(cmdline),    intent(inout) :: cline
        integer, optional, intent(in)    :: which_iter
        integer               :: s, find4eoavg
        real                  :: res05s(params_glob%nstates), res0143s(params_glob%nstates)
        character(len=STDLEN) :: pprocvol
        character(len=32)     :: resmskname
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
                call build_glob%eorecvols(s)%write_eos('recvol_state'//int2str_pad(s,2)//'_part'//int2str_pad(params_glob%part,params_glob%numlen))
            else
                if( present(which_iter) )then
                    params_glob%vols(s) = 'recvol_state'//int2str_pad(s,2)//'_iter'//int2str_pad(which_iter,3)//params_glob%ext
                else
                    params_glob%vols(s) = 'startvol_state'//int2str_pad(s,2)//params_glob%ext
                endif
                params_glob%vols_even(s) = add2fbody(params_glob%vols(s), params_glob%ext, '_even')
                params_glob%vols_odd(s)  = add2fbody(params_glob%vols(s), params_glob%ext, '_odd')
                resmskname         = 'resmask'//params_glob%ext
                call build_glob%eorecvols(s)%sum_eos
                call build_glob%eorecvols(s)%sampl_dens_correct_eos(s, params_glob%vols_even(s), params_glob%vols_odd(s), resmskname, find4eoavg)
                call gen_projection_frcs(cline,  params_glob%vols_even(s), params_glob%vols_odd(s), resmskname, s, build_glob%projfrcs)
                call build_glob%projfrcs%write('frcs_state'//int2str_pad(s,2)//'.bin')
                call gen_anisotropic_optlp(build_glob%vol2, build_glob%projfrcs, build_glob%eulspace_red, s, params_glob%pgrp, params_glob%hpind_fsc, params_glob%l_phaseplate)
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
    end subroutine eonorm_struct_facts

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
                ! masking
                if( cline%defined('mskfile') )then
                    ! mask provided
                    call mskvol%read(resmskname)
                    call vol%zero_background
                    call vol%mul(mskvol)
                else
                    ! circular masking
                    if( params_glob%l_innermsk )then
                        call vol%mask(params_glob%msk, 'soft', inner=params_glob%inner, width=params_glob%width)
                    else
                        call vol%mask(params_glob%msk, 'soft')
                    endif
                endif
        end subroutine prepeovol

    end subroutine gen_projection_frcs

end module simple_strategy2D3D_common
