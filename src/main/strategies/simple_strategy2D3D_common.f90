!@descr: common routines used by the high-level strategy 2D and 3D matchers
module simple_strategy2D3D_common
use simple_pftc_srch_api
use simple_timer
use simple_builder,            only: build_glob
use simple_butterworth,        only: butterworth_filter
use simple_discrete_stack_io,  only: dstack_io
use simple_nanoparticle_utils, only: phasecorr_one_atom
use simple_opt_filter,         only: estimate_lplim
use simple_projector,          only: projector
use simple_strategy2D_utils,   only: calc_cavg_offset
implicit none

public :: prepimgbatch, killimgbatch, read_imgbatch, discrete_read_imgbatch
public :: set_bp_range, set_bp_range2D, sample_ptcls4update, sample_ptcls4fillin, prepimg4align, prepimg4align_bench
public :: prep2Dref, build_batch_particles, prepare_refs_sigmas_ptcls, calc_3Drec
public :: calc_2Dref_offset
private
#include "simple_local_flags.inc"

interface read_imgbatch
    module procedure read_imgbatch_1
    module procedure read_imgbatch_2
    module procedure read_imgbatch_3
end interface read_imgbatch

type(stack_io)  :: stkio_r

contains

    !>  \brief  prepares a batch of image
    subroutine prepimgbatch( batchsz, box )
        integer,           intent(in) :: batchsz
        integer, optional, intent(in) :: box
        integer :: currsz, ibatch, box_here
        logical :: doprep
        if( .not. allocated(build_glob%imgbatch) )then
            doprep = .true.
        else
            currsz = size(build_glob%imgbatch)
            if( batchsz > currsz )then
                call killimgbatch
                doprep = .true.
            else
                doprep = .false.
            endif
        endif
        if( doprep )then
            box_here = params_glob%box
            if( present(box) ) box_here = box
            allocate(build_glob%imgbatch(batchsz))
            !$omp parallel do default(shared) private(ibatch) schedule(static) proc_bind(close)
            do ibatch = 1,batchsz
                call build_glob%imgbatch(ibatch)%new([box_here, box_here, 1], params_glob%smpd, wthreads=.false.)
            end do
            !$omp end parallel do
        endif
    end subroutine prepimgbatch

    subroutine killimgbatch
        integer :: ibatch
        if( allocated(build_glob%imgbatch) )then
            do ibatch=1,size(build_glob%imgbatch)
                call build_glob%imgbatch(ibatch)%kill
            end do
            deallocate(build_glob%imgbatch)
        endif
    end subroutine killimgbatch

    subroutine read_imgbatch_1( fromptop )
        integer, intent(in) :: fromptop(2)
        type(string) :: stkname
        integer :: iptcl, ind_in_batch, ind_in_stk
        do iptcl=fromptop(1),fromptop(2)
            ind_in_batch = iptcl - fromptop(1) + 1
            call build_glob%spproj%get_stkname_and_ind(params_glob%oritype, iptcl, stkname, ind_in_stk)
            if( .not. stkio_r%stk_is_open() )then
                call stkio_r%open(stkname, params_glob%smpd, 'read')
            else if( .not. stkio_r%same_stk(stkname, [params_glob%box,params_glob%box,1]) )then
                call stkio_r%close
                call stkio_r%open(stkname, params_glob%smpd, 'read')
            endif
            call stkio_r%read(ind_in_stk, build_glob%imgbatch(ind_in_batch))
        end do
        call stkio_r%close
    end subroutine read_imgbatch_1

    subroutine read_imgbatch_2( n, pinds, batchlims )
        integer, intent(in)  :: n, pinds(n), batchlims(2)
        type(string) :: stkname
        integer :: ind_in_stk, i, ii
        do i=batchlims(1),batchlims(2)
            ii = i - batchlims(1) + 1
            call build_glob%spproj%get_stkname_and_ind(params_glob%oritype, pinds(i), stkname, ind_in_stk)
            if( .not. stkio_r%stk_is_open() )then
                call stkio_r%open(stkname, params_glob%smpd, 'read')
            else if( .not. stkio_r%same_stk(stkname, [params_glob%box,params_glob%box,1]) )then
                call stkio_r%close
                call stkio_r%open(stkname, params_glob%smpd, 'read')
            endif
            call stkio_r%read(ind_in_stk, build_glob%imgbatch(ii))
        end do
        call stkio_r%close
    end subroutine read_imgbatch_2

    subroutine read_imgbatch_3( iptcl, img )
        integer,      intent(in)    :: iptcl
        class(image), intent(inout) :: img
        type(string) :: stkname
        integer :: ind_in_stk
        call build_glob%spproj%get_stkname_and_ind(params_glob%oritype, iptcl, stkname, ind_in_stk)
        if( .not. stkio_r%stk_is_open() )then
            call stkio_r%open(stkname, params_glob%smpd, 'read')
        else if( .not. stkio_r%same_stk(stkname, [params_glob%box,params_glob%box,1]) )then
            call stkio_r%close
            call stkio_r%open(stkname, params_glob%smpd, 'read')
        endif
        call stkio_r%read(ind_in_stk, img)
        call stkio_r%close
    end subroutine read_imgbatch_3

    subroutine discrete_read_imgbatch( n, pinds, batchlims )
        integer, intent(in)  :: n, pinds(n), batchlims(2)
        type(dstack_io) :: dstkio_r
        type(string) :: stkname
        integer :: ind_in_stk, i, ii
        call dstkio_r%new(params_glob%smpd, params_glob%box)
        do i=batchlims(1),batchlims(2)
            ii = i - batchlims(1) + 1
            call build_glob%spproj%get_stkname_and_ind(params_glob%oritype, pinds(i), stkname, ind_in_stk)
            call dstkio_r%read(stkname, ind_in_stk, build_glob%imgbatch(ii))
        end do
        call dstkio_r%kill
    end subroutine discrete_read_imgbatch

    subroutine set_bp_range( cline )
        class(cmdline), intent(in) :: cline
        real, allocatable :: resarr(:), fsc_arr(:)
        real              :: fsc0143, fsc05
        real              :: mapres(params_glob%nstates)
        integer           :: s, loc(1), lp_ind, arr_sz, fsc_sz, nyqcrop_ind
        type(string)      :: fsc_fname
        logical           :: fsc_bin_exists(params_glob%nstates), all_fsc_bin_exist
        if( params_glob%l_lpset )then
            ! set Fourier index range
            params_glob%kfromto(2) = calc_fourier_index(params_glob%lp, params_glob%box, params_glob%smpd)
            if( cline%defined('lpstop') )then
                params_glob%kfromto(2) = min(params_glob%kfromto(2),&
                    &calc_fourier_index(params_glob%lpstop, params_glob%box, params_glob%smpd))
            endif
            ! FSC values are read anyway
            do s=1,params_glob%nstates
                fsc_fname = FSC_FBODY//int2str_pad(s,2)//BIN_EXT
                if( file_exists(fsc_fname) )then
                    fsc_arr = file2rarr(fsc_fname)
                    fsc_sz  = size(build_glob%fsc(s,:))
                    arr_sz  = size(fsc_arr)
                    if( fsc_sz == arr_sz )then
                        build_glob%fsc(s,:) = fsc_arr(:)
                    else if( fsc_sz > arr_sz )then
                        ! padding
                        build_glob%fsc(s,:arr_sz)   = fsc_arr(:)
                        build_glob%fsc(s,arr_sz+1:) = 0.
                    else
                        ! clipping
                        build_glob%fsc(s,:fsc_sz)   = fsc_arr(:fsc_sz)
                    endif
                    deallocate(fsc_arr)
                endif
            enddo
        else
            ! check all fsc_state*.bin exist
            all_fsc_bin_exist = .true.
            fsc_bin_exists    = .false.
            do s=1,params_glob%nstates
                fsc_fname = FSC_FBODY//int2str_pad(s,2)//BIN_EXT
                fsc_bin_exists( s ) = file_exists(fsc_fname)
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
                        fsc_fname = FSC_FBODY//int2str_pad(s,2)//BIN_EXT
                        fsc_arr = file2rarr(fsc_fname)
                        fsc_sz  = size(build_glob%fsc(s,:))
                        arr_sz  = size(fsc_arr)
                        if( fsc_sz == arr_sz )then
                            build_glob%fsc(s,:) = fsc_arr(:)
                        else if( fsc_sz > arr_sz )then
                            ! padding
                            build_glob%fsc(s,:arr_sz)   = fsc_arr(:)
                            build_glob%fsc(s,arr_sz+1:) = 0.
                        else
                            ! clipping
                            build_glob%fsc(s,:fsc_sz)   = fsc_arr(:fsc_sz)
                        endif
                        deallocate(fsc_arr)
                        call get_resolution(build_glob%fsc(s,:), resarr, fsc05, fsc0143)
                        mapres(s) = fsc0143
                    else
                        ! empty state
                        mapres(s)           = 0.
                        build_glob%fsc(s,:) = 0.
                    endif
                end do
                loc = minloc(mapres) ! best resolved
                if( params_glob%nstates == 1 )then
                    lp_ind = get_find_at_crit(build_glob%fsc(1,:), params_glob%lplim_crit, incrreslim=params_glob%l_incrreslim)
                else
                    lp_ind = get_find_at_crit(build_glob%fsc(loc(1),:), params_glob%lplim_crit)
                endif
                ! interpolation limit is NOT Nyqvist in correlation search
                params_glob%kfromto(2) = calc_fourier_index(resarr(lp_ind), params_glob%box, params_glob%smpd)
            else if( build_glob%spproj_field%isthere(params_glob%fromp,'lp') )then
                params_glob%kfromto(2) = calc_fourier_index(&
                    build_glob%spproj_field%get(params_glob%fromp,'lp'), params_glob%box, params_glob%smpd)
            else
                THROW_HARD('no method available for setting the low-pass limit. Need fsc file or lp find; set_bp_range')
            endif
            ! lpstop overrides any other method for setting the low-pass limit
            if( cline%defined('lpstop') )then
                params_glob%kfromto(2) = min(params_glob%kfromto(2), &
                    calc_fourier_index(params_glob%lpstop, params_glob%box, params_glob%smpd))
            endif
            ! re-set the low-pass limit
            params_glob%lp = calc_lowpass_lim(params_glob%kfromto(2), params_glob%box, params_glob%smpd)
        endif
        ! making sure the resolution limit does not exceed limits of box_crop
        nyqcrop_ind            = calc_fourier_index(2.*params_glob%smpd_crop, params_glob%box_crop, params_glob%smpd_crop)
        params_glob%kfromto(2) = min(params_glob%kfromto(2), nyqcrop_ind)
        params_glob%lp         = calc_lowpass_lim(params_glob%kfromto(2), params_glob%box, params_glob%smpd)
        ! update low-pass limit in project
        call build_glob%spproj_field%set_all2single('lp',params_glob%lp)
    end subroutine set_bp_range

    subroutine set_bp_range2D( cline, which_iter, frac_srch_space )
        class(cmdline), intent(inout) :: cline
        integer,        intent(in)    :: which_iter
        real,           intent(in)    :: frac_srch_space
        real    :: lplim
        integer :: lpstart_find
        params_glob%kfromto(1) = max(2,calc_fourier_index(params_glob%hp, params_glob%box, params_glob%smpd))
        if( params_glob%l_lpset )then
            lplim = params_glob%lp
            params_glob%kfromto(2) = calc_fourier_index(lplim, params_glob%box_crop, params_glob%smpd_crop)
        else
            if( trim(params_glob%stream).eq.'yes' )then
                if( file_exists(params_glob%frcs) )then
                    lplim = build_glob%clsfrcs%estimate_lp_for_align()
                else
                    lplim = params_glob%lplims2D(3)
                endif
                if( cline%defined('lpstop') ) lplim = max(lplim, params_glob%lpstop)
            else
                if( file_exists(params_glob%frcs) .and. which_iter >= LPLIM1ITERBOUND )then
                    lplim = build_glob%clsfrcs%estimate_lp_for_align()
                else
                    if( which_iter < LPLIM1ITERBOUND )then
                        lplim = params_glob%lplims2D(1)
                    else if( frac_srch_space >= FRAC_SH_LIM .and. which_iter > LPLIM3ITERBOUND )then
                        lplim = params_glob%lplims2D(3)
                    else
                        lplim = params_glob%lplims2D(2)
                    endif
                endif
            endif
            params_glob%kfromto(2) = calc_fourier_index(lplim, params_glob%box_crop, params_glob%smpd_crop)
            ! to avoid pathological cases, fall-back on lpstart
            lpstart_find = calc_fourier_index(params_glob%lpstart, params_glob%box_crop, params_glob%smpd_crop)
            if( lpstart_find > params_glob%kfromto(2) ) params_glob%kfromto(2) = lpstart_find
            lplim = calc_lowpass_lim(params_glob%kfromto(2), params_glob%box_crop, params_glob%smpd_crop)
        endif
        ! update low-pas limit in project
        call build_glob%spproj_field%set_all2single('lp',lplim)
    end subroutine set_bp_range2D

    subroutine sample_ptcls4update( pfromto, l_incr_sampl, nptcls2update, pinds )
        integer,              intent(in)    :: pfromto(2)
        logical,              intent(in)    :: l_incr_sampl
        integer,              intent(inout) :: nptcls2update
        integer, allocatable, intent(inout) :: pinds(:)
        type(class_sample),   allocatable   :: clssmp(:)
        type(string) :: fname
        fname = CLASS_SAMPLING_FILE
        if( params_glob%l_update_frac )then
            if( trim(params_glob%balance).eq.'yes' )then
                if( file_exists(fname) )then
                    call read_class_samples(clssmp, fname)
                else
                    THROW_HARD('File for class-biased sampling in fractional update: '//CLASS_SAMPLING_FILE//' does not exist!')
                endif
                ! balanced class sampling
                if( params_glob%l_frac_best )then
                    call build_glob%spproj_field%sample4update_class(clssmp, pfromto, params_glob%update_frac,&
                    nptcls2update, pinds, l_incr_sampl, params_glob%l_greedy_smpl, frac_best=params_glob%frac_best)
                else
                    call build_glob%spproj_field%sample4update_class(clssmp, pfromto, params_glob%update_frac,&
                    nptcls2update, pinds, l_incr_sampl, params_glob%l_greedy_smpl)
                endif
                call deallocate_class_samples(clssmp)
            else
                call build_glob%spproj_field%sample4update_cnt(pfromto, params_glob%update_frac,&
                    nptcls2update, pinds, l_incr_sampl)
            endif
        else
            ! we sample all state > 0
            call build_glob%spproj_field%sample4update_all(pfromto, nptcls2update, pinds, l_incr_sampl)
        endif
        call fname%kill
    end subroutine sample_ptcls4update

    subroutine sample_ptcls4fillin( pfromto, l_incr_sampl, nptcls2update, pinds )
        integer,              intent(in)    :: pfromto(2)
        logical,              intent(in)    :: l_incr_sampl
        integer,              intent(inout) :: nptcls2update
        integer, allocatable, intent(inout) :: pinds(:)
        call build_glob%spproj_field%sample4update_fillin(pfromto, nptcls2update, pinds, l_incr_sampl)
    end subroutine sample_ptcls4fillin

    !>  \brief  prepares one particle image for alignment
    !          serial routine
    subroutine prepimg4align( iptcl, img, img_out, img_out_pd )
        integer,      intent(in)    :: iptcl
        class(image), intent(inout) :: img
        class(image), intent(inout) :: img_out
        class(image), intent(inout) :: img_out_pd
        type(ctf)       :: tfun
        type(ctfparams) :: ctfparms
        real            :: shvec(2), crop_factor
        crop_factor = real(params_glob%box_crop) / real(params_glob%box)
        shvec(1)    = -build_glob%spproj_field%get(iptcl, 'x') * crop_factor
        shvec(2)    = -build_glob%spproj_field%get(iptcl, 'y') * crop_factor
        ! Phase-flipping
        ctfparms = build_glob%spproj%get_ctfparams(params_glob%oritype, iptcl)
        select case(ctfparms%ctfflag)
            case(CTFFLAG_NO, CTFFLAG_FLIP)
                ! fused noise normalization, FFT, clip & shift
                call img%norm_noise_fft_clip_shift(build_glob%lmsk, img_out, shvec)
            case(CTFFLAG_YES)
                ctfparms%smpd = ctfparms%smpd / crop_factor != smpd_crop
                tfun          = ctf(ctfparms%smpd, ctfparms%kv, ctfparms%cs, ctfparms%fraca)
                ! fused noise normalization, FFT, clip, shift & CTF flip
                call img%norm_noise_fft_clip_shift_ctf_flip(build_glob%lmsk, img_out, shvec, tfun, ctfparms)
            case DEFAULT
                THROW_HARD('unsupported CTF flag: '//int2str(ctfparms%ctfflag)//' prepimg4align')
        end select
        ! fused IFFT, mask, FFT
        call img_out%ifft_mask_pad_fft(params_glob%msk_crop, img_out_pd)
    end subroutine prepimg4align

    subroutine prepimg4align_bench( iptcl, img, img_out, img_out_pd, rt_prep1, rt_prep2, rt_prep )
        integer,              intent(in)    :: iptcl
        class(image),         intent(inout) :: img
        class(image),         intent(inout) :: img_out
        class(image),         intent(inout) :: img_out_pd
        real(timer_int_kind), intent(inout) :: rt_prep1, rt_prep2, rt_prep
        integer(timer_int_kind) :: t_prep1, t_prep2, t_prep
        type(ctf)       :: tfun
        type(ctfparams) :: ctfparms
        real            :: shvec(2), crop_factor
        t_prep1 = tic()
        t_prep   = t_prep1
        crop_factor = real(params_glob%box_crop) / real(params_glob%box)
        shvec(1)    = -build_glob%spproj_field%get(iptcl, 'x') * crop_factor
        shvec(2)    = -build_glob%spproj_field%get(iptcl, 'y') * crop_factor
        ! Phase-flipping
        ctfparms = build_glob%spproj%get_ctfparams(params_glob%oritype, iptcl)
        select case(ctfparms%ctfflag)
            case(CTFFLAG_NO, CTFFLAG_FLIP)
                ! fused noise normalization, FFT, clip & shift
                call img%norm_noise_fft_clip_shift(build_glob%lmsk, img_out, shvec)
            case(CTFFLAG_YES)
                ctfparms%smpd = ctfparms%smpd / crop_factor != smpd_crop
                tfun          = ctf(ctfparms%smpd, ctfparms%kv, ctfparms%cs, ctfparms%fraca)
                ! fused noise normalization, FFT, clip, shift & CTF flip
                call img%norm_noise_fft_clip_shift_ctf_flip(build_glob%lmsk, img_out, shvec, tfun, ctfparms)
            case DEFAULT
                THROW_HARD('unsupported CTF flag: '//int2str(ctfparms%ctfflag)//' prepimg4align')
        end select
        rt_prep1 = rt_prep1 + toc(t_prep1)
        t_prep2 = tic()
        ! fused IFFT, mask, FFT
        call img_out%ifft_mask_pad_fft(params_glob%msk_crop, img_out_pd)
        rt_prep2 = rt_prep2 + toc(t_prep2)
        rt_prep   = rt_prep   + toc(t_prep)
    end subroutine prepimg4align_bench

    !>  \brief  Calculates the centering offset of the input cavg
    !>          cavg & particles centering is not performed
    subroutine calc_2Dref_offset( img, icls, centype, xyz )
        class(image), intent(inout) :: img
        integer,      intent(in)    :: icls
        integer,      intent(in)    :: centype
        real,         intent(out)   :: xyz(3)
        real :: xy_cavg(2), crop_factor
        crop_factor = real(params_glob%box_crop) / real(params_glob%box)
        select case(centype)
            case(PARAMS_CEN)
                call build_glob%spproj_field%calc_avg_offset2D(icls, xy_cavg)
                if( arg(xy_cavg) < CENTHRESH )then
                    xyz = 0.
                else if( arg(xy_cavg) > MAXCENTHRESH2D )then
                    xyz(1:2) = xy_cavg * crop_factor
                    xyz(3)   = 0.
                else
                    xyz = img%calc_shiftcen_serial(params_glob%cenlp, params_glob%msk_crop)
                    if( arg(xyz(1:2)/crop_factor - xy_cavg) > MAXCENTHRESH2D ) xyz = 0.
                endif
            case(SEG_CEN)
                call calc_cavg_offset(img, params_glob%cenlp, params_glob%msk_crop, xy_cavg)
                xyz = [xy_cavg(1), xy_cavg(2), 0.]
            case(MASS_CEN)
                xyz = img%calc_shiftcen_serial(params_glob%cenlp, params_glob%msk_crop)
        end select
        if( arg(xyz) < CENTHRESH ) xyz = 0.0
    end subroutine calc_2Dref_offset

    !>  \brief  Prepares one cluster centre image for alignment
    subroutine prep2Dref( img, icls, xyz, img_pd )
        class(image), intent(inout) :: img, img_pd
        integer,      intent(in)    :: icls
        real,         intent(in)    :: xyz(3)
        integer :: filtsz
        real    :: frc(img%get_filtsz()), filter(img%get_filtsz()), crop_factor
        ! Cavg & particles centering
        if( arg(xyz) > CENTHRESH )then
            call img%fft()
            call img%shift2Dserial(xyz(1:2))
            crop_factor = real(params_glob%box_crop) / real(params_glob%box)
            call build_glob%spproj_field%add_shift2class(icls, -xyz(1:2) / crop_factor)
        endif
        ! Filtering
        if( params_glob%l_ml_reg )then
            ! no filtering
        else
            if( params_glob%l_lpset.and.params_glob%l_icm )then
                ! ICM filter only applied when lp is set and performed below, FRC filtering turned off
            else
                ! FRC-based filtering
                call build_glob%clsfrcs%frc_getter(icls, frc)
                if( any(frc > 0.143) )then
                    filtsz = img%get_filtsz()
                    call fsc2optlp_sub(filtsz, frc, filter, merged=params_glob%l_lpset)
                    call img%fft() ! needs to be here in case the shift was never applied (above)
                    call img%apply_filter_serial(filter)
                endif
            endif
        endif
        ! ensure we are in real-space
        call img%ifft()
        ! ICM filter
        if( params_glob%l_lpset.and.params_glob%l_icm )then
            call img%ICM2D( params_glob%lambda, verbose=.false. )
        endif
        ! mask, pad, FFT
        call img%mask_pad_fft(params_glob%msk_crop, img_pd)
    end subroutine prep2Dref

    !>  \brief  initializes all volumes for reconstruction
    subroutine preprecvols
        integer, allocatable :: pops(:)
        integer :: istate
        call build_glob%spproj_field%get_pops(pops, 'state', maxn=params_glob%nstates)
        do istate = 1, params_glob%nstates
            if( pops(istate) > 0)then
                call build_glob%eorecvols(istate)%new(build_glob%spproj)
                call build_glob%eorecvols(istate)%reset_all
            endif
        end do
        deallocate(pops)
    end subroutine preprecvols

    !>  \brief  destructs all volumes for reconstruction
    subroutine killrecvols
        integer :: istate
        do istate = 1, params_glob%nstates
            call build_glob%eorecvols(istate)%kill
        end do
    end subroutine killrecvols

    !>  \brief  determines the reference volume shift and map shifts back to particles
    !>          reference volume shifting is performed in shift_and_mask_refvol
    subroutine calcrefvolshift_and_mapshifts2ptcls(cline, s, volfname, do_center, xyz, map_shift )
        class(cmdline), intent(in)  :: cline
        integer,        intent(in)  :: s
        class(string),  intent(in)  :: volfname
        logical,        intent(out) :: do_center
        real,           intent(out) :: xyz(3)
        logical,        intent(in)  :: map_shift
        real    :: crop_factor
        logical :: has_been_searched
        do_center = .true.
        ! centering
        if( params_glob%center .eq. 'no' .or. params_glob%nstates > 1 .or. &
            .not. params_glob%l_doshift .or. params_glob%pgrp(:1) .ne. 'c' .or. &
            params_glob%l_filemsk .or. params_glob%l_update_frac )then
            do_center = .false.
            xyz       = 0.
            return
        endif
        ! taking care of volume dimensions
        call build_glob%vol%read_and_crop(volfname, params_glob%smpd, params_glob%box_crop, params_glob%smpd_crop)
        ! offset
        xyz = build_glob%vol%calc_shiftcen(params_glob%cenlp,params_glob%msk_crop)
        if( params_glob%pgrp .ne. 'c1' ) xyz(1:2) = 0.     ! shifts only along z-axis for C2 and above
        if( arg(xyz) <= CENTHRESH )then
            do_center = .false.
            xyz       = 0.
            return
        endif
        if( map_shift )then
            ! map back to particle oritentations
            has_been_searched = .not.build_glob%spproj%is_virgin_field(params_glob%oritype)
            if( has_been_searched )then
                crop_factor = real(params_glob%box) / real(params_glob%box_crop)
                call build_glob%spproj_field%map3dshift22d(-xyz(:)*crop_factor, state=s)
            endif
        endif
    end subroutine calcrefvolshift_and_mapshifts2ptcls

    subroutine estimate_lp_refvols( lpfromto )
        real, optional, intent(in)    :: lpfromto(2)
        type(string)       :: vol_even, vol_odd, tmp
        real, allocatable  :: res(:)
        logical, parameter :: DEBUG=.false.
        type(image) :: mskvol
        integer     :: npix, s, loc
        real        :: lpest(params_glob%nstates), lpopt, res_fsc05, res_fsc0143, llpfromto(2)
        601 format(A,1X,F12.3)
        602 format(A,1X,F12.3,1X,F12.3)
        if( params_glob%l_filemsk )then
            ! read 3D envelope mask
            call mskvol%read_and_crop(params_glob%mskfile, params_glob%smpd, params_glob%box_crop, params_glob%smpd_crop)
            call mskvol%remove_edge
        else
            ! finding optimal lp over all states
            call mskvol%disc([params_glob%box_crop,  params_glob%box_crop, params_glob%box_crop],&
                         &params_glob%smpd_crop, params_glob%msk_crop, npix )
        endif
        ! find best resolved state
        res = get_resarr(params_glob%box_crop, params_glob%smpd_crop)
        do s = 1, params_glob%nstates
            call get_resolution_at_fsc(build_glob%fsc(s,:), res, 0.5, lpest(s))
        end do
        loc = minloc(lpest,dim=1)
        ! set range for low-pass limit estimation
        if( present(lpfromto) )then
            llpfromto = lpfromto
        else
            call get_resolution_at_fsc(build_glob%fsc(loc,:), res, 0.6, llpfromto(1))
            call get_resolution_at_fsc(build_glob%fsc(loc,:), res, 0.1, llpfromto(2))
        endif
        ! read volumes
        vol_even = params_glob%vols_even(loc)
        vol_odd  = params_glob%vols_odd(loc)
        if( params_glob%l_ml_reg )then
            ! estimate low-pass limit from unfiltered volumes
            tmp = add2fbody(vol_even, params_glob%ext,'_unfil')
            if( file_exists(tmp) ) vol_even = tmp
            tmp = add2fbody(vol_odd, params_glob%ext,'_unfil')
            if( file_exists(tmp) ) vol_even = tmp
        endif
        call build_glob%vol%read_and_crop(    vol_even, params_glob%smpd, params_glob%box_crop, params_glob%smpd_crop)
        call build_glob%vol_odd%read_and_crop(vol_odd,  params_glob%smpd, params_glob%box_crop, params_glob%smpd_crop)
        ! estimate low-pass limit
        call estimate_lplim(build_glob%vol_odd, build_glob%vol, mskvol, llpfromto, lpopt)
        ! generate output
        call get_resolution_at_fsc(build_glob%fsc(loc,:), res, 0.50,  res_fsc05)
        call get_resolution_at_fsc(build_glob%fsc(loc,:), res, 0.143, res_fsc0143)
        call build_glob%spproj_field%set_all2single('res', res_fsc0143)
        if( params_glob%part == 1 .and. DEBUG )then
        write(logfhandle,601) '>>> RESOLUTION @ FSC=0.5:                     ', res_fsc05
        write(logfhandle,601) '>>> RESOLUTION @ FSC=0.143:                   ', res_fsc0143
        write(logfhandle,602) '>>> LOW-PASS LIMIT RANGE FOR OPTIMIZATION:    ', llpfromto(1), llpfromto(2)
        write(logfhandle,601) '>>> OPTIMAL LOW-PASS LIMIT:                   ', lpopt
        endif
        ! update low-pass limit in project
        call build_glob%spproj_field%set_all2single('lp_est', lpopt)
        if( params_glob%l_lpauto )then
            ! re-set the low-pass limit
            params_glob%lp = lpopt
            ! update the Fourier index limit
            params_glob%kfromto(2) = calc_fourier_index(params_glob%lp, params_glob%box_crop, params_glob%smpd_crop)
            ! update low-pass limit in project
            call build_glob%spproj_field%set_all2single('lp',params_glob%lp)
        endif
        ! destruct
        call mskvol%kill
    end subroutine estimate_lp_refvols

    subroutine read_mask_and_filter_refvols( s )
        integer,  intent(in) :: s
        type(string)         :: vol_even, vol_odd, vol_avg
        logical, allocatable :: l_msk(:,:,:)
        real    :: cur_fil(params_glob%box_crop)
        integer :: filtsz
        ! READ
        vol_even = params_glob%vols_even(s)
        vol_odd  = params_glob%vols_odd(s)
        vol_avg  = params_glob%vols(s)
        call build_glob%vol%read_and_crop(   vol_even, params_glob%smpd, params_glob%box_crop, params_glob%smpd_crop)
        call build_glob%vol_odd%read_and_crop(vol_odd, params_glob%smpd, params_glob%box_crop, params_glob%smpd_crop)
        if( s == 1 .and. params_glob%l_filemsk )then
            ! read 3D envelope mask
            call build_glob%mskvol%read_and_crop(params_glob%mskfile, params_glob%smpd, params_glob%box_crop, params_glob%smpd_crop)
        endif
        ! noise regularization
        if( params_glob%l_noise_reg )then
            call build_glob%vol%add_gauran(params_glob%eps)
            call build_glob%vol_odd%add_gauran(params_glob%eps)
        endif
        ! MASK
        if( params_glob%l_filemsk )then
            ! envelope masking
            call build_glob%vol%zero_env_background(build_glob%mskvol)
            call build_glob%vol_odd%zero_env_background(build_glob%mskvol)
            call build_glob%vol%mul(build_glob%mskvol)
            call build_glob%vol_odd%mul(build_glob%mskvol)
        else
            ! circular masking
            call build_glob%vol%mask3D_soft(params_glob%msk_crop, backgr=0.0)
            call build_glob%vol_odd%mask3D_soft(params_glob%msk_crop, backgr=0.0)
        endif
        ! FILTER
        if( params_glob%l_icm )then
            if( params_glob%l_filemsk )then
                l_msk = build_glob%mskvol%bin2logical()
                call build_glob%vol%ICM3D_eo(build_glob%vol_odd, params_glob%lambda, l_msk)
            else
                call build_glob%vol%ICM3D_eo(build_glob%vol_odd, params_glob%lambda)
            endif
            if( params_glob%l_lpset .or. trim(params_glob%combine_eo).eq.'yes' )then ! no independent volume registration, so average eo pairs
                call build_glob%vol%add(build_glob%vol_odd)
                call build_glob%vol%mul(0.5)
                call build_glob%vol_odd%copy(build_glob%vol)
            endif
        else if( params_glob%l_lpset )then
            ! read average volume that will occupy both even and odd
            call build_glob%vol%read_and_crop(vol_avg, params_glob%smpd, params_glob%box_crop, params_glob%smpd_crop)
            ! noise regularization
            if( params_glob%l_noise_reg )then
                call build_glob%vol%add_gauran(params_glob%eps)
            endif
            ! mask again, BP filter performed below
            if( params_glob%l_filemsk )then
                ! envelope masking
                call build_glob%vol%zero_env_background(build_glob%mskvol)
                call build_glob%vol%mul(build_glob%mskvol)
            else
                ! circular masking
                call build_glob%vol%mask3D_soft(params_glob%msk_crop, backgr=0.0)
            endif
            ! odd <- even
            call build_glob%vol_odd%copy(build_glob%vol)
        endif
        call build_glob%vol%fft
        call build_glob%vol_odd%fft
        if( params_glob%l_ml_reg )then
            ! filtering done when volumes are assembled
        else if( params_glob%l_icm )then
            ! filtering done above
        else if( params_glob%l_lpset )then
            if( trim(params_glob%gauref)=='yes' )then
                ! Gaussian filter
                call build_glob%vol%bpgau3D(0., params_glob%gaufreq)
                call build_glob%vol_odd%bpgau3D(0., params_glob%gaufreq)
            else
                ! Cosine low-pass filter, works best for nanoparticles
                call build_glob%vol%bp(0., params_glob%lp)
                call build_glob%vol_odd%bp(0., params_glob%lp)
            endif
        else
            ! Optimal filter
            filtsz = build_glob%vol%get_filtsz()
            if( any(build_glob%fsc(s,:) > 0.143) )then
                call fsc2optlp_sub(filtsz,build_glob%fsc(s,:),cur_fil)
                call build_glob%vol%apply_filter(cur_fil)
                call build_glob%vol_odd%apply_filter(cur_fil)
            endif
        endif
    end subroutine read_mask_and_filter_refvols

    !>  \brief  grids one particle image to the volume
    subroutine grid_ptcl( fpl, se, o )
        class(fplane_type), intent(in)    :: fpl
        class(sym),         intent(inout) :: se
        class(ori),         intent(inout) :: o
        real    :: pw
        integer :: s, eo
        ! state flag
        s = o%get_state()
        if( s == 0 ) return
        ! eo flag
        eo = o%get_eo()
        ! particle-weight
        pw = 1.0
        if( o%isthere('w') ) pw = o%get('w')
        if( pw > TINY ) call build_glob%eorecvols(s)%grid_plane(se, o, fpl, eo, pw)
    end subroutine grid_ptcl

    !> volumetric 3d reconstruction
    subroutine calc_3Drec( cline, nptcls2update, pinds )
        use simple_imgarr_utils, only: alloc_imgarr, dealloc_imgarr
        class(cmdline),    intent(inout) :: cline
        integer,           intent(in)    :: nptcls2update
        integer,           intent(in)    :: pinds(nptcls2update)
        type(fplane_type), allocatable   :: fpls(:)
        type(ctfparams),   allocatable   :: ctfparms(:)
        type(ori) :: orientation, o_thres, o_mirr
        real      :: shift(2), euls(3), euls_mirr(3)
        integer   :: batchlims(2), iptcl, i, i_batch, ibatch, nptcls_eff, ithr
        logical   :: l_rec_state
        logical   :: DEBUG = .false.
        integer(timer_int_kind) :: t
        real(timer_int_kind)    :: t_ini, t_2dprep, t_grid
        if( trim(params_glob%recthres).eq.'yes' )then
            euls = [params_glob%e1, params_glob%e2, params_glob%e3]
            call euler_mirror(euls, euls_mirr)
            call o_thres%new(is_ptcl=.true.)
            call o_thres%set_euler(euls)
            call o_mirr%new(is_ptcl=.true.)
            call o_mirr%set_euler(euls_mirr)
        endif
        l_rec_state = cline%defined('state')
        ! init volumes
        call preprecvols
        ! prep batch imgs
        call prepimgbatch(MAXIMGBATCHSZ)
        ! allocate array
        allocate(fpls(MAXIMGBATCHSZ),ctfparms(MAXIMGBATCHSZ))
        if( DEBUG )then
            t_ini    = 0.
            t_2dprep = 0.
            t_grid   = 0.
        endif
        ! logical/physical adress mapping
        call memoize_ft_maps([params_glob%boxpd, params_glob%boxpd, 1], build_glob%imgbatch(1)%get_smpd())
        ! heap of padded images
        call alloc_imgarr(nthr_glob, [params_glob%boxpd, params_glob%boxpd, 1], build_glob%imgbatch(1)%get_smpd(), build_glob%img_pad_heap)
        ! gridding batch loop
        nptcls_eff = 0
        do i_batch=1,nptcls2update,MAXIMGBATCHSZ
            if( DEBUG ) t = tic()
            batchlims = [i_batch,min(nptcls2update,i_batch + MAXIMGBATCHSZ - 1)]
            call discrete_read_imgbatch( nptcls2update, pinds, batchlims)
            if( DEBUG )then
                t_ini = t_ini + toc(t)
                t = tic()
            endif
            !$omp parallel do default(shared) private(i,ithr,iptcl,ibatch,shift) schedule(static) proc_bind(close)
            do i=batchlims(1),batchlims(2)
                ithr   = omp_get_thread_num() + 1
                iptcl  = pinds(i)
                ibatch = i - batchlims(1) + 1
                call build_glob%imgbatch(ibatch)%norm_noise_taper_edge_pad_fft(build_glob%lmsk, build_glob%img_pad_heap(ithr))
                ctfparms(ibatch) = build_glob%spproj%get_ctfparams(params_glob%oritype, iptcl)
                shift = build_glob%spproj_field%get_2Dshift(iptcl)
                call build_glob%img_pad_heap(ithr)%gen_fplane4rec(ctfparms(ibatch), shift, iptcl, fpls(ibatch))
            end do
            !$omp end parallel do
            if( DEBUG )then
                t_2dprep = t_2dprep + toc(t)
                t = tic()
            endif
            ! gridding
            do i=batchlims(1),batchlims(2)
                iptcl  = pinds(i)
                ibatch = i - batchlims(1) + 1
                call build_glob%spproj_field%get_ori(iptcl, orientation)
                if( orientation%isstatezero() ) cycle
                if( trim(params_glob%recthres).eq.'yes' )then
                    if( rad2deg(orientation .euldist. o_thres) < params_glob%rec_athres .or.&
                       &rad2deg(orientation .euldist. o_mirr ) < params_glob%rec_athres ) cycle
                endif
                nptcls_eff = nptcls_eff + 1
                call grid_ptcl(fpls(ibatch), build_glob%pgrpsyms, orientation)
            end do
            if( DEBUG )then
                t_grid = t_grid + toc(t)
            endif
        end do
        call dealloc_imgarr(build_glob%img_pad_heap)
        if( DEBUG ) print *,'timing: ',t_ini, t_2dprep, t_grid
        if( trim(params_glob%recthres).eq.'yes' ) print *, 'nptcls in 3D reconstruction = ', nptcls_eff
        ! normalise structure factors
        call norm_struct_facts( cline )
        ! destruct
        call forget_ft_maps
        call killrecvols()
        do ibatch=1,MAXIMGBATCHSZ
            if( allocated(fpls(ibatch)%cmplx_plane) ) deallocate(fpls(ibatch)%cmplx_plane)
            if( allocated(fpls(ibatch)%ctfsq_plane) ) deallocate(fpls(ibatch)%ctfsq_plane)
            ! call fpls(ibatch)%kill
        end do
        deallocate(fpls,ctfparms)
        call orientation%kill
    end subroutine calc_3Drec

    subroutine norm_struct_facts( cline )
        use simple_gridding, only: prep3D_inv_instrfun4mul
        class(cmdline),    intent(inout) :: cline
        type(string) :: recname, volname, volname_prev, volname_prev_even
        type(string) :: volname_prev_odd, str_state, str_iter, fsc_txt_file, eonames(2)
        real, allocatable :: fsc(:)
        type(image) :: vol_prev_even, vol_prev_odd, gridcorr_img
        integer     :: s, find4eoavg, ldim_pd(3), ldim(3)
        real        :: res05s(params_glob%nstates), res0143s(params_glob%nstates)
        real        :: weight_prev, update_frac_trail_rec
        ! init
        ldim    = [params_glob%box_crop,  params_glob%box_crop,  params_glob%box_crop]
        ldim_pd = [params_glob%box_croppd,params_glob%box_croppd,params_glob%box_croppd]
        call build_glob%vol%new(ldim,params_glob%smpd_crop)
        call build_glob%vol2%new(ldim,params_glob%smpd_crop)
        res0143s = 0.
        res05s   = 0.
        ! read in previous reconstruction when trail_rec==yes
        update_frac_trail_rec = 1.0
        if( .not. params_glob%l_distr_exec .and. params_glob%l_trail_rec )then
            if( cline%defined('ufrac_trec') )then
                update_frac_trail_rec = params_glob%ufrac_trec
            else
                update_frac_trail_rec = build_glob%spproj_field%get_update_frac()
            endif
        endif
        ! Prep for correction of the shape of the interpolator
        gridcorr_img = prep3D_inv_instrfun4mul(ldim, ldim_pd, params_glob%smpd_crop)
        ! cycle through states
        do s=1,params_glob%nstates
            if( build_glob%spproj_field%get_pop(s, 'state') == 0 )then
                ! empty state
                build_glob%fsc(s,:) = 0.
                cycle
            endif
            call build_glob%eorecvols(s)%compress_exp
            if( params_glob%l_distr_exec )then
                call build_glob%eorecvols(s)%write_eos(string(VOL_FBODY)//int2str_pad(s,2)//'_part'//int2str_pad(params_glob%part,params_glob%numlen))
            else
                ! global volume name update
                recname = VOL_FBODY//int2str_pad(s,2)
                volname = recname//params_glob%ext
                if( params_glob%l_filemsk .and. params_glob%l_envfsc )then
                    call build_glob%eorecvols(s)%set_automsk(.true.)
                endif
                eonames(1) = recname//'_even'//params_glob%ext%to_char()
                eonames(2) = recname//'_odd'//params_glob%ext%to_char()
                if( params_glob%l_ml_reg )then
                    ! the sum is done after regularization
                else
                    call build_glob%eorecvols(s)%sum_eos
                endif
                if( params_glob%l_trail_rec )then
                    if( .not. cline%defined('vol'//int2str(s)) ) THROW_HARD('vol'//int2str(s)//'required in norm_struct_facts cline when trail_rec==yes')
                    volname_prev      = cline%get_carg('vol'//int2str(s))
                    volname_prev_even = add2fbody(volname_prev, params_glob%ext, '_even')
                    volname_prev_odd  = add2fbody(volname_prev, params_glob%ext, '_odd')
                    if( .not. file_exists(volname_prev_even) ) THROW_HARD('File: '//volname_prev_even%to_char()//' does not exist!')
                    if( .not. file_exists(volname_prev_odd)  ) THROW_HARD('File: '//volname_prev_odd%to_char()//' does not exist!')
                    call vol_prev_even%read_and_crop(volname_prev_even, params_glob%smpd, params_glob%box_crop, params_glob%smpd_crop)
                    call vol_prev_odd %read_and_crop(volname_prev_odd,  params_glob%smpd, params_glob%box_crop, params_glob%smpd_crop)
                    if( allocated(fsc) ) deallocate(fsc)
                    call build_glob%eorecvols(s)%calc_fsc4sampl_dens_correct(vol_prev_even, vol_prev_odd, fsc)
                    call build_glob%eorecvols(s)%sampl_dens_correct_eos(s, eonames(1), eonames(2), find4eoavg, fsc)
                else 
                    call build_glob%eorecvols(s)%sampl_dens_correct_eos(s, eonames(1), eonames(2), find4eoavg)
                endif
                str_state = int2str_pad(s,2)
                if( cline%defined('which_iter') )then
                    str_iter     = int2str_pad(params_glob%which_iter,3)
                    fsc_txt_file = 'RESOLUTION_STATE'//str_state%to_char()//'_ITER'//str_iter%to_char()
                else
                    fsc_txt_file = 'RESOLUTION_STATE'//str_state%to_char()
                endif
                call build_glob%eorecvols(s)%write_fsc2txt(fsc_txt_file)
                if( params_glob%l_ml_reg )then
                    call build_glob%eorecvols(s)%sum_eos
                endif
                call build_glob%eorecvols(s)%get_res(res05s(s), res0143s(s))
                call build_glob%eorecvols(s)%sampl_dens_correct_sum(build_glob%vol)
                ! need to put the sum back at lowres for the eo pairs
                call build_glob%vol%fft
                call build_glob%vol2%zero_and_unflag_ft
                call build_glob%vol2%read(eonames(1))
                call build_glob%vol2%fft()
                call build_glob%vol2%insert_lowres(build_glob%vol, find4eoavg)
                call build_glob%vol2%ifft()
                call build_glob%vol2%mul(gridcorr_img)
                call build_glob%vol2%write(eonames(1), del_if_exists=.true.)
                call build_glob%vol2%zero_and_unflag_ft
                call build_glob%vol2%read(eonames(2))
                call build_glob%vol2%fft()
                call build_glob%vol2%insert_lowres(build_glob%vol, find4eoavg)
                call build_glob%vol2%ifft()
                call build_glob%vol2%mul(gridcorr_img)
                call build_glob%vol2%write(eonames(2), del_if_exists=.true.)
                ! merged volume
                call build_glob%vol%ifft
                call build_glob%vol%mul(gridcorr_img)
                call build_glob%vol%write(volname, del_if_exists=.true.)
                if( params_glob%l_trail_rec .and. update_frac_trail_rec < 0.99 )then
                    call build_glob%vol%read(eonames(1))  ! even current
                    call build_glob%vol2%read(eonames(2)) ! odd current
                    weight_prev = 1. - update_frac_trail_rec
                    call vol_prev_even%mul(weight_prev)
                    call vol_prev_odd%mul (weight_prev)
                    call build_glob%vol%mul(update_frac_trail_rec)
                    call build_glob%vol2%mul(update_frac_trail_rec)
                    call build_glob%vol%add(vol_prev_even)
                    call build_glob%vol2%add(vol_prev_odd)
                    call build_glob%vol%write(eonames(1))  ! even trailed
                    call build_glob%vol2%write(eonames(2)) ! odd trailed
                    call vol_prev_even%kill
                    call vol_prev_odd%kill
                endif
                call build_glob%vol%fft()
                call build_glob%vol2%fft()
                ! updating command-line accordingly (needed in multi-stage wflows)
                call cline%set('vol'//int2str(s), volname)
                ! updating the global parameter object accordingly (needed in multi-stage wflows)
                params_glob%vols(s) = volname
            endif
        end do
        call build_glob%vol2%kill
        call gridcorr_img%kill
    end subroutine norm_struct_facts

    subroutine prepare_refs_sigmas_ptcls( pftc, cline, eucl_sigma, ptcl_imgs, ptcl_imgs_pad, batchsz, which_iter, do_polar )
        class(polarft_calc),      intent(inout) :: pftc
        class(cmdline),           intent(in)    :: cline !< command line
        class(euclid_sigma2),     intent(inout) :: eucl_sigma
        type(image), allocatable, intent(inout) :: ptcl_imgs(:)
        type(image), allocatable, intent(inout) :: ptcl_imgs_pad(:)
        integer,                  intent(in)    :: batchsz
        integer,                  intent(in)    :: which_iter
        logical,     optional,    intent(in)    :: do_polar
        real, allocatable :: gaufilter(:)
        type(string)      :: fname
        integer           :: ithr, iproj, nrefs, filtsz
        logical           :: l_polar, l_filtrefs
        l_polar = .false.
        if( present(do_polar) ) l_polar = do_polar
        ! PREPARATION OF pftc AND REFERENCES
        if( l_polar )then
            ! PREPARATION OF pftc AND REFERENCES
            nrefs = params_glob%nspace * params_glob%nstates
            call pftc%new(nrefs, [1,batchsz], params_glob%kfromto)
            ! Read polar references
            call pftc%polar_cavger_new(.true.)
            call pftc%polar_cavger_read_all(string(POLAR_REFS_FBODY//BIN_EXT))
            call build_glob%clsfrcs%read(string(FRCS_FILE))
            ! prepare filter
            l_filtrefs = .false.
            if(trim(params_glob%gauref).eq.'yes')then
                l_filtrefs = .true.
                filtsz = build_glob%clsfrcs%get_filtsz()
                allocate(gaufilter(filtsz),source=0.)
                call gaussian_filter(params_glob%gaufreq, params_glob%smpd, params_glob%box, gaufilter)
            endif
            ! PREPARATION OF REFERENCES IN pftc
            !$omp parallel do default(shared) private(iproj)&
            !$omp schedule(static) proc_bind(close)
            do iproj = 1,params_glob%nspace
                if( l_filtrefs ) call pftc%polar_filterrefs(iproj, gaufilter)
                ! transfer to pftc
                if( params_glob%l_lpset )then
                    ! put the merged class average in both even and odd positions
                    call pftc%polar_cavger_set_ref_pft(iproj, 'merged')
                    call pftc%cp_even2odd_ref(iproj)
                else
                    ! transfer e/o refs to pftc
                    call pftc%polar_cavger_set_ref_pft(iproj, 'even')
                    call pftc%polar_cavger_set_ref_pft(iproj, 'odd')
                endif
            end do
            !$omp end parallel do
            ! Center prior to memoization
            if( (trim(params_glob%center)=='yes') .and. (trim(params_glob%center_type)=='params') .and.&
                &(params_glob%pgrp(:1)=='c1') .and. (.not.params_glob%l_update_frac) .and.&
                &(params_glob%nstates==1) .and. params_glob%l_doshift )then
                call pftc%center_3Dpolar_refs(build_glob%spproj_field, build_glob%eulspace)
            endif
            ! Memoize
            call pftc%memoize_refs
        else
            ! (if needed) estimating lp (over all states) and reseting params_glob%lp and params_glob%kfromto
            call prepare_polar_references(pftc, cline, batchsz)
        endif
        ! PREPARATION OF SIGMAS
        if( params_glob%l_needs_sigma )then
            fname = SIGMA2_FBODY//int2str_pad(params_glob%part,params_glob%numlen)//'.dat'
            call eucl_sigma%new(fname, params_glob%box)
            call eucl_sigma%read_part(  build_glob%spproj_field)
            call eucl_sigma%read_groups(build_glob%spproj_field)
        end if
        ! PREPARATION OF PARTICLES
        call prepimgbatch(batchsz)
        allocate(ptcl_imgs(nthr_glob), ptcl_imgs_pad(nthr_glob))
        !$omp parallel do default(shared) private(ithr) schedule(static) proc_bind(close)
        do ithr = 1,nthr_glob
            call ptcl_imgs(ithr)%new(    [params_glob%box_crop,  params_glob%box_crop,  1], params_glob%smpd_crop, wthreads=.false.)
            call ptcl_imgs_pad(ithr)%new([params_glob%box_croppd,params_glob%box_croppd,1], params_glob%smpd_crop, wthreads=.false.)
        enddo
        !$omp end parallel do
        call build_glob%vol%kill
        call build_glob%vol_odd%kill
        call build_glob%vol2%kill
        if( allocated(gaufilter) ) deallocate(gaufilter)
    end subroutine prepare_refs_sigmas_ptcls

    subroutine prepare_polar_references( pftc, cline, batchsz )
        class(polarft_calc), intent(inout) :: pftc
        class(cmdline),      intent(in)    :: cline !< command line
        integer,             intent(in)    :: batchsz
        type(ori) :: o_tmp
        real      :: xyz(3)
        integer   :: s, iproj, iref, nrefs
        logical   :: do_center
        if( cline%defined('lpstart') .and. cline%defined('lpstop') )then
            call estimate_lp_refvols([params_glob%lpstart,params_glob%lpstop])
        else
            call estimate_lp_refvols
        endif
        ! pftc
        nrefs = params_glob%nspace * params_glob%nstates
        call pftc%new(nrefs, [1,batchsz], params_glob%kfromto)
        ! read reference volumes and create polar projections
        do s=1,params_glob%nstates
            if( str_has_substr(params_glob%refine, 'prob') )then
                ! already mapping shifts in prob_tab with shared-memory execution
                call calcrefvolshift_and_mapshifts2ptcls( cline, s, params_glob%vols(s), do_center, xyz, map_shift=l_distr_exec_glob)
            else
                call calcrefvolshift_and_mapshifts2ptcls( cline, s, params_glob%vols(s), do_center, xyz, map_shift=.true.)
            endif
            call read_mask_and_filter_refvols(s)
            ! PREPARE E/O VOLUMES
            ! do even/odd in separate passes to reduce memory usage when padding is large (2X)

            ! PREPARE EVEN REFERENCES
            call build_glob%vol_pad%new([params_glob%box_croppd, params_glob%box_croppd, params_glob%box_croppd],&
            &params_glob%smpd_crop, wthreads=.true.)
            if( do_center )then
                call build_glob%vol%fft()
                call build_glob%vol%shift([xyz(1),xyz(2),xyz(3)])
            endif
            ! back to real space
            call build_glob%vol%ifft()
            ! FT volume
            call build_glob%vol%pad_fft(build_glob%vol_pad)
            ! expand for fast interpolation & correct for norm when clipped
            call build_glob%vol_pad%expand_cmat
            !$omp parallel do default(shared) private(iproj,o_tmp,iref) schedule(static) proc_bind(close)
            do iproj=1,params_glob%nspace
                iref = (s - 1) * params_glob%nspace + iproj
                call build_glob%eulspace%get_ori(iproj, o_tmp)
                call build_glob%vol_pad%fproject_polar_oversamp(iref, o_tmp, pftc, iseven=.true., mask=build_glob%l_resmsk)
                call o_tmp%kill
            end do
            !$omp end parallel do
            call build_glob%vol_pad%kill
            call build_glob%vol_pad%kill_expanded

            ! PREPARE ODD REFERENCES
            call build_glob%vol_odd_pad%new([params_glob%box_croppd, params_glob%box_croppd, params_glob%box_croppd],&
            &params_glob%smpd_crop, wthreads=.true.)
            if( do_center )then
                call build_glob%vol_odd%fft()
                call build_glob%vol_odd%shift([xyz(1),xyz(2),xyz(3)])
            endif
            ! back to real space
            call build_glob%vol_odd%ifft()
            ! FT volume
            call build_glob%vol_odd%pad_fft(build_glob%vol_odd_pad)
            ! expand for fast interpolation & correct for norm when clipped
            call build_glob%vol_odd_pad%expand_cmat
            !$omp parallel do default(shared) private(iproj,o_tmp,iref) schedule(static) proc_bind(close)
            do iproj=1,params_glob%nspace
                iref = (s - 1) * params_glob%nspace + iproj
                call build_glob%eulspace%get_ori(iproj, o_tmp)
                call build_glob%vol_odd_pad%fproject_polar_oversamp(iref, o_tmp, pftc, iseven=.false., mask=build_glob%l_resmsk)
                call o_tmp%kill
            end do
            !$omp end parallel do
            call build_glob%vol_odd_pad%kill
            call build_glob%vol_odd_pad%kill_expanded
        end do
        call pftc%memoize_refs
    end subroutine prepare_polar_references

    subroutine build_batch_particles( pftc, nptcls_here, pinds_here, tmp_imgs, tmp_imgs_pad )
        class(polarft_calc), intent(inout) :: pftc
        integer,             intent(in)    :: nptcls_here
        integer,             intent(in)    :: pinds_here(nptcls_here)
        type(image),         intent(inout) :: tmp_imgs(params_glob%nthr), tmp_imgs_pad(params_glob%nthr)
        complex, allocatable :: pft(:,:)
        integer     :: iptcl_batch, iptcl, ithr
        ! reassign particles indices & associated variables
        call pftc%reallocate_ptcls(nptcls_here, pinds_here)
        call discrete_read_imgbatch( nptcls_here, pinds_here, [1,nptcls_here])
        ! mask memoization for prepimg4align
        call tmp_imgs(1)%memoize_mask_coords
        ! memoize CTF stuff
        call memoize_ft_maps(tmp_imgs(1)%get_ldim(), tmp_imgs(1)%get_smpd())   
        ! memoize for polarize_oversamp
        call tmp_imgs_pad(1)%memoize4polarize_oversamp(pftc%get_pdim())     
        !$omp parallel do default(shared) private(iptcl,iptcl_batch,ithr,pft) schedule(static) proc_bind(close)
        do iptcl_batch = 1,nptcls_here
            ithr  = omp_get_thread_num() + 1
            iptcl = pinds_here(iptcl_batch)
            ! prep
            call prepimg4align(iptcl, build_glob%imgbatch(iptcl_batch), tmp_imgs(ithr), tmp_imgs_pad(ithr))
            ! transfer to polar coordinates
            pft = pftc%allocate_pft()
            call tmp_imgs_pad(ithr)%polarize_oversamp(pft, mask=build_glob%l_resmsk)
            call pftc%set_ptcl_pft(iptcl, pft)
            deallocate(pft)
            ! e/o flags
            call pftc%set_eo(iptcl, nint(build_glob%spproj_field%get(iptcl,'eo'))<=0 )
        end do
        !$omp end parallel do
        call pftc%create_polar_absctfmats(build_glob%spproj, 'ptcl3D')
        ! Memoize particles FFT parameters
        call pftc%memoize_ptcls
        ! destruct
        call forget_ft_maps
    end subroutine build_batch_particles

end module simple_strategy2D3D_common
