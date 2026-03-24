!@descr: common routines used by the high-level strategy 2D and 3D matchers
module simple_strategy2D3D_common
use simple_pftc_srch_api
use simple_timer
use simple_builder,            only: builder
use simple_butterworth,        only: butterworth_filter
use simple_discrete_stack_io,  only: dstack_io
use simple_opt_filter,         only: estimate_lplim3D
use simple_polarft_calc,       only: vol_pad2ref_pfts
use simple_projector,          only: projector
use simple_strategy2D_utils,   only: calc_cavg_offset
implicit none

! Particle image I/O
public :: prepimgbatch, killimgbatch, read_imgbatch, discrete_read_imgbatch
! Resolutions limits for alignment
public :: set_bp_range, set_bp_range2D
! Particle sampling
public :: sample_ptcls4update, sample_ptcls4fillin
! Particle image processing for alignment
public :: prepimg4align, prepimg4align_bench, build_batch_particles
! Reference processing for alignment
public :: prep2Dref, prepare_refs_sigmas_ptcls, calc_2Dref_offset
! Offline reconstruction
public :: calc_3Drec
! Online reconstruction
public :: init_rec, prep_imgs4rec, update_rec, finalize_rec
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
    subroutine prepimgbatch( params, build, batchsz, box )
        class(parameters), intent(in)    :: params
        class(builder),    intent(inout) :: build
        integer,           intent(in)    :: batchsz
        integer, optional, intent(in)    :: box
        integer :: currsz, ibatch, box_here
        logical :: doprep
        if( .not. allocated(build%imgbatch) )then
            doprep = .true.
        else
            currsz = size(build%imgbatch)
            if( batchsz > currsz )then
                call killimgbatch(build)
                doprep = .true.
            else
                doprep = .false.
            endif
        endif
        if( doprep )then
            box_here = params%box
            if( present(box) ) box_here = box
            allocate(build%imgbatch(batchsz))
            !$omp parallel do default(shared) private(ibatch) schedule(static) proc_bind(close)
            do ibatch = 1,batchsz
                call build%imgbatch(ibatch)%new([box_here, box_here, 1], params%smpd, wthreads=.false.)
            end do
            !$omp end parallel do
        endif
    end subroutine prepimgbatch

    subroutine killimgbatch( build )
        type(builder), intent(inout) :: build
        integer :: ibatch
        if( allocated(build%imgbatch) )then
            do ibatch=1,size(build%imgbatch)
                call build%imgbatch(ibatch)%kill
            end do
            deallocate(build%imgbatch)
        endif
    end subroutine killimgbatch

    subroutine read_imgbatch_1( params, build, fromptop )
        class(parameters), intent(in)    :: params
        class(builder),    intent(inout) :: build
        integer,           intent(in)    :: fromptop(2)
        type(string) :: stkname
        integer :: iptcl, ind_in_batch, ind_in_stk
        do iptcl=fromptop(1),fromptop(2)
            ind_in_batch = iptcl - fromptop(1) + 1
            call build%spproj%get_stkname_and_ind(params%oritype, iptcl, stkname, ind_in_stk)
            if( .not. stkio_r%stk_is_open() )then
                call stkio_r%open(stkname, params%smpd, 'read')
            else if( .not. stkio_r%same_stk(stkname, [params%box,params%box,1]) )then
                call stkio_r%close
                call stkio_r%open(stkname, params%smpd, 'read')
            endif
            call stkio_r%read(ind_in_stk, build%imgbatch(ind_in_batch))
        end do
        call stkio_r%close
    end subroutine read_imgbatch_1

    subroutine read_imgbatch_2( params, build, n, pinds, batchlims )
        class(parameters), intent(in)    :: params
        class(builder),    intent(inout) :: build
        integer,           intent(in)    :: n, pinds(n), batchlims(2)
        type(string) :: stkname
        integer :: ind_in_stk, i, ii
        do i=batchlims(1),batchlims(2)
            ii = i - batchlims(1) + 1
            call build%spproj%get_stkname_and_ind(params%oritype, pinds(i), stkname, ind_in_stk)
            if( .not. stkio_r%stk_is_open() )then
                call stkio_r%open(stkname, params%smpd, 'read')
            else if( .not. stkio_r%same_stk(stkname, [params%box,params%box,1]) )then
                call stkio_r%close
                call stkio_r%open(stkname, params%smpd, 'read')
            endif
            call stkio_r%read(ind_in_stk, build%imgbatch(ii))
        end do
        call stkio_r%close
    end subroutine read_imgbatch_2

    subroutine read_imgbatch_3( params, build, iptcl, img )
        class(parameters), intent(in)    :: params
        class(builder),    intent(inout) :: build
        integer,           intent(in)    :: iptcl
        class(image),      intent(inout) :: img
        type(string) :: stkname
        integer :: ind_in_stk
        call build%spproj%get_stkname_and_ind(params%oritype, iptcl, stkname, ind_in_stk)
        if( .not. stkio_r%stk_is_open() )then
            call stkio_r%open(stkname, params%smpd, 'read')
        else if( .not. stkio_r%same_stk(stkname, [params%box,params%box,1]) )then
            call stkio_r%close
            call stkio_r%open(stkname, params%smpd, 'read')
        endif
        call stkio_r%read(ind_in_stk, img)
        call stkio_r%close
    end subroutine read_imgbatch_3

    subroutine discrete_read_imgbatch( params, build, n, pinds, batchlims )
        class(parameters), intent(in)    :: params
        class(builder),    intent(inout) :: build
        integer,           intent(in)    :: n, pinds(n), batchlims(2)
        type(dstack_io) :: dstkio_r
        type(string) :: stkname
        integer :: ind_in_stk, i, ii
        call dstkio_r%new(params%smpd, params%box)
        do i=batchlims(1),batchlims(2)
            ii = i - batchlims(1) + 1
            call build%spproj%get_stkname_and_ind(params%oritype, pinds(i), stkname, ind_in_stk)
            call dstkio_r%read(stkname, ind_in_stk, build%imgbatch(ii))
        end do
        call dstkio_r%kill
    end subroutine discrete_read_imgbatch

    subroutine set_bp_range( params, build, cline )
        class(parameters), intent(inout) :: params
        class(builder),    intent(inout) :: build
        class(cmdline),    intent(in)    :: cline
        real, allocatable :: resarr(:), fsc_arr(:)
        real              :: fsc0143, fsc05
        real              :: mapres(params%nstates)
        integer           :: s, loc(1), lp_ind, arr_sz, fsc_sz, nyqcrop_ind
        type(string)      :: fsc_fname
        logical           :: fsc_bin_exists(params%nstates), all_fsc_bin_exist
        if( params%l_lpset )then
            ! set Fourier index range
            params%kfromto(2) = calc_fourier_index(params%lp, params%box, params%smpd)
            if( cline%defined('lpstop') )then
                params%kfromto(2) = min(params%kfromto(2),&
                    &calc_fourier_index(params%lpstop, params%box, params%smpd))
            endif
            ! FSC values are read anyway
            do s=1,params%nstates
                fsc_fname = FSC_FBODY//int2str_pad(s,2)//BIN_EXT
                if( file_exists(fsc_fname) )then
                    fsc_arr = file2rarr(fsc_fname)
                    fsc_sz  = size(build%fsc(s,:))
                    arr_sz  = size(fsc_arr)
                    if( fsc_sz == arr_sz )then
                        build%fsc(s,:) = fsc_arr(:)
                    else if( fsc_sz > arr_sz )then
                        ! padding
                        build%fsc(s,:arr_sz)   = fsc_arr(:)
                        build%fsc(s,arr_sz+1:) = 0.
                    else
                        ! clipping
                        build%fsc(s,:fsc_sz)   = fsc_arr(:fsc_sz)
                    endif
                    deallocate(fsc_arr)
                endif
            enddo
        else
            ! check all fsc_state*.bin exist
            all_fsc_bin_exist = .true.
            fsc_bin_exists    = .false.
            do s=1,params%nstates
                fsc_fname = FSC_FBODY//int2str_pad(s,2)//BIN_EXT
                fsc_bin_exists( s ) = file_exists(fsc_fname)
                if( build%spproj_field%get_pop(s, 'state') > 0 .and. .not.fsc_bin_exists(s))&
                    & all_fsc_bin_exist = .false.
            enddo
            if(build%spproj%is_virgin_field(params%oritype)) &
                all_fsc_bin_exist = (count(fsc_bin_exists)==params%nstates)
            ! set low-pass Fourier index limit
            if( all_fsc_bin_exist )then
                resarr = build%img%get_res()
                do s=1,params%nstates
                    if( fsc_bin_exists(s) )then
                        fsc_fname = FSC_FBODY//int2str_pad(s,2)//BIN_EXT
                        fsc_arr = file2rarr(fsc_fname)
                        fsc_sz  = size(build%fsc(s,:))
                        arr_sz  = size(fsc_arr)
                        if( fsc_sz == arr_sz )then
                            build%fsc(s,:) = fsc_arr(:)
                        else if( fsc_sz > arr_sz )then
                            ! padding
                            build%fsc(s,:arr_sz)   = fsc_arr(:)
                            build%fsc(s,arr_sz+1:) = 0.
                        else
                            ! clipping
                            build%fsc(s,:fsc_sz)   = fsc_arr(:fsc_sz)
                        endif
                        deallocate(fsc_arr)
                        call get_resolution(build%fsc(s,:), resarr, fsc05, fsc0143)
                        mapres(s) = fsc0143
                    else
                        ! empty state
                        mapres(s)           = 0.
                        build%fsc(s,:) = 0.
                    endif
                end do
                loc = minloc(mapres) ! best resolved
                if( params%nstates == 1 )then
                    lp_ind = get_find_at_crit(build%fsc(1,:), params%lplim_crit, incrreslim=params%l_incrreslim)
                else
                    lp_ind = get_find_at_crit(build%fsc(loc(1),:), params%lplim_crit)
                endif
                ! interpolation limit is NOT Nyqvist in correlation search
                params%kfromto(2) = calc_fourier_index(resarr(lp_ind), params%box, params%smpd)
            else if( build%spproj_field%isthere(params%fromp,'lp') )then
                params%kfromto(2) = calc_fourier_index(&
                    build%spproj_field%get(params%fromp,'lp'), params%box, params%smpd)
            else
                THROW_HARD('no method available for setting the low-pass limit. Need fsc file or lp find; set_bp_range')
            endif
            ! lpstop overrides any other method for setting the low-pass limit
            if( cline%defined('lpstop') )then
                params%kfromto(2) = min(params%kfromto(2), &
                    calc_fourier_index(params%lpstop, params%box, params%smpd))
            endif
            ! re-set the low-pass limit
            params%lp = calc_lowpass_lim(params%kfromto(2), params%box, params%smpd)
        endif
        ! making sure the resolution limit does not exceed limits of box_crop
        nyqcrop_ind       = calc_fourier_index(2.*params%smpd_crop, params%box_crop, params%smpd_crop)
        params%kfromto(2) = min(params%kfromto(2), nyqcrop_ind)
        params%lp         = calc_lowpass_lim(params%kfromto(2), params%box, params%smpd)
        ! update low-pass limit in project
        call build%spproj_field%set_all2single('lp',params%lp)
    end subroutine set_bp_range

    subroutine set_bp_range2D( params, build, cline, which_iter, frac_srch_space )
        class(parameters), intent(inout) :: params
        class(builder),    intent(inout) :: build
        class(cmdline),    intent(inout) :: cline
        integer,           intent(in)    :: which_iter
        real,              intent(in)    :: frac_srch_space
        real    :: lplim
        integer :: lpstart_find
        params%kfromto(1) = max(2,calc_fourier_index(params%hp, params%box, params%smpd))
        if( params%l_lpset )then
            lplim = params%lp
            params%kfromto(2) = calc_fourier_index(lplim, params%box_crop, params%smpd_crop)
        else
            if( trim(params%stream2d).eq.'yes' )then
                if( file_exists(params%frcs) )then
                    lplim = build%clsfrcs%estimate_lp_for_align()
                else
                    lplim = params%lplims2D(3)
                endif
                if( cline%defined('lpstop') ) lplim = max(lplim, params%lpstop)
            else
                if( file_exists(params%frcs) .and. which_iter >= LPLIM1ITERBOUND )then
                    lplim = build%clsfrcs%estimate_lp_for_align()
                else
                    if( which_iter < LPLIM1ITERBOUND )then
                        lplim = params%lplims2D(1)
                    else if( frac_srch_space >= FRAC_SH_LIM .and. which_iter > LPLIM3ITERBOUND )then
                        lplim = params%lplims2D(3)
                    else
                        lplim = params%lplims2D(2)
                    endif
                endif
            endif
            params%kfromto(2) = calc_fourier_index(lplim, params%box_crop, params%smpd_crop)
            ! to avoid pathological cases, fall-back on lpstart
            lpstart_find = calc_fourier_index(params%lpstart, params%box_crop, params%smpd_crop)
            if( lpstart_find > params%kfromto(2) ) params%kfromto(2) = lpstart_find
            lplim     = calc_lowpass_lim(params%kfromto(2), params%box_crop, params%smpd_crop)
            params%lp = lplim
        endif
        ! update low-pas limit in project
        call build%spproj_field%set_all2single('lp',lplim)
    end subroutine set_bp_range2D

    subroutine sample_ptcls4update( params, build, pfromto, l_incr_sampl, nptcls2update, pinds )
        class(parameters),    intent(in)    :: params
        class(builder),       intent(inout) :: build
        integer,              intent(in)    :: pfromto(2)
        logical,              intent(in)    :: l_incr_sampl
        integer,              intent(inout) :: nptcls2update
        integer, allocatable, intent(inout) :: pinds(:)
        type(class_sample),   allocatable   :: clssmp(:)
        type(string) :: fname
        fname = CLASS_SAMPLING_FILE
        if( params%l_update_frac )then
            if( trim(params%balance).eq.'yes' )then
                if( file_exists(fname) )then
                    call read_class_samples(clssmp, fname)
                else
                    THROW_HARD('File for class-biased sampling in fractional update: '//CLASS_SAMPLING_FILE//' does not exist!')
                endif
                ! balanced class sampling
                if( params%l_frac_best )then
                    call build%spproj_field%sample4update_class(clssmp, pfromto, params%update_frac,&
                    nptcls2update, pinds, l_incr_sampl, params%l_greedy_smpl, frac_best=params%frac_best)
                else
                    call build%spproj_field%sample4update_class(clssmp, pfromto, params%update_frac,&
                    nptcls2update, pinds, l_incr_sampl, params%l_greedy_smpl)
                endif
                call deallocate_class_samples(clssmp)
            else
                call build%spproj_field%sample4update_cnt(pfromto, params%update_frac,&
                    nptcls2update, pinds, l_incr_sampl)
            endif
        else
            ! we sample all state > 0
            call build%spproj_field%sample4update_all(pfromto, nptcls2update, pinds, l_incr_sampl)
        endif
        call fname%kill
    end subroutine sample_ptcls4update

    subroutine sample_ptcls4fillin( params, build, pfromto, l_incr_sampl, nptcls2update, pinds )
        class(parameters),    intent(in)    :: params
        class(builder),       intent(inout) :: build
        integer,              intent(in)    :: pfromto(2)
        logical,              intent(in)    :: l_incr_sampl
        integer,              intent(inout) :: nptcls2update
        integer, allocatable, intent(inout) :: pinds(:)
        call build%spproj_field%sample4update_fillin(pfromto, params%update_frac, nptcls2update, pinds, l_incr_sampl)
    end subroutine sample_ptcls4fillin

    !>  \ brief  prepares one particle image for alignment
    !          serial routine
    subroutine prepimg4align( params, build, iptcl, img, img_out, img_out_pd )
        class(parameters), intent(in)    :: params
        class(builder),    intent(inout) :: build
        integer,           intent(in)    :: iptcl
        class(image),      intent(inout) :: img
        class(image),      intent(inout) :: img_out
        class(image),      intent(inout) :: img_out_pd
        type(ctf)       :: tfun
        type(ctfparams) :: ctfparms
        real            :: shvec(2), crop_factor
        crop_factor = real(params%box_crop) / real(params%box)
        shvec(1)    = -build%spproj_field%get(iptcl, 'x') * crop_factor
        shvec(2)    = -build%spproj_field%get(iptcl, 'y') * crop_factor
        ! Phase-flipping
        ctfparms = build%spproj%get_ctfparams(params%oritype, iptcl)
        select case(ctfparms%ctfflag)
            case(CTFFLAG_NO, CTFFLAG_FLIP)
                ! fused noise normalization, FFT, clip & shift
                call img%norm_noise_fft_clip_shift(build%lmsk, img_out, shvec)
            case(CTFFLAG_YES)
                ctfparms%smpd = ctfparms%smpd / crop_factor != smpd_crop
                tfun          = ctf(ctfparms%smpd, ctfparms%kv, ctfparms%cs, ctfparms%fraca)
                ! fused noise normalization, FFT, clip, shift & CTF flip
                call img%norm_noise_fft_clip_shift_ctf_flip(build%lmsk, img_out, shvec, tfun, ctfparms)
            case DEFAULT
                THROW_HARD('unsupported CTF flag: '//int2str(ctfparms%ctfflag)//' prepimg4align')
        end select
        ! fused IFFT, mask, FFT
        call img_out%ifft_mask_pad_fft(params%msk_crop, img_out_pd)
    end subroutine prepimg4align

    subroutine prepimg4align_bench( params, build, iptcl, img, img_out, img_out_pd, rt_prep1, rt_prep2, rt_prep )
        class(parameters),    intent(in)     :: params
        class(builder),       intent(inout) :: build
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
        crop_factor = real(params%box_crop) / real(params%box)
        shvec(1)    = -build%spproj_field%get(iptcl, 'x') * crop_factor
        shvec(2)    = -build%spproj_field%get(iptcl, 'y') * crop_factor
        ! Phase-flipping
        ctfparms = build%spproj%get_ctfparams(params%oritype, iptcl)
        select case(ctfparms%ctfflag)
            case(CTFFLAG_NO, CTFFLAG_FLIP)
                ! fused noise normalization, FFT, clip & shift
                call img%norm_noise_fft_clip_shift(build%lmsk, img_out, shvec)
            case(CTFFLAG_YES)
                ctfparms%smpd = ctfparms%smpd / crop_factor != smpd_crop
                tfun          = ctf(ctfparms%smpd, ctfparms%kv, ctfparms%cs, ctfparms%fraca)
                ! fused noise normalization, FFT, clip, shift & CTF flip
                call img%norm_noise_fft_clip_shift_ctf_flip(build%lmsk, img_out, shvec, tfun, ctfparms)
            case DEFAULT
                THROW_HARD('unsupported CTF flag: '//int2str(ctfparms%ctfflag)//' prepimg4align')
        end select
        rt_prep1 = rt_prep1 + toc(t_prep1)
        t_prep2 = tic()
        ! fused IFFT, mask, FFT
        call img_out%ifft_mask_pad_fft(params%msk_crop, img_out_pd)
        rt_prep2 = rt_prep2 + toc(t_prep2)
        rt_prep   = rt_prep   + toc(t_prep)
    end subroutine prepimg4align_bench

    !>  \brief  Calculates the centering offset of the input cavg
    !>          cavg & particles centering is not performed
    subroutine calc_2Dref_offset( params, build, img, icls, centype, xyz )
        class(parameters), intent(in)    :: params
        class(builder),    intent(inout) :: build
        class(image),      intent(inout) :: img
        integer,           intent(in)    :: icls
        integer,           intent(in)    :: centype
        real,              intent(out)   :: xyz(3)
        real :: xy_cavg(2), crop_factor
        crop_factor = real(params%box_crop) / real(params%box)
        select case(centype)
            case(PARAMS_CEN)
                call build%spproj_field%calc_avg_offset2D(icls, xy_cavg)
                if( arg(xy_cavg) < CENTHRESH )then
                    xyz = 0.
                else if( arg(xy_cavg) > MAXCENTHRESH2D )then
                    xyz(1:2) = xy_cavg * crop_factor
                    xyz(3)   = 0.
                else
                    xyz = img%calc_shiftcen_serial(params%cenlp, params%msk_crop)
                    if( arg(xyz(1:2)/crop_factor - xy_cavg) > MAXCENTHRESH2D ) xyz = 0.
                endif
            case(SEG_CEN)
                call calc_cavg_offset(img, params%cenlp, params%msk_crop, xy_cavg)
                xyz = [xy_cavg(1), xy_cavg(2), 0.]
            case(MASS_CEN)
                xyz = img%calc_shiftcen_serial(params%cenlp, params%msk_crop)
        end select
        if( arg(xyz) < CENTHRESH ) xyz = 0.0
    end subroutine calc_2Dref_offset

    !>  \brief  Prepares one cluster centre image for alignment
    subroutine prep2Dref( params, build, img, icls, xyz, img_pd )
        class(parameters), intent(in)    :: params
        class(builder),    intent(inout) :: build
        class(image),      intent(inout) :: img, img_pd
        integer,           intent(in)    :: icls
        real,              intent(in)    :: xyz(3)
        integer :: filtsz
        real    :: frc(img%get_filtsz()), filter(img%get_filtsz()), crop_factor
        ! Cavg & particles centering
        if( arg(xyz) > CENTHRESH )then
            call img%fft()
            call img%shift2Dserial(xyz(1:2))
            crop_factor = real(params%box_crop) / real(params%box)
            call build%spproj_field%add_shift2class(icls, -xyz(1:2) / crop_factor)
        endif
        ! Filtering
        if( params%l_ml_reg )then
            ! no filtering
        else
            if( params%l_lpset.and.params%l_gauref )then
                ! Gaussian filter only applied when lp is set, FRC filtering turned off
                call img%fft
                call img%bpgau2d(0., params%gaufreq)
            else
                ! FRC-based filtering
                call build%clsfrcs%frc_getter(icls, frc)
                if( any(frc > 0.143) )then
                    filtsz = img%get_filtsz()
                    call fsc2optlp_sub(filtsz, frc, filter, merged=params%l_lpset)
                    call img%fft    ! in case the shift was not applied above
                    call img%apply_filter_serial(filter)
                endif
            endif
        endif
        ! ensure we are in real-space
        call img%ifft()
        ! mask, pad, FFT
        call img%mask_pad_fft(params%msk_crop, img_pd)
    end subroutine prep2Dref

    !>  \brief  initializes all volumes for reconstruction
    subroutine preprecvols( params, build )
        class(parameters), intent(in)    :: params
        class(builder),    intent(inout) :: build
        integer, allocatable :: pops(:)
        integer :: istate
        call build%spproj_field%get_pops(pops, 'state', maxn=params%nstates)
        do istate = 1, params%nstates
            if( pops(istate) > 0)then
                call build%eorecvols(istate)%new(params, build%spproj)
                call build%eorecvols(istate)%reset_all
            endif
        end do
        deallocate(pops)
    end subroutine preprecvols

    !>  \brief  destructs all volumes for reconstruction
    subroutine killrecvols( params, build )
        class(parameters), intent(in) :: params
        class(builder),    intent(inout) :: build
        integer :: istate
        do istate = 1, params%nstates
            call build%eorecvols(istate)%kill
        end do
    end subroutine killrecvols

    !>  \brief  determines the reference volume shift and map shifts back to particles
    !>          reference volume shifting is performed in shift_and_mask_refvol
    subroutine calcrefvolshift_and_mapshifts2ptcls( params, build, s, volfname, do_center, xyz, map_shift )
        class(parameters), intent(in) :: params
        class(builder),    intent(inout) :: build
        integer,           intent(in)    :: s
        class(string),     intent(in)    :: volfname
        logical,           intent(out)   :: do_center
        real,              intent(out)   :: xyz(3)
        logical,           intent(in)    :: map_shift
        real    :: crop_factor
        logical :: has_been_searched
        do_center = .true.
        ! centering
        if( params%center .eq. 'no' .or. params%nstates > 1 .or. &
            .not. params%l_doshift .or. params%pgrp(:1) .ne. 'c' .or. &
            params%l_filemsk .or. params%l_update_frac )then
            do_center = .false.
            xyz       = 0.
            return
        endif
        ! taking care of volume dimensions
        call build%vol%read_and_crop(volfname, params%smpd, params%box_crop, params%smpd_crop)
        ! offset
        xyz = build%vol%calc_shiftcen(params%cenlp,params%msk_crop)
        if( params%pgrp .ne. 'c1' ) xyz(1:2) = 0.     ! shifts only along z-axis for C2 and above
        if( arg(xyz) <= CENTHRESH )then
            do_center = .false.
            xyz       = 0.
            return
        endif
        if( map_shift )then
            ! map back to particle oritentations
            has_been_searched = .not.build%spproj%is_virgin_field(params%oritype)
            if( has_been_searched )then
                crop_factor = real(params%box) / real(params%box_crop)
                call build%spproj_field%map3dshift22d(-xyz(:)*crop_factor, state=s)
            endif
        endif
    end subroutine calcrefvolshift_and_mapshifts2ptcls

    subroutine report_resolution( params, build, state )
        class(parameters), intent(in)  :: params
        class(builder),    intent(in)  :: build
        integer,           intent(out) :: state
        real :: res_fsc05, res_fsc0143
        real, allocatable  :: res(:)
        integer :: s
        real    :: lpest(params%nstates)
        res = get_resarr(params%box_crop, params%smpd_crop)
        do s = 1, params%nstates
            call get_resolution_at_fsc(build%fsc(s,:), res, 0.5, lpest(s))
        end do
        state = minloc(lpest,dim=1)
        call get_resolution_at_fsc(build%fsc(state,:), res, 0.50,  res_fsc05)
        call get_resolution_at_fsc(build%fsc(state,:), res, 0.143, res_fsc0143)
        call build%spproj_field%set_all2single('res', res_fsc0143)
    end subroutine report_resolution

    subroutine estimate_lp_refvols( params, build, cline, lpstart, lpstop, state )
        use simple_polarft_calc, only: polarft_estimate_lplim3D
        class(parameters), intent(inout) :: params
        class(builder),    intent(inout) :: build
        class(cmdline),    intent(in)    :: cline
        real,              intent(in)    :: lpstart, lpstop
        integer,           intent(in)    :: state
        type(string)       :: vol_even, vol_odd
        logical, parameter :: DEBUG=.true.
        type(image) :: mskvol
        integer     :: npix
        real        :: lpopt
        if( lpstart < lpstop )then
            THROW_HARD('Invalid low-pass range ordering: lpstart must be >= lpstop')
        endif
        if( (trim(params%polar).eq.'yes').and.(.not.cline%defined('vol1')) )then
            ! POLAR REPRESENTATION
            call polarft_estimate_lplim3D(params%box_crop, params%smpd_crop, lpstart, lpstop, lpopt)
        else
            ! CARTSESIAN REPRESENTATION
            if( params%l_filemsk )then
                ! read 3D envelope mask
                call mskvol%read_and_crop(params%mskfile, params%smpd, params%box_crop, params%smpd_crop)
                call mskvol%remove_edge
            else
                ! finding optimal lp over all states
                call mskvol%disc([params%box_crop,  params%box_crop, params%box_crop], params%smpd_crop,&
                    &params%msk_crop, npix )
            endif
            ! read volumes
            vol_even = params%vols_even(state)
            vol_odd  = params%vols_odd(state)
            call build%vol%read_and_crop(    vol_even, params%smpd, params%box_crop, params%smpd_crop)
            call build%vol_odd%read_and_crop(vol_odd,  params%smpd, params%box_crop, params%smpd_crop)
            ! estimate low-pass limit
            call estimate_lplim3D(build%vol_odd, build%vol, mskvol, lpstart, lpstop, lpopt)
        endif        
        ! update low-pass limit in project
        call build%spproj_field%set_all2single('lp_est', lpopt)
        if( params%l_lpauto )then
            ! re-set the low-pass limit
            params%lp = lpopt
            ! update the Fourier index limit
            params%kfromto(2) = calc_fourier_index(params%lp, params%box_crop, params%smpd_crop)
            ! update low-pass limit in project
            call build%spproj_field%set_all2single('lp',params%lp)
        endif
        ! destruct
        call mskvol%kill
    end subroutine estimate_lp_refvols

    subroutine read_mask_and_filter_refvols( params, build, s )
        class(parameters), intent(in)    :: params
        class(builder),    intent(inout) :: build
        integer,           intent(in)    :: s
        type(string)         :: vol_even, vol_odd, vol_avg
        logical, allocatable :: l_msk(:,:,:)
        real    :: cur_fil(params%box_crop)
        integer :: filtsz
        ! READ
        vol_even = params%vols_even(s)
        vol_odd  = params%vols_odd(s)
        vol_avg  = params%vols(s)
        call build%vol%read_and_crop(   vol_even, params%smpd, params%box_crop, params%smpd_crop)
        call build%vol_odd%read_and_crop(vol_odd, params%smpd, params%box_crop, params%smpd_crop)
        if( s == 1 .and. params%l_filemsk )then
            ! read 3D envelope mask
            call build%mskvol%read_and_crop(params%mskfile, params%smpd, params%box_crop, params%smpd_crop)
        endif
        ! noise regularization
        if( params%l_noise_reg )then
            call build%vol%add_gauran(params%eps)
            call build%vol_odd%add_gauran(params%eps)
        endif
        ! MASK
        if( params%l_filemsk )then
            ! envelope masking
            call build%vol%zero_env_background(build%mskvol)
            call build%vol_odd%zero_env_background(build%mskvol)
            call build%vol%mul(build%mskvol)
            call build%vol_odd%mul(build%mskvol)
        else
            ! circular masking
            call build%vol%mask3D_soft(params%msk_crop, backgr=0.0)
            call build%vol_odd%mask3D_soft(params%msk_crop, backgr=0.0)
        endif
        ! FILTER
        if( params%l_icm )then
            if( params%l_filemsk )then
                l_msk = build%mskvol%bin2logical()
                call build%vol%ICM3D_eo(build%vol_odd, params%lambda, l_msk)
            else
                call build%vol%ICM3D_eo(build%vol_odd, params%lambda)
            endif
            if( params%l_lpset .or. trim(params%combine_eo).eq.'yes' )then ! no independent volume registration, so average eo pairs
                call build%vol%add(build%vol_odd)
                call build%vol%mul(0.5)
                call build%vol_odd%copy_fast(build%vol)
            endif
            ! FT
            call build%vol%fft
            call build%vol_odd%fft
        else if( params%l_lpset )then
            ! read average volume that will occupy both even and odd
            call build%vol%read_and_crop(vol_avg, params%smpd, params%box_crop, params%smpd_crop)
            ! noise regularization
            if( params%l_noise_reg )then
                call build%vol%add_gauran(params%eps)
            endif
            ! mask again, BP filter performed below
            if( params%l_filemsk )then
                ! envelope masking
                call build%vol%zero_env_background(build%mskvol)
                call build%vol%mul(build%mskvol)
            else
                ! circular masking
                call build%vol%mask3D_soft(params%msk_crop, backgr=0.0)
            endif
            ! FT & odd <- even
            call build%vol%fft
            call build%vol_odd%copy_fast(build%vol)
        else
            ! FT
            call build%vol%fft
            call build%vol_odd%fft
        endif
        if( params%l_ml_reg )then
            ! filtering done when volumes are assembled
        else if( params%l_icm )then
            ! filtering done above
        else if( params%l_lpset )then
            if( trim(params%gauref)=='yes' )then
                ! Gaussian filter
                call build%vol%bpgau3D(0., params%gaufreq)
                call build%vol_odd%bpgau3D(0., params%gaufreq)
            else
                ! Cosine low-pass filter, works best for nanoparticles
                call build%vol%bp(0., params%lp)
                call build%vol_odd%bp(0., params%lp)
            endif
        else
            ! Optimal filter
            filtsz = build%vol%get_filtsz()
            if( any(build%fsc(s,:) > 0.143) )then
                call fsc2optlp_sub(filtsz,build%fsc(s,:),cur_fil)
                call build%vol%apply_filter(cur_fil)
                call build%vol_odd%apply_filter(cur_fil)
            endif
        endif
    end subroutine read_mask_and_filter_refvols

    !>  \brief  grids one particle image to the volume
    subroutine grid_ptcl( build, fpl, se, o )
        class(builder),     intent(inout) :: build
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
        if( pw > TINY ) call build%eorecvols(s)%grid_plane(se, o, fpl, eo, pw)
    end subroutine grid_ptcl

    !> volumetric 3d reconstruction
    subroutine calc_3Drec( params, build, cline, nptcls, pinds )
        use simple_imgarr_utils, only: alloc_imgarr, dealloc_imgarr
        class(parameters), intent(inout) :: params
        class(builder),    intent(inout) :: build
        class(cmdline),    intent(inout) :: cline
        integer,           intent(in)    :: nptcls
        integer,           intent(in)    :: pinds(nptcls)
        type(fplane_type), allocatable   :: fpls(:)
        integer :: batchlims(2), ibatch, batchsz
        logical :: DEBUG = .false.
        integer(timer_int_kind) :: t, t0
        real(timer_int_kind)    :: t_init, t_read, t_prep, t_grid, t_tot
        if( DEBUG ) t0 = tic()
        ! Initialize objects for recontruction
        if( DEBUG ) t = tic()
        call init_rec(params, build, MAXIMGBATCHSZ, fpls)
        ! Prep batch image objects
        call prepimgbatch(params, build, MAXIMGBATCHSZ)
        if( DEBUG ) t_init = toc(t)
        ! gridding batch loop
        if( DEBUG ) then
            t_read = 0.d0
            t_prep = 0.d0
            t_grid = 0.d0
        endif
        do ibatch = 1,nptcls,MAXIMGBATCHSZ
            batchlims = [ibatch, min(nptcls, ibatch+MAXIMGBATCHSZ-1)]
            batchsz   = batchlims(2) - batchlims(1) + 1
            ! read images
            if( DEBUG ) t = tic()
            call discrete_read_imgbatch(params, build, nptcls, pinds, batchlims)
            if( DEBUG ) t_read = t_read + toc(t)
            ! preprocess images into padded objects
            if( DEBUG ) t = tic()
            call prep_imgs4rec(params, build, batchsz, build%imgbatch(:batchsz),&
                                &pinds(batchlims(1):batchlims(2)), fpls(:batchsz))
            if( DEBUG ) t_prep = t_prep + toc(t)
            ! insert padded slices into lattice
            if( DEBUG ) t = tic()
            call update_rec(params, build, batchsz, pinds(batchlims(1):batchlims(2)), fpls(:batchsz))
            if( DEBUG ) t_grid = t_grid + toc(t)
        end do
        ! Normalize volumes and deallocate temporary objects
        call finalize_rec(params, build, cline, fpls)
        if( DEBUG .and. (params%part==1) )then
            t_tot = toc(t0)
            print *,'Init          : ', t_init
            print *,'Read          : ', t_read
            print *,'Prep          : ', t_prep
            print *,'Grid          : ', t_grid
            print *,'Total rec time: ', t_tot
        endif
    end subroutine calc_3Drec

    !>  Initiates objects required for online volumetric 3d reconstruction
    !>  Does not read images
    subroutine init_rec( params, build, maxbatchsz, fplanes )
        use simple_imgarr_utils, only: alloc_imgarr
        class(parameters),              intent(in)    :: params
        class(builder),                 intent(inout) :: build
        integer,                        intent(in)    :: maxbatchsz
        type(fplane_type), allocatable, intent(inout) :: fplanes(:)
        ! sanity check for ml_reg
        if( params%l_ml_reg )then
            if( .not. allocated(build%esig%sigma2_noise) )then
                THROW_HARD('build%esig%sigma2_noise is not allocated while ml_reg is enabled; calc_3Drec')
            endif
        endif
        ! init volumes
        call preprecvols(params, build)
        ! allocate convenience CTF & memory aligned objects
        if( allocated(fplanes) )  deallocate(fplanes)
        allocate(fplanes(maxbatchsz))
        ! heap of padded images
        call alloc_imgarr(nthr_glob, [params%boxpd, params%boxpd, 1], params%smpd, build%img_pad_heap)
    end subroutine init_rec

    !> Preprocess particle images for online volumetric 3d reconstruction
    subroutine prep_imgs4rec( params, build, nptcls, ptcl_imgs, pinds, fplanes )
        class(parameters), intent(in)    :: params
        class(builder),    intent(inout) :: build
        integer,           intent(in)    :: nptcls
        class(image),      intent(inout) :: ptcl_imgs(nptcls)
        integer,           intent(in)    :: pinds(nptcls)
        type(fplane_type), intent(inout) :: fplanes(nptcls)
        type(ctfparams) :: ctfparms(nthr_glob)
        real      :: shift(2)
        integer   :: iptcl, i, ithr, kfromto(2)
        ! logical/physical adress mapping because overwritten by polarize
        call memoize_ft_maps([params%boxpd, params%boxpd, 1], params%smpd)
        ! gridding batch loop
        kfromto = build%esig%get_kfromto()
        !$omp parallel do default(shared) private(i,ithr,iptcl,shift) schedule(static) proc_bind(close)
        do i = 1,nptcls
            ithr   = omp_get_thread_num() + 1
            iptcl  = pinds(i)
            call ptcl_imgs(i)%norm_noise_taper_edge_pad_fft(build%lmsk, build%img_pad_heap(ithr))
            ctfparms(ithr) = build%spproj%get_ctfparams(params%oritype, iptcl)
            shift = build%spproj_field%get_2Dshift(iptcl)
            if( params%l_ml_reg )then
                call build%img_pad_heap(ithr)%gen_fplane4rec(kfromto, params%smpd_crop, ctfparms(ithr),&
                &shift, iptcl, fplanes(i), build%esig%sigma2_noise(kfromto(1):kfromto(2),iptcl))
            else
                call build%img_pad_heap(ithr)%gen_fplane4rec(kfromto, params%smpd_crop, ctfparms(ithr),&
                &shift, iptcl, fplanes(i))
            endif
        end do
        !$omp end parallel do
    end subroutine prep_imgs4rec

    !> volumetric 3d reconstruction
    subroutine update_rec( params, build, nptcls, pinds, fplanes )
        class(parameters), intent(in)    :: params
        class(builder),    intent(inout) :: build
        integer,           intent(in)    :: nptcls
        integer,           intent(in)    :: pinds(nptcls)
        type(fplane_type), intent(inout) :: fplanes(nptcls)
        type(ori) :: orientation
        integer   :: iptcl, i
        call memoize_ft_maps([params%boxpd, params%boxpd, 1], params%smpd)
        ! gridding
        do i = 1,nptcls
            iptcl  = pinds(i)
            call build%spproj_field%get_ori(iptcl, orientation)
            if( orientation%isstatezero() ) cycle
            call grid_ptcl(build, fplanes(i), build%pgrpsyms, orientation)
        end do
        call orientation%kill
    end subroutine update_rec

    !> volumetric 3d reconstruction
    subroutine finalize_rec( params, build, cline, fplanes )
        use simple_imgarr_utils, only: alloc_imgarr, dealloc_imgarr
        class(parameters),              intent(inout) :: params
        class(builder),                 intent(inout) :: build
        class(cmdline),                 intent(inout) :: cline
        type(fplane_type), allocatable, intent(inout) :: fplanes(:)
        integer :: i
        ! deallocate convenience objects
        do i = 1,size(fplanes)
            if( allocated(fplanes(i)%cmplx_plane) ) deallocate(fplanes(i)%cmplx_plane)
            if( allocated(fplanes(i)%ctfsq_plane) ) deallocate(fplanes(i)%ctfsq_plane)
        end do
        deallocate(fplanes)
        call dealloc_imgarr(build%img_pad_heap)
        ! normalise structure factors
        call memoize_ft_maps([params%boxpd, params%boxpd, 1], params%smpd)
        call norm_struct_facts(params, build, cline)
        ! destruct
        call forget_ft_maps
        call killrecvols(params, build)
    end subroutine finalize_rec

    subroutine norm_struct_facts( params, build, cline )
        use simple_gridding, only: prep3D_inv_instrfun4mul
        class(parameters), intent(inout) :: params
        class(builder),    intent(inout) :: build
        class(cmdline),    intent(inout) :: cline
        type(string) :: recname, volname, volname_prev, volname_prev_even
        type(string) :: volname_prev_odd, str_state, str_iter, fsc_txt_file, eonames(2)
        real, allocatable :: fsc(:)
        type(image) :: vol_prev_even, vol_prev_odd, gridcorr_img
        integer     :: s, find4eoavg, ldim_pd(3), ldim(3)
        real        :: res05s(params%nstates), res0143s(params%nstates)
        real        :: weight_prev, update_frac_trail_rec
        ! init
        ldim    = [params%box_crop,  params%box_crop,  params%box_crop]
        ldim_pd = [params%box_croppd,params%box_croppd,params%box_croppd]
        call build%vol%new(ldim,params%smpd_crop)
        call build%vol2%new(ldim,params%smpd_crop)
        res0143s = 0.
        res05s   = 0.
        ! read in previous reconstruction when trail_rec==yes
        update_frac_trail_rec = 1.0
        if( .not. params%l_distr_worker .and. params%l_trail_rec )then
            if( cline%defined('ufrac_trec') )then
                update_frac_trail_rec = params%ufrac_trec
            else
                update_frac_trail_rec = build%spproj_field%get_update_frac()
            endif
        endif
        ! Prep for correction of the shape of the interpolator
        gridcorr_img = prep3D_inv_instrfun4mul(ldim, ldim_pd, params%smpd_crop)
        ! cycle through states
        do s=1,params%nstates
            if( build%spproj_field%get_pop(s, 'state') == 0 )then
                ! empty state
                build%fsc(s,:) = 0.
                cycle
            endif
            call build%eorecvols(s)%compress_exp
            if( params%l_distr_worker )then
                call build%eorecvols(s)%write_eos(string(VOL_FBODY)//int2str_pad(s,2)//'_part'//int2str_pad(params%part,params%numlen))
            else
                ! global volume name update
                recname = VOL_FBODY//int2str_pad(s,2)
                volname = recname//MRC_EXT
                eonames(1) = recname//'_even'//MRC_EXT
                eonames(2) = recname//'_odd'//MRC_EXT
                if( params%l_ml_reg )then
                    ! the sum is done after regularization
                else
                    call build%eorecvols(s)%sum_eos
                endif
                if( params%l_trail_rec )then
                    if( .not. cline%defined('vol'//int2str(s)) ) THROW_HARD('vol'//int2str(s)//'required in norm_struct_facts cline when trail_rec==yes')
                    volname_prev      = cline%get_carg('vol'//int2str(s))
                    volname_prev_even = add2fbody(volname_prev, MRC_EXT, '_even')
                    volname_prev_odd  = add2fbody(volname_prev, MRC_EXT, '_odd')
                    if( .not. file_exists(volname_prev_even) ) THROW_HARD('File: '//volname_prev_even%to_char()//' does not exist!')
                    if( .not. file_exists(volname_prev_odd)  ) THROW_HARD('File: '//volname_prev_odd%to_char()//' does not exist!')
                    call vol_prev_even%read_and_crop(volname_prev_even, params%smpd, params%box_crop, params%smpd_crop)
                    call vol_prev_odd %read_and_crop(volname_prev_odd,  params%smpd, params%box_crop, params%smpd_crop)
                    if( allocated(fsc) ) deallocate(fsc)
                    call build%eorecvols(s)%calc_fsc4sampl_dens_correct(vol_prev_even, vol_prev_odd, fsc)
                    call build%eorecvols(s)%sampl_dens_correct_eos(s, eonames(1), eonames(2), find4eoavg, fsc)
                else 
                    call build%eorecvols(s)%sampl_dens_correct_eos(s, eonames(1), eonames(2), find4eoavg)
                endif
                str_state = int2str_pad(s,2)
                if( cline%defined('which_iter') )then
                    str_iter     = int2str_pad(params%which_iter,3)
                    fsc_txt_file = 'RESOLUTION_STATE'//str_state%to_char()//'_ITER'//str_iter%to_char()
                else
                    fsc_txt_file = 'RESOLUTION_STATE'//str_state%to_char()
                endif
                call build%eorecvols(s)%write_fsc2txt(fsc_txt_file)
                if( params%l_ml_reg )then
                    call build%eorecvols(s)%sum_eos
                endif
                call build%eorecvols(s)%get_res(res05s(s), res0143s(s))
                call build%eorecvols(s)%sampl_dens_correct_sum(build%vol)
                ! need to put the sum back at lowres for the eo pairs
                call build%vol%fft
                call build%vol2%zero_and_unflag_ft
                call build%vol2%read(eonames(1))
                call build%vol2%fft()
                call build%vol2%insert_lowres(build%vol, find4eoavg)
                call build%vol2%ifft()
                call build%vol2%mul(gridcorr_img)
                call build%vol2%write(eonames(1), del_if_exists=.true.)
                call build%vol2%zero_and_unflag_ft
                call build%vol2%read(eonames(2))
                call build%vol2%fft()
                call build%vol2%insert_lowres(build%vol, find4eoavg)
                call build%vol2%ifft()
                call build%vol2%mul(gridcorr_img)
                call build%vol2%write(eonames(2), del_if_exists=.true.)
                ! merged volume
                call build%vol%ifft
                call build%vol%mul(gridcorr_img)
                call build%vol%write(volname, del_if_exists=.true.)
                if( params%l_trail_rec .and. update_frac_trail_rec < 0.99 )then
                    call build%vol%read(eonames(1))  ! even current
                    call build%vol2%read(eonames(2)) ! odd current
                    weight_prev = 1. - update_frac_trail_rec
                    call vol_prev_even%mul(weight_prev)
                    call vol_prev_odd%mul (weight_prev)
                    call build%vol%mul(update_frac_trail_rec)
                    call build%vol2%mul(update_frac_trail_rec)
                    call build%vol%add(vol_prev_even)
                    call build%vol2%add(vol_prev_odd)
                    call build%vol%write(eonames(1))  ! even trailed
                    call build%vol2%write(eonames(2)) ! odd trailed
                    call vol_prev_even%kill
                    call vol_prev_odd%kill
                endif
                call build%vol%fft()
                call build%vol2%fft()
                ! updating command-line and parameters objects accordingly (needed in multi-stage wflows)
                params%vols(s) = volname
                call cline%set('vol'//int2str(s), params%vols(s))
            endif
        end do
        call build%vol2%kill
        call gridcorr_img%kill
    end subroutine norm_struct_facts

    subroutine prepare_refs_sigmas_ptcls( params, build, cline, ptcl_imgs, ptcl_imgs_pad, batchsz, which_iter, do_polar )
        class(parameters),        intent(inout) :: params
        class(builder),           intent(inout) :: build
        class(cmdline),           intent(in)    :: cline !< command line
        type(image), allocatable, intent(inout) :: ptcl_imgs(:)
        type(image), allocatable, intent(inout) :: ptcl_imgs_pad(:)
        integer,                  intent(in)    :: batchsz
        integer,                  intent(in)    :: which_iter
        logical,     optional,    intent(in)    :: do_polar
        real, allocatable :: gaufilter(:)
        type(string)      :: fname
        integer           :: ithr, iproj, nrefs, filtsz, state
        logical           :: l_polar, l_filtrefs
        l_polar = .false.
        if( present(do_polar) ) l_polar = do_polar
        ! PREPARATION OF pftc AND REFERENCES
        if( l_polar )then
            ! Resolution limit estimation
            call report_resolution(params, build, state)
            if( cline%defined('lpstart') .and. cline%defined('lpstop') )then
                call estimate_lp_refvols(params, build, cline, params%lpstart, params%lpstop, state)
            endif
            ! Calculator init
            nrefs     = params%nspace * params%nstates
            call build%pftc%new(params, nrefs, [1,batchsz], params%kfromto)
            ! Read polar references
            call build%pftc%polar_cavger_new(.true.)
            call build%pftc%polar_cavger_read_all(string(POLAR_REFS_FBODY//BIN_EXT))
            call build%clsfrcs%read(string(FRCS_FILE))
            ! prepare filter
            l_filtrefs = .false.
            if(trim(params%gauref).eq.'yes')then
                l_filtrefs = .true.
                filtsz = build%clsfrcs%get_filtsz()
                allocate(gaufilter(filtsz),source=0.)
                call gaussian_filter(params%gaufreq, params%smpd, params%box, gaufilter)
            endif
            ! PREPARATION OF REFERENCES IN pftc
            !$omp parallel do default(shared) private(iproj)&
            !$omp schedule(static) proc_bind(close)
            do iproj = 1,params%nspace
                if( l_filtrefs ) call build%pftc%polar_filterrefs(iproj, gaufilter)
                ! transfer to pftc
                if( params%l_lpset )then
                    ! put the merged class average in both even and odd positions
                    call build%pftc%polar_cavger_set_ref_pft(iproj, 'merged')
                    call build%pftc%cp_even2odd_ref(iproj)
                else
                    ! transfer e/o refs to pftc
                    call build%pftc%polar_cavger_set_ref_pft(iproj, 'even')
                    call build%pftc%polar_cavger_set_ref_pft(iproj, 'odd')
                endif
            end do
            !$omp end parallel do
            ! Center prior to memoization
            ! if( (trim(params%center)=='yes')   .and. (trim(params%center_type)=='params') .and.&
            !         &(params%pgrp(:2).eq.'c1') .and. (.not.params%l_update_frac)          .and.&
            !         &(params%nstates==1)       .and.       params%l_doshift )then
            !     call build%pftc%center_3Dpolar_refs(build%spproj_field, build%eulspace)
            ! endif
            ! Memoize
            call build%pftc%memoize_refs
        else
            ! (if needed) estimating lp (over all states) and reseting params%lp and params%kfromto
            call prepare_polar_references(params, build, cline, batchsz)
        endif
        ! PREPARATION OF SIGMAS
        if( params%cc_objfun == OBJFUN_EUCLID )then
            fname = SIGMA2_FBODY//int2str_pad(params%part,params%numlen)//'.dat'
            call build%esig%new(params, build%pftc, fname, params%box)
            call build%esig%read_part(  build%spproj_field)
            call build%esig%read_groups(build%spproj_field)
        end if
        ! PREPARATION OF PARTICLES
        call prepimgbatch(params, build, batchsz)
        allocate(ptcl_imgs(nthr_glob), ptcl_imgs_pad(nthr_glob))
        !$omp parallel do default(shared) private(ithr) schedule(static) proc_bind(close)
        do ithr = 1,nthr_glob
            call ptcl_imgs(ithr)%new(    [params%box_crop,  params%box_crop,  1], params%smpd_crop, wthreads=.false.)
            call ptcl_imgs_pad(ithr)%new([params%box_croppd,params%box_croppd,1], params%smpd_crop, wthreads=.false.)
        enddo
        !$omp end parallel do
        call build%vol%kill
        call build%vol_odd%kill
        call build%vol2%kill
        if( allocated(gaufilter) ) deallocate(gaufilter)
    end subroutine prepare_refs_sigmas_ptcls

    subroutine prepare_polar_references( params, build, cline, batchsz )
        class(parameters),   intent(inout) :: params
        class(builder),      intent(inout) :: build
        class(cmdline),      intent(in)    :: cline !< command line
        integer,             intent(in)    :: batchsz
        real      :: xyz(3)
        integer   :: s, nrefs, state
        logical   :: do_center, l_prob_align_mode
        call report_resolution(params, build, state)
        if( cline%defined('lpstart') .and. cline%defined('lpstop') )then
            call estimate_lp_refvols(params, build, cline, params%lpstart, params%lpstop, state)
        endif
        ! pftc
        nrefs = params%nspace * params%nstates
        call build%pftc%new(params, nrefs, [1,batchsz], params%kfromto)
        select case(trim(params%refine))
            case('prob','prob_state','prob_neigh')
                l_prob_align_mode = .true.
            case DEFAULT
                l_prob_align_mode = .false.
        end select
        ! read reference volumes and create polar projections
        do s=1,params%nstates
            if( l_prob_align_mode )then
                ! already mapping shifts in prob_tab with shared-memory execution
                call calcrefvolshift_and_mapshifts2ptcls(params, build, s, params%vols(s), do_center, xyz, map_shift=l_distr_worker_glob)
            else
                call calcrefvolshift_and_mapshifts2ptcls(params, build, s, params%vols(s), do_center, xyz, map_shift=.true.)
            endif
            call read_mask_and_filter_refvols(params, build, s)
            ! PREPARE E/O VOLUMES
            ! do even/odd in separate passes to reduce memory usage when padding is large (2X)
            ! PREPARE EVEN REFERENCES
            call build%vol_pad%new([params%box_croppd, params%box_croppd, params%box_croppd],&
            &params%smpd_crop, wthreads=.true.)
            if( do_center )then
                call build%vol%fft()
                call build%vol%shift([xyz(1),xyz(2),xyz(3)])
            endif
            ! back to real space
            call build%vol%ifft()
            ! FT volume
            call build%vol%pad_fft(build%vol_pad)
            ! expand for fast interpolation & correct for norm when clipped
            call build%vol_pad%expand_cmat(params%box)
            call vol_pad2ref_pfts(build%pftc, build%vol_pad, build%eulspace, s, .true., build%l_resmsk)
            call build%vol_pad%kill
            call build%vol_pad%kill_expanded
            ! PREPARE ODD REFERENCES
            call build%vol_odd_pad%new([params%box_croppd, params%box_croppd, params%box_croppd],&
            &params%smpd_crop, wthreads=.true.)
            if( do_center )then
                call build%vol_odd%fft()
                call build%vol_odd%shift([xyz(1),xyz(2),xyz(3)])
            endif
            ! back to real space
            call build%vol_odd%ifft()
            ! FT volume
            call build%vol_odd%pad_fft(build%vol_odd_pad)
            ! expand for fast interpolation & correct for norm when clipped
            call build%vol_odd_pad%expand_cmat(params%box)
            call vol_pad2ref_pfts(build%pftc, build%vol_odd_pad, build%eulspace, s, .false., build%l_resmsk)
            call build%vol_odd_pad%kill
            call build%vol_odd_pad%kill_expanded
        end do
        call build%pftc%memoize_refs
    end subroutine prepare_polar_references

    subroutine build_batch_particles( params, build, nptcls_here, pinds_here, tmp_imgs, tmp_imgs_pad, imgs4rec )
        class(parameters),      intent(in)    :: params
        class(builder),         intent(inout) :: build
        integer,                intent(in)    :: nptcls_here
        integer,                intent(in)    :: pinds_here(nptcls_here)
        class(image),           intent(inout) :: tmp_imgs(params%nthr), tmp_imgs_pad(params%nthr)
        class(image), optional, intent(inout) :: imgs4rec(nptcls_here)
        complex, allocatable :: pft(:,:)
        integer :: iptcl_batch, iptcl, ithr
        logical :: l_backup_imgs
        l_backup_imgs = present(imgs4rec)
        ! reassign particles indices & associated variables
        call build%pftc%reallocate_ptcls(nptcls_here, pinds_here)
        call discrete_read_imgbatch(params, build, nptcls_here, pinds_here, [1,nptcls_here])
        ! mask memoization for prepimg4align
        call tmp_imgs(1)%memoize_mask_coords
        ! memoize CTF stuff
        call memoize_ft_maps(tmp_imgs(1)%get_ldim(), tmp_imgs(1)%get_smpd())   
        ! memoize for polarize_oversamp
        call tmp_imgs_pad(1)%memoize4polarize_oversamp(build%pftc%get_pdim_interp())
        !$omp parallel do default(shared) private(iptcl,iptcl_batch,ithr,pft) schedule(static) proc_bind(close)
        do iptcl_batch = 1,nptcls_here
            ithr  = omp_get_thread_num() + 1
            iptcl = pinds_here(iptcl_batch)
            ! prep
            if( l_backup_imgs ) call imgs4rec(iptcl_batch)%copy_fast(build%imgbatch(iptcl_batch))
            call prepimg4align(params, build, iptcl, build%imgbatch(iptcl_batch), tmp_imgs(ithr), tmp_imgs_pad(ithr))
            ! transfer to polar coordinates
            pft = build%pftc%allocate_ptcl_pft()
            call tmp_imgs_pad(ithr)%polarize_oversamp(pft, mask=build%l_resmsk)
            call build%pftc%set_ptcl_pft(iptcl, pft)
            deallocate(pft)
            ! e/o flags
            call build%pftc%set_eo(iptcl, nint(build%spproj_field%get(iptcl,'eo'))<=0 )
        end do
        !$omp end parallel do
        call build%pftc%create_polar_absctfmats(build%spproj, 'ptcl3D')
        ! Memoize particles FFT parameters
        call build%pftc%memoize_ptcls
        ! destruct
        call forget_ft_maps
    end subroutine build_batch_particles

end module simple_strategy2D3D_common
