! common PRIME2D/PRIME3D routines used primarily by the Hadamard matchers
module simple_strategy2D3D_common
include 'simple_lib.f08'
use simple_image,             only: image
use simple_cmdline,           only: cmdline
use simple_builder,           only: build_glob
use simple_parameters,        only: params_glob
use simple_stack_io,          only: stack_io
use simple_discrete_stack_io, only: dstack_io
use simple_polarft_corrcalc,  only: pftcc_glob
implicit none

public :: prepimgbatch, killimgbatch, read_imgbatch, set_bp_range, set_bp_range2D, sample_ptcls4update, prepimg4align,&
&prep2Dref, preprecvols, killrecvols, calcrefvolshift_and_mapshifts2ptcls, read_and_filter_refvols,&
&preprefvol, grid_ptcl, calc_3Drec, calc_projdir3Drec, norm_struct_facts, discrete_read_imgbatch
private
#include "simple_local_flags.inc"

interface read_imgbatch
    module procedure read_imgbatch_1
    module procedure read_imgbatch_2
    module procedure read_imgbatch_3
end interface read_imgbatch

real, parameter :: SHTHRESH = 0.001
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
                    if( .not. stkio_r%stk_is_open() )then
                        call stkio_r%open(stkname, params_glob%smpd, 'read')
                    else if( .not. stkio_r%same_stk(stkname, [params_glob%box,params_glob%box,1]) )then
                        call stkio_r%close
                        call stkio_r%open(stkname, params_glob%smpd, 'read')
                    endif
                    call stkio_r%read(ind_in_stk, build_glob%imgbatch(ind_in_batch))
                endif
            end do
        else
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
        endif
        call stkio_r%close
    end subroutine read_imgbatch_1

    subroutine read_imgbatch_2( n, pinds, batchlims )
        integer,          intent(in)  :: n, pinds(n), batchlims(2)
        character(len=:), allocatable :: stkname
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
        integer,          intent(in)    :: iptcl
        type(image),      intent(inout) :: img
        character(len=:), allocatable   :: stkname
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

    subroutine discrete_read_imgbatch( n, pinds, batchlims, use_denoised )
        integer,          intent(in)  :: n, pinds(n), batchlims(2)
        type(dstack_io)               :: dstkio_r
        character(len=:), allocatable :: stkname, stkname_den
        logical,          optional    :: use_denoised
        integer :: ind_in_stk, i, ii
        logical :: uuse_denoised
        uuse_denoised = .false.
        if( present(use_denoised) ) uuse_denoised = use_denoised
        call dstkio_r%new(params_glob%smpd, params_glob%box)
        do i=batchlims(1),batchlims(2)
            ii = i - batchlims(1) + 1
            if( uuse_denoised )then
                call build_glob%spproj%get_stkname_and_ind(params_glob%oritype, pinds(i), stkname, ind_in_stk, stkname_den)
                call dstkio_r%read(stkname_den, ind_in_stk, build_glob%imgbatch(ii))
            else
                call build_glob%spproj%get_stkname_and_ind(params_glob%oritype, pinds(i), stkname, ind_in_stk)
                call dstkio_r%read(stkname, ind_in_stk, build_glob%imgbatch(ii))
            endif
        end do
        call dstkio_r%kill
    end subroutine discrete_read_imgbatch

    subroutine set_bp_range( cline )
        class(cmdline), intent(in) :: cline
        real, allocatable     :: resarr(:), fsc_arr(:)
        real                  :: fsc0143, fsc05
        real                  :: mapres(params_glob%nstates)
        integer               :: s, loc(1), lp_ind, arr_sz, fsc_sz
        character(len=STDLEN) :: fsc_fname
        logical               :: fsc_bin_exists(params_glob%nstates), all_fsc_bin_exist
        if( params_glob%l_lpset )then
            ! set Fourier index range
            params_glob%kfromto(2) = calc_fourier_index(params_glob%lp, params_glob%box, params_glob%smpd)
            if( cline%defined('lpstop') )then
                params_glob%kfromto(2) = min(params_glob%kfromto(2),&
                    &calc_fourier_index(params_glob%lpstop, params_glob%box, params_glob%smpd))
            endif
            ! FSC values are read anyway
            do s=1,params_glob%nstates
                fsc_fname = trim(FSC_FBODY)//int2str_pad(s,2)//BIN_EXT
                if( file_exists(fsc_fname) )then
                    fsc_arr = file2rarr(trim(adjustl(fsc_fname)))
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
                fsc_fname = trim(FSC_FBODY)//int2str_pad(s,2)//BIN_EXT
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
                        fsc_fname = trim(FSC_FBODY)//int2str_pad(s,2)//BIN_EXT
                        fsc_arr = file2rarr(trim(adjustl(fsc_fname)))
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
                loc = maxloc(mapres) ! worst resolved
                if( params_glob%nstates == 1 )then
                    ! get median updatecnt
                    if( build_glob%spproj_field%median('updatecnt') > 1.0 )then ! more than half have been updated
                        lp_ind = get_find_at_corr(build_glob%fsc(1,:), params_glob%lplim_crit, incrreslim=params_glob%l_incrreslim)
                    else
                        lp_ind = get_find_at_corr(build_glob%fsc(1,:), 0.5, incrreslim=params_glob%l_incrreslim) ! more conservative limit @ start
                    endif
                else
                    lp_ind = get_find_at_corr(build_glob%fsc(loc(1),:), params_glob%lplim_crit)
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
        ! update low-pas limit in project
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

    subroutine sample_ptcls4update( pfromto, l_incr_sampl, nptcls2update, pinds, ptcl_mask )
        integer,              intent(in)    :: pfromto(2)
        logical,              intent(in)    :: l_incr_sampl
        integer,              intent(inout) :: nptcls2update
        integer, allocatable, intent(inout) :: pinds(:)
        logical, allocatable, intent(inout) :: ptcl_mask(:)
        type(class_sample),   allocatable   :: clssmp(:)
        integer :: icls
        if( params_glob%l_update_frac )then
            if( trim(params_glob%balance).eq.'yes' )then
                if( file_exists(CLASS_SAMPLING_FILE) )then
                    call read_class_samples(clssmp, CLASS_SAMPLING_FILE)
                else
                    THROW_HARD('File for class-biased sampling in fractional update: '//CLASS_SAMPLING_FILE//' does not exists!')
                endif
                ! balanced class sampling 
                call build_glob%spproj_field%sample4update_class(clssmp, params_glob%greediness,&
                &pfromto, params_glob%update_frac, nptcls2update, pinds, ptcl_mask, l_incr_sampl)
                call deallocate_class_samples(clssmp)
            else
                call build_glob%spproj_field%sample4update_rnd(pfromto,&
                &params_glob%update_frac, nptcls2update, pinds, ptcl_mask, l_incr_sampl)
            endif
        else
            ! we sample all state > 0
            call build_glob%spproj_field%sample4update_all(pfromto, nptcls2update, pinds, ptcl_mask, l_incr_sampl)
        endif
        if( l_incr_sampl )then
            ! increment update counter
            call build_glob%spproj_field%incr_updatecnt(pfromto, ptcl_mask)
        endif
    end subroutine sample_ptcls4update

    !>  \brief  prepares one particle image for alignment
    !!          serial routine
    subroutine prepimg4align( iptcl, img, img_out )
        use simple_ctf, only: ctf
        integer,      intent(in)    :: iptcl
        class(image), intent(inout) :: img
        class(image), intent(inout) :: img_out
        type(ctf)       :: tfun
        type(ctfparams) :: ctfparms
        real            :: x, y, sdev_noise, crop_factor
        ! Normalise
        call img%norm_noise(build_glob%lmsk, sdev_noise)
        ! Fourier cropping
        call img%fft()
        call img%clip(img_out)
        ! Shift image to rotational origin
        crop_factor = real(params_glob%box_crop) / real(params_glob%box)
        x = build_glob%spproj_field%get(iptcl, 'x') * crop_factor
        y = build_glob%spproj_field%get(iptcl, 'y') * crop_factor
        if(abs(x) > SHTHRESH .or. abs(y) > SHTHRESH)then
            call img_out%shift2Dserial([-x,-y])
        endif
        ! Phase-flipping
        ctfparms = build_glob%spproj%get_ctfparams(params_glob%oritype, iptcl)
        select case(ctfparms%ctfflag)
            case(CTFFLAG_NO, CTFFLAG_FLIP)
                ! nothing to do
            case(CTFFLAG_YES)
                ctfparms%smpd = ctfparms%smpd / crop_factor != smpd_crop
                tfun = ctf(ctfparms%smpd, ctfparms%kv, ctfparms%cs, ctfparms%fraca)
                call tfun%apply_serial(img_out, 'flip', ctfparms)
            case DEFAULT
                THROW_HARD('unsupported CTF flag: '//int2str(ctfparms%ctfflag)//' prepimg4align')
        end select
        ! Back to real space
        call img_out%ifft
        ! Soft-edged mask
        if( params_glob%l_focusmsk )then
            call img_out%mask(params_glob%focusmsk*crop_factor, 'soft')
        else
            if( params_glob%l_needs_sigma )then
                call img_out%mask(params_glob%msk_crop, 'softavg')
            else
                call img_out%mask(params_glob%msk_crop, 'soft')
            endif
        endif
        ! gridding prep
        if( params_glob%gridding.eq.'yes' ) call build_glob%img_crop_polarizer%div_by_instrfun(img_out)
        ! return to Fourier space
        call img_out%fft()
    end subroutine prepimg4align

    !>  \brief  prepares one cluster centre image for alignment
    subroutine prep2Dref( img_in, img_out, icls, iseven, center, xyz_in, xyz_out )
        class(image),      intent(inout) :: img_in
        class(image),      intent(inout) :: img_out
        integer,           intent(in)    :: icls
        logical,           intent(in)    :: iseven
        logical, optional, intent(in)    :: center
        real,    optional, intent(in)    :: xyz_in(3)
        real,    optional, intent(out)   :: xyz_out(3)
        integer :: filtsz
        real    :: frc(img_out%get_filtsz()), filter(img_out%get_filtsz())
        real    :: xy_cavg(2), xyz(3), sharg, crop_factor
        logical :: do_center
        filtsz = img_in%get_filtsz()
        crop_factor = real(params_glob%box_crop) / real(params_glob%box)
        ! centering only performed if params_glob%center.eq.'yes'
        do_center = (params_glob%center .eq. 'yes')
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
                if( trim(params_glob%masscen).ne.'yes' )then
                    call build_glob%spproj_field%calc_avg_offset2D(icls, xy_cavg)
                    if( arg(xy_cavg) < CENTHRESH )then
                        xyz = 0.
                    else if( arg(xy_cavg) > MAXCENTHRESH2D )then
                        xyz(1:2) = xy_cavg * crop_factor
                        xyz(3)   = 0.
                    else
                        xyz = img_in%calc_shiftcen_serial(params_glob%cenlp, params_glob%msk_crop, iter_center=(params_glob%iter_center .eq. 'yes'))
                        if( arg(xyz(1:2)/crop_factor - xy_cavg) > MAXCENTHRESH2D ) xyz = 0.
                    endif
                else
                    xyz = img_in%calc_shiftcen_serial(params_glob%cenlp, params_glob%msk_crop, iter_center=(params_glob%iter_center .eq. 'yes'))
                endif
                sharg = arg(xyz)
                if( sharg > CENTHRESH )then
                    ! apply shift and update the corresponding class parameters
                    call img_in%fft()
                    call img_in%shift2Dserial(xyz(1:2))
                    call build_glob%spproj_field%add_shift2class(icls, -xyz(1:2) / crop_factor)
                else
                    xyz = 0.
                endif
                if( present(xyz_out) ) xyz_out = xyz
            endif
        endif
        if( params_glob%l_ml_reg )then
            ! no filtering
        else
            if( params_glob%l_lpset.and.params_glob%l_icm )then
                ! ICM filter only applied when lp is set and performed below, FRC filtering turned off
            else
                ! FRC-based filtering
                call build_glob%clsfrcs%frc_getter(icls, params_glob%hpind_fsc, params_glob%l_phaseplate, frc)
                if( any(frc > 0.143) )then
                    call fsc2optlp_sub(filtsz, frc, filter, merged=params_glob%l_lpset)
                    call img_in%fft() ! needs to be here in case the shift was never applied (above)
                    call img_in%apply_filter_serial(filter)
                endif
            endif
        endif
        ! ensure we are in real-space
        call img_in%ifft()
        ! clip image if needed
        call img_in%clip(img_out)
        ! ICM filter
        if( params_glob%l_lpset.and.params_glob%l_icm )then
            call img_out%ICM2D( params_glob%lambda, verbose=.false. )
        endif
        ! apply mask
        call img_out%mask(params_glob%msk_crop, 'soft', backgr=0.0)
        ! gridding prep
        if( params_glob%gridding.eq.'yes' ) call build_glob%img_crop_polarizer%div_by_instrfun(img_out)
        ! move to Fourier space
        call img_out%fft()
    end subroutine prep2Dref

    !>  \brief  initializes all volumes for reconstruction
    subroutine preprecvols
        character(len=:), allocatable :: part_str, fbody
        integer,          allocatable :: pops(:)
        integer :: istate
        allocate(part_str, source=int2str_pad(params_glob%part,params_glob%numlen))
        call build_glob%spproj_field%get_pops(pops, 'state')
        do istate = 1, params_glob%nstates
            if( pops(istate) > 0)then
                call build_glob%eorecvols(istate)%new(build_glob%spproj)
                call build_glob%eorecvols(istate)%reset_all
                if( params_glob%l_trail_rec )then
                    fbody = trim(VOL_FBODY)//int2str_pad(istate,2)//'_part'//part_str
                    if( build_glob%eorecvols(istate)%ldim_even_match(fbody) )then
                        call build_glob%eorecvols(istate)%read_eos(fbody)
                        call build_glob%eorecvols(istate)%apply_weight(1.-params_glob%update_frac)
                        call build_glob%eorecvols(istate)%expand_exp
                    endif
                endif
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
        class(cmdline),   intent(in)  :: cline
        integer,          intent(in)  :: s
        character(len=*), intent(in)  :: volfname
        logical,          intent(out) :: do_center
        real,             intent(out) :: xyz(3)
        logical,          intent(in)  :: map_shift
        real    :: crop_factor
        logical :: has_been_searched
        do_center   = .true.
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
        xyz = build_glob%vol%calc_shiftcen(params_glob%cenlp,params_glob%msk_crop, iter_center=(params_glob%iter_center .eq. 'yes'))
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

    subroutine estimate_lp_refvols( s, lpopt )
        use simple_opt_filter, only: estimate_lplim
        integer, intent(in)  :: s
        real,    intent(out) :: lpopt
        type(image) :: mskvol
        integer     :: npix
        logical     :: l_update_lp
        call mskvol%disc([params_glob%box_crop,params_glob%box_crop,params_glob%box_crop],&
                            &params_glob%smpd_crop, params_glob%msk_crop, npix )
        if( params_glob%lp_auto.eq.'fsc' )then
            lpopt = calc_lowpass_lim(get_find_at_corr(build_glob%fsc(s,:), params_glob%lplim_crit),&
            &params_glob%box_crop, params_glob%smpd_crop)
        else
            call estimate_lplim(build_glob%vol_odd, build_glob%vol, mskvol, [params_glob%lpstart,params_glob%lpstop], lpopt)
        endif
        l_update_lp = .false.
        if( s == 1 )then
            l_update_lp = .true. ! always update for state == 1       
        else
            if( lpopt > params_glob%lp )then
                l_update_lp = .true. ! the limit for the state with the worst resolution wins
            endif
        endif
        if( l_update_lp )then
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

    subroutine read_and_filter_refvols( s )
        integer, intent(in) :: s
        character(len=:), allocatable :: vol_even, vol_even_unfil, vol_odd, vol_odd_unfil, vol_avg
        real    :: cur_fil(params_glob%box_crop), lpopt, lpest
        integer :: filtsz
        vol_even = params_glob%vols_even(s)
        vol_odd  = params_glob%vols_odd(s)
        vol_avg  = params_glob%vols(s)
        if( params_glob%l_lpauto .and. params_glob%l_ml_reg )then
            vol_even_unfil = add2fbody(vol_even,params_glob%ext,'_unfil')
            vol_odd_unfil  = add2fbody(vol_odd,params_glob%ext,'_unfil')
            ! estimate low-pass limit from unifiltered volumes
            call build_glob%vol%read_and_crop(vol_even_unfil,    params_glob%smpd, params_glob%box_crop, params_glob%smpd_crop)
            call build_glob%vol_odd%read_and_crop(vol_odd_unfil, params_glob%smpd, params_glob%box_crop, params_glob%smpd_crop)
            call estimate_lp_refvols(s, lpopt)
            ! read in filtered volumes (ML-reg)
            call build_glob%vol%read_and_crop(vol_even,          params_glob%smpd, params_glob%box_crop, params_glob%smpd_crop)
            call build_glob%vol_odd%read_and_crop(vol_odd,       params_glob%smpd, params_glob%box_crop, params_glob%smpd_crop)
        else
            call build_glob%vol%read_and_crop(vol_even,          params_glob%smpd, params_glob%box_crop, params_glob%smpd_crop)
            call build_glob%vol_odd%read_and_crop(vol_odd,       params_glob%smpd, params_glob%box_crop, params_glob%smpd_crop)
            if( params_glob%l_lpauto ) call estimate_lp_refvols(s, lpopt)
        endif
        if( params_glob%l_icm )then
            call build_glob%vol%ICM3D_eo(build_glob%vol_odd, params_glob%lambda)
            if( params_glob%l_lpset )then ! no independent volume registration, so average eo pairs
                call build_glob%vol%add(build_glob%vol_odd)
                call build_glob%vol%mul(0.5)
                call build_glob%vol_odd%copy(build_glob%vol)
            endif
        else if( params_glob%l_lpset )then
            ! the average volume occupies both even and odd
            call build_glob%vol%read_and_crop(vol_avg, params_glob%smpd, params_glob%box_crop, params_glob%smpd_crop)
            call build_glob%vol_odd%copy(build_glob%vol)
        endif
        call build_glob%vol%fft
        call build_glob%vol_odd%fft
        lpest = params_glob%lp
        if( params_glob%l_lpauto ) lpest = lpopt
        if( params_glob%l_ml_reg )then
            ! filtering done when volumes are assembled
        else if( params_glob%l_icm )then
            ! filtering done above
        else if( params_glob%l_lpset )then
            ! Cosine low-pass filter, works best for nanoparticles
            call build_glob%vol%bp(0., lpest)
        else
            filtsz = build_glob%vol%get_filtsz()
            if( any(build_glob%fsc(s,:) > 0.143) )then
                call fsc2optlp_sub(filtsz,build_glob%fsc(s,:),cur_fil)
                call build_glob%vol%apply_filter(cur_fil)
            endif
        endif
    end subroutine read_and_filter_refvols

    !>  \brief  prepares one volume for references extraction
    subroutine preprefvol( cline, s, do_center, xyz, iseven )
        use simple_projector,          only: projector
        use simple_butterworth,        only: butterworth_filter
        use simple_nanoparticle_utils, only: phasecorr_one_atom
        class(cmdline), intent(in) :: cline
        integer,        intent(in) :: s
        logical,        intent(in) :: do_center
        real,           intent(in) :: xyz(3)
        logical,        intent(in) :: iseven
        type(projector), pointer :: vol_ptr => null()
        type(image)              :: mskvol
        if( iseven )then
            vol_ptr => build_glob%vol
        else
            vol_ptr => build_glob%vol_odd
        endif
        if( do_center )then
            call vol_ptr%fft()
            call vol_ptr%shift([xyz(1),xyz(2),xyz(3)])
        endif
        ! back to real space
        call vol_ptr%ifft()
        ! noise regularization
        if( params_glob%l_noise_reg )then
            call vol_ptr%add_gauran(params_glob%eps)
        endif
        ! masking
        if( params_glob%l_filemsk )then
            ! envelope masking
            call mskvol%new([params_glob%box_crop,params_glob%box_crop,params_glob%box_crop],params_glob%smpd_crop)
            call mskvol%read(params_glob%mskfile)
            call vol_ptr%zero_env_background(mskvol)
            call vol_ptr%mul(mskvol)
            call mskvol%kill
        else
            ! circular masking
            call vol_ptr%mask(params_glob%msk_crop, 'soft', backgr=0.0)
        endif
        ! gridding prep
        if( params_glob%gridding.eq.'yes' )then
            call vol_ptr%div_w_instrfun(params_glob%interpfun, alpha=params_glob%alpha)
        endif
        ! FT volume
        call vol_ptr%fft()
        ! expand for fast interpolation & correct for norm when clipped
        call vol_ptr%expand_cmat(params_glob%alpha,norm4proj=.true.)
    end subroutine preprefvol

    !>  \brief  grids one particle image to the volume
    subroutine grid_ptcl( fpl, se, o )
        use simple_fplane,      only   : fplane
        class(fplane),   intent(in)    :: fpl
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
        if( pw > TINY ) call build_glob%eorecvols(s)%grid_plane(se, o, fpl, eo, pw)
    end subroutine grid_ptcl

    !> volumetric 3d reconstruction
    subroutine calc_3Drec( cline, nptcls2update, pinds )
        use simple_fplane, only: fplane
        class(cmdline),    intent(inout) :: cline
        integer,           intent(in)    :: nptcls2update
        integer,           intent(in)    :: pinds(nptcls2update)
        type(fplane),      allocatable   :: fpls(:)
        type(ctfparams),   allocatable   :: ctfparms(:)
        type(ori)        :: orientation
        real             :: shift(2), sdev_noise
        integer          :: batchlims(2), iptcl, i, i_batch, ibatch
        ! init volumes
        call preprecvols
        ! prep batch imgs
        call prepimgbatch(MAXIMGBATCHSZ)
        ! allocate array
        allocate(fpls(MAXIMGBATCHSZ),ctfparms(MAXIMGBATCHSZ))
        ! gridding batch loop
        do i_batch=1,nptcls2update,MAXIMGBATCHSZ
            batchlims = [i_batch,min(nptcls2update,i_batch + MAXIMGBATCHSZ - 1)]
            call discrete_read_imgbatch( nptcls2update, pinds, batchlims)
            !$omp parallel do default(shared) private(i,iptcl,ibatch,shift,sdev_noise) schedule(static) proc_bind(close)
            do i=batchlims(1),batchlims(2)
                iptcl  = pinds(i)
                ibatch = i - batchlims(1) + 1
                if( .not.fpls(ibatch)%does_exist() ) call fpls(ibatch)%new(build_glob%imgbatch(1))
                call build_glob%imgbatch(ibatch)%norm_noise(build_glob%lmsk, sdev_noise)
                call build_glob%imgbatch(ibatch)%fft
                ctfparms(ibatch) = build_glob%spproj%get_ctfparams(params_glob%oritype, iptcl)
                shift = build_glob%spproj_field%get_2Dshift(iptcl)
                call fpls(ibatch)%gen_planes(build_glob%imgbatch(ibatch), ctfparms(ibatch), shift, iptcl)
            end do
            !$omp end parallel do
            ! gridding
            do i=batchlims(1),batchlims(2)
                iptcl       = pinds(i)
                ibatch      = i - batchlims(1) + 1
                call build_glob%spproj_field%get_ori(iptcl, orientation)
                if( orientation%isstatezero() ) cycle
                call grid_ptcl(fpls(ibatch), build_glob%pgrpsyms, orientation)
            end do
        end do
        ! normalise structure factors
        call norm_struct_facts( cline )
        ! destruct
        call killrecvols()
        do ibatch=1,MAXIMGBATCHSZ
            call fpls(ibatch)%kill
        end do
        deallocate(fpls,ctfparms)
        call orientation%kill
    end subroutine calc_3Drec

    !> Volumetric 3d reconstruction from summed projection directions
    !> Only supports single state
    subroutine calc_projdir3Drec( cline, nptcls2update, pinds )
        !$ use omp_lib
        !$ use omp_lib_kinds
        use simple_fplane,   only: fplane
        use simple_gridding, only: gen_instrfun_img
        use simple_timer
        class(cmdline),    intent(inout) :: cline
        integer,           intent(in)    :: nptcls2update
        integer,           intent(in)    :: pinds(nptcls2update)
        type(fplane),      allocatable   :: fpls(:), projdirs(:,:)
        integer,           allocatable   :: eopops(:,:), tmp(:,:)
        type(ctfparams) :: ctfparms
        type(image)     :: instrimg, numimg, denomimg
        type(ori)       :: orientation
        real            :: shift(2), e3, sdev_noise, w
        integer         :: batchlims(2), iptcl, i,j, i_batch, ibatch, iproj, eo, peo, ithr, pproj
        logical         :: DEBUG    = .false.
        logical         :: BILINEAR = .true.
        integer(timer_int_kind) :: t
        real(timer_int_kind)    :: t_ini, t_pad, t_sum, t_rec
        if( params_glob%nstates /= 1   ) THROW_HARD('PROJREC & NSTATES>1 not supported yet')
        if( DEBUG ) t = tic()
        ! init volumes
        call preprecvols
        ! prep batch imgs
        call prepimgbatch(MAXIMGBATCHSZ)
        ! allocations
        allocate(fpls(MAXIMGBATCHSZ),projdirs(params_glob%nspace,2),eopops(params_glob%nspace,2))
        ! e/o projection directions populations
        eopops = 0
        !$omp parallel default(shared) private(i,iptcl,iproj,eo,ibatch) proc_bind(close)
        !$omp do schedule(static) reduction(+:eopops)
        do i = 1,nptcls2update
            iptcl  = pinds(i)
            iproj  = nint(build_glob%spproj_field%get(iptcl, 'proj'))
            eo     = build_glob%spproj_field%get_eo(iptcl)+1
            eopops(iproj,eo) = eopops(iproj,eo) + 1
        end do
        !$omp end do
        ! projection direction slices to insert into volume
        !$omp do schedule(static)
        do iproj = 1,params_glob%nspace
            if( eopops(iproj,1) > 0 ) call projdirs(iproj,1)%new(build_glob%imgbatch(1),pad=.true., genplane=.false.)
            if( eopops(iproj,2) > 0 ) call projdirs(iproj,2)%new(build_glob%imgbatch(1),pad=.true., genplane=.false.)
        end do
        !$omp end do nowait
        ! particles to be summed into projection direction slices
        !$omp do schedule(static)
        do ibatch = 1,min(nptcls2update,MAXIMGBATCHSZ)
            call fpls(ibatch)%new(build_glob%imgbatch(1),pad=.true.)
        end do
        !$omp end do
        !$omp end parallel
        ! instrument function
        call instrimg%new([params_glob%box,params_glob%box,1], params_glob%smpd)
        if( BILINEAR )then
            call gen_instrfun_img(instrimg, 'linear', padded_dim=params_glob%boxpd, norm=.true.)
        else
            call gen_instrfun_img(instrimg, 'kb', fpls(1)%kbwin, padded_dim=params_glob%boxpd, norm=.true.)
        endif
        if( DEBUG ) t_ini = toc(t)
        ! batch loop
        tmp = eopops
        if( DEBUG )then
            t_pad = 0.
            t_sum = 0.
        endif
        do i_batch=1,nptcls2update,MAXIMGBATCHSZ
            batchlims = [i_batch,min(nptcls2update,i_batch + MAXIMGBATCHSZ - 1)]
            ! particles in-plane transformation
            call discrete_read_imgbatch( nptcls2update, pinds, batchlims)
            if( DEBUG ) t = tic()
            !$omp parallel do default(shared) private(i,iptcl,ibatch,shift,e3,sdev_noise,ctfparms,ithr)&
            !$omp schedule(static) proc_bind(close)
            do i = batchlims(1),batchlims(2)
                iptcl    = pinds(i)
                ibatch   = i - batchlims(1) + 1
                ctfparms = build_glob%spproj%get_ctfparams(params_glob%oritype, iptcl)
                shift    = build_glob%spproj_field%get_2Dshift(iptcl)
                e3       = build_glob%spproj_field%e3get(iptcl)
                call build_glob%imgbatch(ibatch)%norm_noise(build_glob%lmsk, sdev_noise)
                call build_glob%imgbatch(ibatch)%div(instrimg)
                call build_glob%imgbatch(ibatch)%fft
                call fpls(ibatch)%gen_planes_pad(build_glob%imgbatch(ibatch), ctfparms, shift, e3, iptcl, BILINEAR)
            end do
            !$omp end parallel do
            if( DEBUG )then
                t_pad = t_pad + toc(t)
                t = tic()
            endif
            ! particles summation
            ! FAST
            !$omp parallel do default(shared) private(i,j,iproj,pproj,iptcl,ibatch,w,eo,peo)&
            !$omp schedule(dynamic) proc_bind(close)
            do j = 1,2*params_glob%nspace
                ! For better e/o balancing
                if( j <= params_glob%nspace )then
                    iproj = j
                    eo    = 1
                else
                    iproj = j - params_glob%nspace
                    eo    = 2
                endif
                if( tmp(iproj,eo)==0 ) cycle
                do i = batchlims(1),batchlims(2)
                    iptcl  = pinds(i)
                    pproj  = nint(build_glob%spproj_field%get(iptcl, 'proj'))
                    if( iproj /= pproj ) cycle
                    peo = build_glob%spproj_field%get_eo(iptcl)+1
                    if( peo /= eo ) cycle
                    w = build_glob%spproj_field%get(iptcl, 'w')
                    if( w < TINY ) cycle
                    ibatch = i - batchlims(1) + 1
                    projdirs(iproj,eo)%cmplx_plane = projdirs(iproj,eo)%cmplx_plane + w * fpls(ibatch)%cmplx_plane
                    projdirs(iproj,eo)%ctfsq_plane = projdirs(iproj,eo)%ctfsq_plane + w * fpls(ibatch)%ctfsq_plane
                    tmp(iproj,eo) = tmp(iproj,eo) - 1
                enddo
            enddo
            !$omp end parallel do
            if( DEBUG ) t_sum = t_sum + toc(t)
        end do
        ! some cleanup
        !$omp parallel do default(shared) private(ibatch) schedule(static) proc_bind(close)
        do ibatch = 1,size(fpls)
            call fpls(ibatch)%kill
        end do
        !$omp end parallel do
        deallocate(fpls,tmp)
        ! projections directions reconstructon
        if( DEBUG ) t = tic()
        do iproj = 1,params_glob%nspace
            call build_glob%eulspace%get_ori(iproj, orientation)
            call orientation%set('state',1.)
            call orientation%set('w',    1.)
            if( eopops(iproj,1) > 0 )then
                call orientation%set('eo',0.)
                call grid_ptcl(projdirs(iproj,1), build_glob%pgrpsyms, orientation)
            endif
            if( eopops(iproj,2) > 0 )then
                call orientation%set('eo',1.)
                call grid_ptcl(projdirs(iproj,2), build_glob%pgrpsyms, orientation)
            endif
        end do
        if( DEBUG ) t_rec = toc(t)
        if( DEBUG )then
            print *,'timing: ',t_ini, t_pad, t_sum, t_rec
            call build_glob%eulspace%write('reforis.txt')
            call build_glob%spproj_field%write('oris.txt')
            call instrimg%write('instrumentfunction.mrc')
            do iproj = 1,params_glob%nspace
                if( eopops(iproj,1) > 0 )then
                    call projdirs(iproj,1)%convert2img(numimg, denomimg)
                    ! call numimg%ctf_dens_correct(denomimg)
                    ! call numimg%ifft
                    call denomimg%write('projdirs_even.mrc',iproj)
                endif
                if( eopops(iproj,2) > 0 )then
                    call projdirs(iproj,2)%convert2img(numimg, denomimg)
                    ! call numimg%ctf_dens_correct(denomimg)
                    ! call numimg%ifft
                    call denomimg%write('projdirs_odd.mrc',iproj)
                endif
            end do
            call numimg%kill
            call denomimg%kill
        endif
        ! some cleanup
        !$omp parallel do default(shared) private(iproj) schedule(static) proc_bind(close)
        do iproj = 1,params_glob%nspace
            call projdirs(iproj,1)%kill
            call projdirs(iproj,2)%kill
        end do
        !$omp end parallel do
        deallocate(projdirs)
        ! normalise structure factors
        call norm_struct_facts( cline )
        ! cleanup
        call instrimg%kill
        call killrecvols()
        call orientation%kill
    end subroutine calc_projdir3Drec

    subroutine norm_struct_facts( cline )
        use simple_masker, only: masker
        class(cmdline),    intent(inout) :: cline
        character(len=:), allocatable    :: iter_str
        character(len=STDLEN) :: pprocvol, lpvol
        real, allocatable     :: optlp(:), res(:)
        type(masker)          :: envmsk
        integer               :: s, find4eoavg, ldim(3)
        real                  :: res05s(params_glob%nstates), res0143s(params_glob%nstates), lplim, bfac
        ! init
        ldim = [params_glob%box_crop,params_glob%box_crop,params_glob%box_crop]
        call build_glob%vol%new(ldim,params_glob%smpd_crop)
        call build_glob%vol2%new(ldim,params_glob%smpd_crop)
        res0143s = 0.
        res05s   = 0.
        ! cycle through states
        do s=1,params_glob%nstates
            if( build_glob%spproj_field%get_pop(s, 'state') == 0 )then
                ! empty state
                build_glob%fsc(s,:) = 0.
                cycle
            endif
            call build_glob%eorecvols(s)%compress_exp
            if( params_glob%l_distr_exec )then
                call build_glob%eorecvols(s)%write_eos(VOL_FBODY//int2str_pad(s,2)//'_part'//&
                    int2str_pad(params_glob%part,params_glob%numlen))
            else
                ! global volume name update
                params_glob%vols(s) = VOL_FBODY//int2str_pad(s,2)//params_glob%ext
                ! updating command-line according (needed in multi-stage wflows)
                call cline%set('vol'//int2str(s), params_glob%vols(s))
                if( params_glob%l_filemsk .and. params_glob%l_envfsc )then
                    call build_glob%eorecvols(s)%set_automsk(.true.)
                endif
                params_glob%vols_even(s) = add2fbody(params_glob%vols(s), params_glob%ext, '_even')
                params_glob%vols_odd(s)  = add2fbody(params_glob%vols(s), params_glob%ext, '_odd')
                if( params_glob%l_ml_reg )then
                    call build_glob%eorecvols(s)%sampl_dens_correct_eos(s, params_glob%vols_even(s), &
                        &params_glob%vols_odd(s), find4eoavg)
                    call build_glob%eorecvols(s)%get_res(res05s(s), res0143s(s))
                    call build_glob%eorecvols(s)%sum_eos
                else
                    call build_glob%eorecvols(s)%sum_eos
                    call build_glob%eorecvols(s)%sampl_dens_correct_eos(s, params_glob%vols_even(s), &
                        &params_glob%vols_odd(s), find4eoavg)
                    call build_glob%eorecvols(s)%get_res(res05s(s), res0143s(s))
                endif
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
                ! post-process volume
                pprocvol = add2fbody(trim(params_glob%vols(s)), params_glob%ext, PPROC_SUFFIX)
                lpvol    = add2fbody(trim(params_glob%vols(s)), params_glob%ext, LP_SUFFIX)
                build_glob%fsc(s,:) = file2rarr('fsc_state'//int2str_pad(s,2)//'.bin')
                ! low-pass limit
                if( params_glob%l_lpset )then
                    lplim = params_glob%lp
                else
                    lplim = res0143s(s)
                endif
                ! B-factor estimation
                if( cline%defined('bfac') )then
                    bfac = params_glob%bfac
                else
                    if( lplim < 5. )then
                        bfac = build_glob%vol%guinier_bfac(HPLIM_GUINIER, lplim)
                        write(logfhandle,'(A,1X,F8.2)') '>>> B-FACTOR DETERMINED TO:', bfac
                    else
                        bfac = 0.
                    endif
                endif
                ! B-factor application
                call build_glob%vol2%copy(build_glob%vol)
                call build_glob%vol%apply_bfac(bfac)
                ! low-pass filter
                if( params_glob%l_lpset )then
                    call build_glob%vol%bp(0., lplim)
                    call build_glob%vol2%bp(0., lplim)
                else
                    res   = build_glob%vol%get_res()
                    optlp = fsc2optlp(build_glob%fsc(s,:))
                    where( res < TINY ) optlp = 0.
                    lplim = res0143s(s)
                    ! optimal low-pass filter from FSC
                    call build_glob%vol%apply_filter(optlp)
                    call build_glob%vol2%apply_filter(optlp)
                    ! final low-pass filtering for smoothness
                    call build_glob%vol%bp(0., res0143s(s))
                    call build_glob%vol2%bp(0., res0143s(s))
                endif
                call build_glob%vol%ifft()
                call build_glob%vol2%ifft()
                ! write low-pass filtered without B-factor or mask
                call build_glob%vol2%write(lpvol)
                ! masking
                if( params_glob%l_filemsk )then
                    if( params_glob%l_filemsk )then
                        call envmsk%new(ldim, params_glob%smpd_crop)
                        call envmsk%read(params_glob%mskfile)
                    endif
                    call build_glob%vol%zero_background
                    if( cline%defined('lp_backgr') )then
                        call build_glob%vol%lp_background(envmsk,params_glob%lp_backgr)
                    else
                        call build_glob%vol%mul(envmsk)
                    endif
                    call envmsk%kill
                else
                    call build_glob%vol%mask(params_glob%msk_crop, 'soft')
                endif
                ! write
                call build_glob%vol%write(pprocvol)
            endif
        end do
        if((.not.params_glob%l_distr_exec) .and. (.not.params_glob%l_lpset))then
            ! set the resolution limit according to the worst resolved model
            params_glob%lp = min(params_glob%lp,max(params_glob%lpstop,maxval(res0143s)))
        endif
        call build_glob%vol2%kill
    end subroutine norm_struct_facts
    
end module simple_strategy2D3D_common
