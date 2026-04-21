!@descr: search-space and particle-selection policy routines for matcher workflows
module simple_matcher_smpl_and_lplims
use simple_pftc_srch_api
use simple_builder, only: builder
implicit none

public :: set_bp_range3D, set_bp_range2D
public :: sample_ptcls4update3D, sample_ptcls4fillin, sample_ptcls4update2D
private
#include "simple_local_flags.inc"

contains

    subroutine set_bp_range3D( params, build, cline )
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
                        mapres(s)      = 0.
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
                THROW_HARD('no method for low-pass limit: need fsc file or lp field; set_bp_range3D')
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
    end subroutine set_bp_range3D

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

    subroutine sample_ptcls4update3D( params, build, pfromto, l_incr_sampl, nptcls2update, pinds )
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
    end subroutine sample_ptcls4update3D

    subroutine sample_ptcls4fillin( params, build, pfromto, l_incr_sampl, nptcls2update, pinds )
        class(parameters),    intent(in)    :: params
        class(builder),       intent(inout) :: build
        integer,              intent(in)    :: pfromto(2)
        logical,              intent(in)    :: l_incr_sampl
        integer,              intent(inout) :: nptcls2update
        integer, allocatable, intent(inout) :: pinds(:)
        call build%spproj_field%sample4update_fillin(pfromto, params%update_frac, nptcls2update, pinds, l_incr_sampl)
    end subroutine sample_ptcls4fillin

    subroutine sample_ptcls4update2D( params, build, pfromto, l_updatefrac, nptcls, pinds )
        class(parameters),    intent(in)    :: params
        class(builder),       intent(inout) :: build
        logical,              intent(in)    :: l_updatefrac
        integer,              intent(in)    :: pfromto(2)
        integer,              intent(inout) :: nptcls
        integer, allocatable, intent(inout) :: pinds(:)
        logical :: l_has_been_sampled
        l_has_been_sampled = build%spproj_field%has_been_sampled()
        if( l_updatefrac )then
            ! abinitio2D controller policy:
            ! startit==1  -> sticky sampled subset stage
            ! startit>1   -> stochastic resampling biased toward low updatecnt
            write(logfhandle,'(A,2I5,L2,F8.3,L2,A,2I8)') &
                '>>> ABINITIO2D SAMPLING DEBUG: startit/which_iter=', params%startit, params%which_iter, &
                l_updatefrac, params%update_frac, l_has_been_sampled, ' pfromto=', pfromto(1), pfromto(2)
            if( params%startit == 1 )then
                if( l_has_been_sampled )then
                    write(logfhandle,'(A)') '>>> ABINITIO2D SAMPLING DEBUG: using sample4update_reprod'
                    call build%spproj_field%sample4update_reprod(pfromto, nptcls, pinds)
                else
                    write(logfhandle,'(A)') '>>> ABINITIO2D SAMPLING DEBUG: using sample4update_rnd'
                    call build%spproj_field%sample4update_rnd(pfromto, params%update_frac, nptcls, pinds, .true.)
                    call build%spproj_field%set_updatecnt(1, pinds)
                endif
            else
                write(logfhandle,'(A)') '>>> ABINITIO2D SAMPLING DEBUG: using sample4update_cnt'
                call build%spproj_field%sample4update_cnt(pfromto, params%update_frac, nptcls, pinds, .true.)
            endif
        else
            write(logfhandle,'(A,2I5,L2,F8.3,A,2I8)') &
                '>>> ABINITIO2D SAMPLING DEBUG: startit/which_iter=', params%startit, params%which_iter, &
                l_updatefrac, params%update_frac, ' using sample4update_all pfromto=', pfromto(1), pfromto(2)
            call build%spproj_field%sample4update_all(pfromto, nptcls, pinds, .true.)
        endif
    end subroutine sample_ptcls4update2D

end module simple_matcher_smpl_and_lplims
